/*

"yet another smoothed particle hydrodynamics" (YASPH) simulation code

*/
 
// TODO: cut out glue I/O code into a separate header file

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>

#include "textio.h"
#include "argsio.h"
#include "kernels2d.h"
#include "hashindex2d.h"
#include "autotune.h"
#include "clkutil.h"

typedef struct tParticle {
  double m;
  double x;     // dynamic state is (x, y, vx, vy, u)
  double y;
  double vx;
  double vy;
  double u;
  double rho;   // evaluated from state (m, x, y); exact conservation of mass
  double p;     // used to compute accelerations
  double c;     // local speed of sound
  double mumax;
  double vxdot;
  double vydot;
  double udot;
  //int type;
  //int status;
} tParticle;

typedef struct tDotState {
  double xdot;
  double ydot;
  double vxdot;
  double vydot;
  double udot;
} tDotState;

#define MAX_BARRIERS 10
#define MAX_FILENAME_LENGTH 64

enum enum_barrier {
  barrier_wall_type = 1,
  barrier_ball_type = 2
};

typedef struct tBarrierData {
  int type;
  double x0;
  double y0;
  double x1;
  double y1;
  double r;
} tBarrierData;

enum enum_stepper {
  stepper_gpusph = 1,
  stepper_heun   = 2
};

enum enum_viscosity {
  viscosity_off      = 0,
  viscosity_monaghan = 1
};

typedef struct tSimParameters {
  double dt;
  int steps;
  int threads;
  int verbosity;
  int viscosity;
  int kernel_auto_type;
  double hash_load_factor;
  kernel_2d_func_ptr kernel_func;
  double kernel_eta;
  double kernel_h;
  double kernel_width;
  double gamma;
  double alpha;
  double beta;
  double eta;
  double gx;     // gravity vector (gx, gy)
  double gy;
  double epsbt;  // tangent direction barrier loss
  double epsbn;  // normal direction barrier loss
  char stepper_name[16];
  int stepper_type;
  char kernel_name[16];
  int trace_steps;
  char trace_file[MAX_FILENAME_LENGTH];
  int frame_steps;
  char frame_file[MAX_FILENAME_LENGTH];
  char final_file[MAX_FILENAME_LENGTH];
  char param_file[MAX_FILENAME_LENGTH];
  int num_barriers;
  tBarrierData barrier[MAX_BARRIERS];
} tSimParameters;

bool load_particlefile(const char* inputfilename, 
                       tParticle** ptr_ptr_particle,
                       tDotState** ptr_ptr_dotstate,
                       int* ptr_num_particles);

bool offload_particlefile(const char* outputfilename,
                          int num_particles,
                          tParticle** ptr_ptr_particle,
                          tDotState** ptr_ptr_dotstate);

bool setup_parameters(tSimParameters* P,
                      const tArgsio* A);

bool serialize_parameters(const tSimParameters* P,
                          const char* paramsfilename);

/* global simulation parameters */
tSimParameters SimParameters;

/* ------------------------------------------------------------ */
/* Timestep monitoring subprograms (h is constant)              */
/* ------------------------------------------------------------ */

// dt value defined through: 0.5*|acc|*(dt)^2 = h
double calc_timestep_delta_tf(double h,
                              int nump,
                              const tParticle* sp)
{
  double maxval = 0.0;
  for (int i = 0; i < nump; i++) {
    const double accx = sp[i].vxdot;
    const double accy = sp[i].vydot;
    const double acc2 = accx * accx + accy * accy;
    if (acc2 > maxval) maxval = acc2;
  }
  return (maxval > 0.0 ? sqrt(2.0 * h / sqrt(maxval)) : 0.0);
}

// dt value defined in terms of local speed of sound
double calc_timestep_delta_tcv(double h,
                               double alpha,
                               double beta,
                               int nump,
                               const tParticle* sp)
{
  double maxval = 0.0;
  for (int i = 0; i < nump; i++) {
    const double ci = sp[i].c;
    const double val = ci + 0.6 * (alpha * ci + beta * sp[i].mumax);
    if (val > maxval) maxval = val;
  }
  return (maxval > 0.0 ? h / maxval : 0.0);
}

/* ------------------------------------------------------------ */
/* Simulation state diagnostics                                 */
/* ------------------------------------------------------------ */

typedef struct tSimDiagnostics {
  double cmx;
  double cmy;
  double M;
  double Ek;
  double Eu;
  double Eg;
  double E;
  double px;
  double py;
  double Lz;
  double tf;
  double tcv;
} tSimDiagnostics;

void calc_diagnostics(tSimDiagnostics* diag,
                      int nump,
                      const tParticle* sp,
                      const tSimParameters* P)
{
  const double gx = P->gx;
  const double gy = P->gy;
  diag->cmx = 0.0;
  diag->cmy = 0.0;
  diag->M = 0.0;
  diag->Ek = 0.0;
  diag->Eu = 0.0;
  diag->Eg = 0.0;
  diag->E = 0.0;
  diag->px = 0.0;
  diag->py = 0.0;
  diag->Lz = 0.0;
  for (int i = 0; i < nump; i++) {
    const double mi = sp[i].m;
    diag->M   += mi;
    diag->cmx += mi * sp[i].x;
    diag->cmy += mi * sp[i].y;
    const double Ki = 0.5 * mi * (sp[i].vx * sp[i].vx + sp[i].vy * sp[i].vy);
    diag->Ek  += Ki;
    const double Ui = mi * sp[i].u;
    diag->Eu +=  Ui;
    const double Vi = -1.0 * mi * (gx * sp[i].x + gy * sp[i].y);
    diag->Eg  += Vi;
    diag->E   += Ki + Vi + Ui;
    diag->px  += mi * sp[i].vx;
    diag->py  += mi * sp[i].vy;
    diag->Lz  += mi * (sp[i].x * sp[i].vy - sp[i].y * sp[i].vx);  // (x,y,0) cross m*(vx,vy,0)
  }
  diag->cmx /= diag->M;
  diag->cmy /= diag->M;
  diag->tf = calc_timestep_delta_tf(P->kernel_h, nump, sp);
  diag->tcv = (P->viscosity != viscosity_off ? calc_timestep_delta_tcv(P->kernel_h, P->alpha, P->beta, nump, sp) : 0.0);
}

void trace_write_header(FILE* pf) {
  fprintf(pf, "# step,time,mass,cmx,cmy,energy,kin,int,pot,px,py,lz,tf,tcv\n");
}

void trace_write_row(FILE* pf, 
                     int64_t step, 
                     double time, 
                     const tSimDiagnostics* diag,
                     double dt)
{
  fprintf(pf, 
          "%jd,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
          step, time, diag->M, diag->cmx, diag->cmy, diag->E, diag->Ek, diag->Eu, diag->Eg, diag->px, diag->py, diag->Lz, diag->tf / dt, diag->tcv / dt);
}

/* ------------------------------------------------------------ */
/* Simulation state binary file dump                            */
/* ------------------------------------------------------------ */

void frame_write_header(FILE* pf, 
                        int nump)
{
  const int64_t np = (int64_t) nump;
  fwrite((const void*) &np, sizeof(int64_t), 1, pf);
}

void frame_write(FILE* pf, 
                 int64_t step,
                 double time,
                 int nump, 
                 const tParticle* sp)
{
  fwrite((const void*) &step, sizeof(int64_t), 1, pf);
  fwrite((const void*) &time, sizeof(double), 1, pf);
  for (int i = 0; i < nump; i++) {
    // Writes: m, x, y, vx, vy, u (see tParticle definition) 
    fwrite((const void*) &(sp[i].m), sizeof(double), 6, pf);
  }
}

/* ------------------------------------------------------------ */
/* Interaction callbacks                                        */
/* ------------------------------------------------------------ */

void density_summation_callback(int i, int j, void* aux) {
  tParticle* p = (tParticle*) aux;

  double w, wx, wy;
  const double xi = p[i].x;
  const double yi = p[i].y;

  const double xj = p[j].x;
  const double yj = p[j].y;

  (*SimParameters.kernel_func)(xi - xj, yi - yj, SimParameters.kernel_h, &w, &wx, &wy);

  p[i].rho += p[j].m * w;
}

void dot_summation_callback(int i, int j, void* aux) {
  if (i == j) return;

  tParticle* p = (tParticle*) aux;

  double w, wx, wy;
  const double pi = p[i].p;
  const double rhoi = p[i].rho;

  const double xi = p[i].x;
  const double yi = p[i].y;

  const double pj = p[j].p;
  const double rhoj = p[j].rho;

  const double xj = p[j].x;
  const double yj = p[j].y;

  const double dxij = xi - xj;
  const double dyij = yi - yj;

  const double h = SimParameters.kernel_h;

  (*SimParameters.kernel_func)(dxij, dyij, h, &w, &wx, &wy);

  const double mj = p[j].m;
  const double Ci = pi / (rhoi * rhoi);

  const double Aij = mj * (Ci + pj / (rhoj * rhoj));

  p[i].vxdot -= Aij * wx;
  p[i].vydot -= Aij * wy;

  const double vxi = p[i].vx;
  const double vyi = p[i].vy;

  const double vxj = p[j].vx;
  const double vyj = p[j].vy;

  const double dvxij = vxi - vxj;
  const double dvyij = vyi - vyj;

  p[i].udot += Ci * mj * (dvxij * wx + dvyij * wy);

  if (SimParameters.viscosity != viscosity_monaghan) return;

  const double vdotr = dvxij * dxij + dvyij * dyij;

  if (vdotr > 0.0) return;

  const double mu = (h * vdotr) / (dxij * dxij + dyij * dyij + SimParameters.eta * h * h);
  const double cbar = 0.5 * (p[i].c + p[j].c);
  const double rhobar = 0.5 * (rhoi + rhoj);

  const double Bij = mj * (-1.0 * SimParameters.alpha * cbar * mu + SimParameters.beta * mu * mu) / rhobar;

  p[i].vxdot -= Bij * wx;
  p[i].vydot -= Bij * wy;

  if (mu > p[i].mumax) p[i].mumax = mu;
}

/* ------------------------------------------------------------ */
/* Simulator elements                                           */
/* ------------------------------------------------------------ */

bool refresh_hash_index(tHashIndex2D* hti,
                        int nump, 
                        const tParticle* sp,
                        double supradius, 
                        bool testIndex)
{
  for (int i = 0; i < nump; i++) {
    hti->key_[i].xi = (int32_t) floor(sp[i].x / supradius);
    hti->key_[i].yi = (int32_t) floor(sp[i].y / supradius);
  }
  create_HashIndex2D(hti, nump); 
  return (testIndex ? test_HashIndex2D(hti) : true);
}

void refresh_rho_vdot_and_udot(const tHashIndex2D* hti,
                               int nump,
                               tParticle* sp,
                               tDotState* copy)
{
  const double gamma_minus_one = SimParameters.gamma - 1.0;
  const double sqrt_g_gm1 = sqrt(SimParameters.gamma * gamma_minus_one);

  #pragma omp parallel for
  for (int i = 0; i < nump; i++) {
    sp[i].rho = 0.0;
    singleInteract_HashIndex2D(hti, i, &density_summation_callback, sp);
    sp[i].p = gamma_minus_one * sp[i].u * sp[i].rho;
    if (SimParameters.viscosity == viscosity_off) sp[i].c = 0.0;
      else sp[i].c = sqrt_g_gm1 * sqrt(sp[i].u);  // equals sqrt(gamma*P/rho) = speed of sound
  }

  const double gx = SimParameters.gx;
  const double gy = SimParameters.gy;

  #pragma omp parallel for
  for (int i = 0; i < nump; i++) {
    sp[i].vxdot = 0.0;
    sp[i].vydot = 0.0;
    sp[i].udot = 0.0;
    sp[i].mumax = 0.0;
    singleInteract_HashIndex2D(hti, i, &dot_summation_callback, sp);
    sp[i].vxdot += gx;
    sp[i].vydot += gy;
  } 

  if (copy == NULL)
    return;
  
  for (int i = 0; i < nump; i++) {
    copy[i].xdot = sp[i].vx;
    copy[i].ydot = sp[i].vy;
    copy[i].vxdot = sp[i].vxdot;
    copy[i].vydot = sp[i].vydot;
    copy[i].udot = sp[i].udot;
  }
}

/* ------------------------------------------------------------ */
/* Basic predictor-corrector "stepping template"                */
/* ------------------------------------------------------------ */

// heuns:  c0 = dt
// gpusph: c0 = dt/2
void update_particles(double c0,
                      int nump,
                      tParticle* sp)
{
  for (int i = 0; i < nump; i++) {
    sp[i].x  += c0 * sp[i].vx;
    sp[i].y  += c0 * sp[i].vy;
    sp[i].vx += c0 * sp[i].vxdot;
    sp[i].vy += c0 * sp[i].vydot;
    sp[i].u  += c0 * sp[i].udot;
  }
}

// heuns:  c1 = dt/2, c2 = dt/2
// gpusph: c1 = dt,   c2 = -dt/2
void correct_particles(double c1,
                       double c2,
                       int nump,
                       tParticle* sp,
                       tDotState* dcpy)
{
  for (int i = 0; i < nump; i++) {
    sp[i].x  += c1 * sp[i].vx    + c2 * dcpy[i].xdot;
    sp[i].y  += c1 * sp[i].vy    + c2 * dcpy[i].ydot;
    sp[i].vx += c1 * sp[i].vxdot + c2 * dcpy[i].vxdot;
    sp[i].vy += c1 * sp[i].vydot + c2 * dcpy[i].vydot;
    sp[i].u  += c1 * sp[i].udot  + c2 * dcpy[i].udot;
  }
}

void setup_stepper_coeffs(int type,
                          double dt, 
                          double* coeffs)
{
  switch(type) {
    case stepper_heun:
      coeffs[0] = dt;
      coeffs[1] = dt / 2.0;
      coeffs[2] = -1.0 * dt / 2.0;
      break;
    case stepper_gpusph:
    default:
      coeffs[0] = dt / 2.0;
      coeffs[1] = dt;
      coeffs[2] = -1.0 * dt / 2.0;
      break;
  }
}

// generic ODE (u = full state vector):
//   du/dt = D(u)

// 2 evals per step (GPUSPH documentation):
//   uh = u + (dt/2) * D(u)
//   u = u + dt * D(uh)
//     = uh - (dt/2) * D(u) + dt * D(uh)
//     = uh + dt * (D(uh) - 0.5 * D(u))

// 2 evals per step (Heuns method)
//   uf = u + dt * D(u)
//   u = u + (dt/2) * (D(u) + D(uf)) 
//     = uf + (dt/2) * (D(uf) - D(u))

/* ------------------------------------------------------------ */
/* Practical subprogram to handle boundary surface              */
/* ------------------------------------------------------------ */

// boundary defined by two distinct points (x0, y0), (x1, y1)
// the vector (x1 - x0, y1 - y0) = (tx, ty) is tangent
// the normal direction is defined as: n = (ty, -tx)
// points on the outside of the boundary (sign defined by normal)
// and which are moving away from said boundary, have their
// normal velocity component reflected.
//
// {epn = 0, ept = 0} -> no energy loss
// 0 <= ep[t|n] <= 1 
// (implies energy losses in tangential/normal directions)

void adjust_wall_barrier(int nump,
                         tParticle* sp,
                         double x0, 
                         double y0,
                         double x1, 
                         double y1,
                         double ept, 
                         double epn)
{
  const double nfactor = -2.0 + epn;
  const double tfactor = -ept;
  const double tx = x1 - x0;
  const double ty = y1 - y0;
  const double tnorm = sqrt(tx * tx + ty * ty);
  if (tnorm == 0.0)
    return;
  const double thatx = tx / tnorm;
  const double thaty = ty / tnorm;
  const double nhatx = thaty;
  const double nhaty = -thatx;
  for (int i = 0; i < nump; i++) {
    const double xi = sp[i].x;
    const double yi = sp[i].y;
    const double rx = xi - x0;
    const double ry = yi - y0;
    const double beta = rx * nhatx + ry * nhaty;
    if (beta < 0.0)
      continue;
    // r = r0 + that * alfa + beta * nhat
    //const double alfa = rx * thatx + ry * thaty;
    const double gamma = sp[i].vx * nhatx + sp[i].vy * nhaty;
    if (gamma < 0.0)
      continue;
    sp[i].vx += nfactor * gamma * nhatx;
    sp[i].vy += nfactor * gamma * nhaty;
    if (tfactor == 0.0)
      continue;
    const double gammat = sp[i].vx * thatx + sp[i].vy * thaty;
    sp[i].vx += tfactor * gammat * thatx;
    sp[i].vy += tfactor * gammat * thaty;
  }
}

void adjust_ball_barrier(int nump,
                         tParticle* sp,
                         double x0, 
                         double y0,
                         double R,
                         bool outside,
                         double ept,
                         double epn)
{
  if (R == 0.0) return;
  const double normal_flip = (outside ? -1.0 : 1.0);
  const double nfactor = -2.0 + epn;
  const double tfactor = -ept;
  const double R2 = R * R;
  for (int i = 0; i < nump; i++) {
    const double xi = sp[i].x;
    const double yi = sp[i].y;
    const double rx = xi - x0;
    const double ry = yi - y0;
    const double rsq = rx * rx + ry * ry;
    if (!outside && rsq < R2) continue;
    if (outside && rsq > R2) continue;
    const double rnorm = sqrt(rsq);
    const double nhatx = normal_flip * rx / rnorm;
    const double nhaty = normal_flip * ry / rnorm;
    const double thatx = -nhaty;
    const double thaty = nhatx;
    const double gamma = sp[i].vx * nhatx + sp[i].vy * nhaty;
    if (gamma < 0.0) continue;
    sp[i].vx += nfactor * gamma * nhatx;
    sp[i].vy += nfactor * gamma * nhaty;
    if (tfactor == 0.0) continue;
    const double gammat = sp[i].vx * thatx + sp[i].vy * thaty;
    sp[i].vx += tfactor * gammat * thatx;
    sp[i].vy += tfactor * gammat * thaty;
  }
}

void apply_barriers(const tSimParameters* P,
                    int nump,
                    tParticle* sp)
{
  if (P->num_barriers == 0) return;
  for (int b = 0; b < P->num_barriers; b++) {
    switch (P->barrier[b].type) {
      case barrier_wall_type:
        adjust_wall_barrier(nump,
                            sp,
                            P->barrier[b].x0, 
                            P->barrier[b].y0, 
                            P->barrier[b].x1, 
                            P->barrier[b].y1, 
                            P->epsbt, 
                            P->epsbn);
        break;
      case barrier_ball_type:
        adjust_ball_barrier(nump,
                            sp,
                            P->barrier[b].x0, 
                            P->barrier[b].y0, 
                            fabs(P->barrier[b].r), 
                            (P->barrier[b].r < 0), 
                            P->epsbt, 
                            P->epsbn);
        break; 
    }
  }
  return;
}

/* ------------------------------------------------------------ */
/* Main CLI                                                     */
/* ------------------------------------------------------------ */

int main(int argc, const char** argv)
{
  if (argc < 3) {
    printf("usage: %s particlefile configfile [arg1=val1 ...]\n", argv[0]);
    return 1;
  }

  const char* particlefilename = argv[1];
  const char* configfilename = argv[2];

  tArgsio args;
  argsio_init(&args, 1024);
  const int nargs_from_file = argsio_import_file(&args, 
                                                 configfilename);

  if (nargs_from_file == 0) {
    printf("warning: no arguments read from \"%s\" (empty or missing file)\n", 
           configfilename);
  }

  const int num_cli_args = argc - 3;

  if (num_cli_args > 0) {
    const int num_cli_kv_added = argsio_add_cmdargs(&args, 
                                                    num_cli_args, 
                                                    &argv[3], 
                                                    false);
    int probed_verbosity = 0;
    const bool verbosity_specified = argsio_get_int(&args, "verbosity", &probed_verbosity);
    if ((verbosity_specified && probed_verbosity != 0) || !verbosity_specified) {
      printf("appended (or overwrote) %i %s from command line\n", 
             num_cli_kv_added, 
             (num_cli_kv_added == 1 ? "option" : "options"));
    }
  }

  bool has_valid_parameters = setup_parameters(&SimParameters, 
                                               &args);

  if (SimParameters.verbosity > 0) {
    argsio_printf(&args, "provided arguments");
  }

  if ((strcmp(SimParameters.final_file, particlefilename) == 0) ||
      (strcmp(SimParameters.final_file, configfilename) == 0)) {
    SimParameters.final_file[0] = '\0';
    printf("warning: disabled final-file since it equals an input file\n");
  }

  if (!has_valid_parameters) {
    printf("warning: failed to parse provided arguments\n");
  }

  argsio_uninit(&args);

  /*if (SimParameters.verbosity > 0 && SimParameters.num_barriers > 0) {
    printf("using %i barriers\n", SimParameters.num_barriers);
  }*/

  tParticle* ptr_particle = NULL;
  tDotState* ptr_dotstate = NULL;
  int num_particles = 0;

  const bool particles_are_loaded = load_particlefile(particlefilename, 
                                                      &ptr_particle, 
                                                      &ptr_dotstate,
                                                      &num_particles);

  if (!particles_are_loaded) {
    printf("warning: failed to load data from \"%s\"\n", 
           particlefilename);
  }

  if (has_valid_parameters && 
      particles_are_loaded && 
      SimParameters.kernel_h == 0.0)
  {
    autotune_kernel_bandwidth(SimParameters.kernel_eta,
                              SimParameters.kernel_auto_type,
                              SimParameters.kernel_func,
                              SimParameters.kernel_width,
                              num_particles, 
                              (const void *) &(ptr_particle[0].x),
                              sizeof(tParticle),
                              (const void *) &(ptr_particle[0].y),
                              sizeof(tParticle),
                              &(SimParameters.kernel_h),
                              SimParameters.verbosity > 1);

    has_valid_parameters = (has_valid_parameters) && (SimParameters.kernel_h > 0.0);
    if (has_valid_parameters && SimParameters.verbosity > 0) {
      printf("auto-tuned h = %e (eta = %f, type = %i, kernel = %s)\n", 
             SimParameters.kernel_h, 
             SimParameters.kernel_eta, 
             SimParameters.kernel_auto_type,
             SimParameters.kernel_name);
    }
  }

  const int max_threads = omp_get_max_threads();

  if (SimParameters.threads > max_threads) {
    printf("warning: thread usage capped at maximum = %i\n", max_threads);
    SimParameters.threads = max_threads;
  }

  if (strlen(SimParameters.param_file) != 0) {
    if (!serialize_parameters(&SimParameters, SimParameters.param_file)) {
      printf("warning: failed to write parameters file: \"%s\"\n", SimParameters.param_file);
    }
  }

  tHashIndex2D hti;
  const bool index_is_up = (allocate_HashIndex2D(&hti, 
                                                 num_particles, 
                                                 SimParameters.hash_load_factor) == 0);

  if (has_valid_parameters && 
      particles_are_loaded && 
      index_is_up) 
  {
    FILE *ptracefile = NULL;
    int trace_counter = 1;
    bool use_tracefile = (strlen(SimParameters.trace_file) > 0 && SimParameters.trace_steps > 0);
    if (use_tracefile) {
      ptracefile = fopen(SimParameters.trace_file, "w");
	    if (!ptracefile) use_tracefile = false;
        else trace_write_header(ptracefile);
    }

    FILE* pframefile = NULL;
    int frame_counter = 1;
    bool use_framefile = (strlen(SimParameters.frame_file) > 0 && SimParameters.frame_steps > 0);
    if (use_framefile) {
      pframefile = fopen(SimParameters.frame_file, "wb");
      if (!pframefile) use_framefile = false;
        else frame_write_header(pframefile, num_particles);
    }

    struct timespec tic, toc;
    struct timespec dotcalc_tic, dotcalc_toc;
    int64_t dotcalc_ticks = 0;

    tSimDiagnostics state_stats;

    double stepper_coeffs[3];
    setup_stepper_coeffs(SimParameters.stepper_type, SimParameters.dt, stepper_coeffs);

    if (SimParameters.threads != 0)
      omp_set_num_threads(SimParameters.threads);

    const double support_radius = SimParameters.kernel_h * SimParameters.kernel_width;
    double t = 0.0;

    clkutil_stamp(tic);

    for (int64_t k = 0; k < SimParameters.steps; k++)
    {
      if (use_framefile && --frame_counter == 0) {
        frame_write(pframefile, k, t, num_particles, ptr_particle);
        frame_counter = SimParameters.frame_steps;
      }

      refresh_hash_index(&hti,
                         num_particles, 
                         ptr_particle,
                         support_radius, 
                         false);

      clkutil_stamp(dotcalc_tic);

      refresh_rho_vdot_and_udot(&hti, 
                                num_particles,
                                ptr_particle, 
                                ptr_dotstate);

      clkutil_stamp(dotcalc_toc);
      dotcalc_ticks += clkutil_elapsed(dotcalc_tic, dotcalc_toc);

      // Place diagnostics calc after the particle vxdot, vydot, and udot have been computed for the state.
      // But must be before the particles are moved.
      calc_diagnostics(&state_stats,
                       num_particles,
                       ptr_particle,
                       &SimParameters);

      if (use_tracefile && --trace_counter == 0) {
        trace_write_row(ptracefile, k, t, &state_stats, SimParameters.dt);
        trace_counter = SimParameters.trace_steps;
      }

      update_particles(stepper_coeffs[0],
                       num_particles,
                       ptr_particle);
      
      refresh_hash_index(&hti,
                         num_particles, 
                         ptr_particle,
                         support_radius, 
                         false);

      clkutil_stamp(dotcalc_tic);

      refresh_rho_vdot_and_udot(&hti,
                                num_particles, 
                                ptr_particle, 
                                NULL);

      clkutil_stamp(dotcalc_toc);
      dotcalc_ticks += clkutil_elapsed(dotcalc_tic, dotcalc_toc);

      correct_particles(stepper_coeffs[1],
                        stepper_coeffs[2],
                        num_particles,
                        ptr_particle,
                        ptr_dotstate);
      
      apply_barriers(&SimParameters, 
                     num_particles, 
                     ptr_particle);

      t += SimParameters.dt;
    }

    clkutil_stamp(toc);

    const double total_time_seconds = ((double) clkutil_elapsed(tic, toc)) * 1.0e-9;
    const double total_dotcalc_elapsed = ((double) dotcalc_ticks) * 1.0e-9;
    if (SimParameters.verbosity > 0)
      printf("elapsed time = %.2f sec (dotcalc = %.2f sec)\n", total_time_seconds, total_dotcalc_elapsed);

    if (use_tracefile) {
      fclose(ptracefile);
      if (SimParameters.verbosity > 0)
        printf("closed \"%s\"\n", SimParameters.trace_file);
    }

    if (use_framefile) {
      fclose(pframefile);
      if (SimParameters.verbosity > 0)
        printf("closed \"%s\"\n", SimParameters.frame_file);
    }

    if (SimParameters.verbosity > 0)
      printf("stop time = %e (after %i steps with dt = %e)\n", t, SimParameters.steps, SimParameters.dt);
  }

  deallocate_HashIndex2D(&hti);

  offload_particlefile((strlen(SimParameters.final_file) == 0 ? NULL : SimParameters.final_file), 
                       num_particles,
                       &ptr_particle, 
                       &ptr_dotstate);

  return 0;
}

/* ------------------------------------------------------------ */
/* Parameter setup / default                                    */
/* ------------------------------------------------------------ */

bool setup_parameters(tSimParameters* P,
                      const tArgsio* A)
{
  memset(P, 0, sizeof(tSimParameters));

  P->dt = 1.0e-6;
  P->steps = 0;
  P->threads = 0;
  P->verbosity = 1;
  strcat(P->stepper_name, "gpusph");
  strcat(P->kernel_name, "quintic");
  P->kernel_func = NULL;
  P->kernel_auto_type = 1;
  P->hash_load_factor = 0.10;
  P->kernel_width = 0.0; // always defined through kernel_name below
  P->kernel_h = 0.0;     // zero implies auto-tune
  P->kernel_eta = 1.10;  // 1.20 ?
  P->gamma = 5.0 / 3.0;
  P->viscosity = viscosity_off;
  P->alpha = 1.0;
  P->beta = 2.0;
  P->eta = 0.01;  // alpha, beta, eta : parameters for viscosity feature
  P->epsbn = 0.100;
  P->epsbt = 0.025;
  P->trace_steps = 10;
  P->frame_steps = 100;
  P->gx = 0.0;
  P->gy = 0.0;

  if (!argsio_all_unique(A)) {
    printf("error: at least one argument name is specified multiple times\n");
    return false;
  }

  if (argsio_get_int(A, "verbosity", &(P->verbosity))) {
    if (P->verbosity < 0) return false;
  }

  if (argsio_get_int(A, "viscosity", &(P->viscosity))) {
    if (P->viscosity != viscosity_off && P->viscosity != viscosity_monaghan) return false;
  }

  if (argsio_get_real(A, "dt", &(P->dt))) {
    if (P->dt <= 0.0) return false;
  }

  // int via cast (from double) allows specifying large values snugly: e.g. "1.0e6"
  if (argsio_get_int_via_cast(A, "steps", &(P->steps))) {
    if (P->steps < 0) return false;
  }

  if (argsio_get_int(A, "threads", &(P->threads))) {
    if (P->threads < 0) return false;
  }

  if (argsio_get_int(A, "kernel-auto-type", &(P->kernel_auto_type))) {
    if (P->kernel_auto_type < 1 || P->kernel_auto_type > 2) return false;
  }

  if (argsio_get_real(A, "kernel-h", &(P->kernel_h))) {
    if (P->kernel_h < 0.0) return false; // 0 is allowed (and default) and means "autotune"
  }

  argsio_get_value(A, "stepper-name", P->stepper_name);
  P->stepper_type = (strcmp(P->stepper_name, "heun") == 0 ? stepper_heun : stepper_gpusph);

  if (P->stepper_type != stepper_gpusph && P->stepper_type != stepper_heun)
    return false;

  argsio_get_value(A, "kernel-name", P->kernel_name);
  for (int i = 0; i < NUM_AVAILABLE_KERNELS; i++) {
    if (strcmp(kernel_name_func_data[i].name, P->kernel_name) == 0) {
      P->kernel_func = kernel_name_func_data[i].fptr;
      P->kernel_width = kernel_name_func_data[i].width;
      break;
    }
  }

  if (P->kernel_func == NULL || P->kernel_width <= 0.0) {
    return false;
  }

  if (argsio_get_real(A, "gamma", &(P->gamma))) {
    if (P->gamma < 1.0 || P->gamma > 3.0) return false;
  }

  if (argsio_get_real(A, "kernel-eta", &(P->kernel_eta))) {
    if (P->kernel_eta < 0.5 || P->kernel_eta > 2.0) return false; // FIXME: check *correct* bounds 
  }

  // FIXME: this one could use a safety check (buffer length?)
  argsio_get_value(A, "trace-file", P->trace_file);

  if (argsio_get_int(A, "trace-steps", &(P->trace_steps))) {
    if (P->trace_steps <= 0 && strlen(P->trace_file) > 0 && P->verbosity > 0) {
      printf("warning: trace-steps <= 0 disables the provided trace-file\n");
    }
  }

  // FIXME: this one could use a safety check (buffer length?)
  argsio_get_value(A, "frame-file", P->frame_file);

  if (argsio_get_int(A, "frame-steps", &(P->frame_steps))) {
    if (P->frame_steps <= 0 && strlen(P->frame_file) > 0 && P->verbosity > 0) {
      printf("warning: frame-steps <= 0 disables the provided frame-file\n");
    }
  }

  if (argsio_get_real(A, "alpha", &(P->alpha))) {
    if (P->alpha < 0.0) return false;
  }

  if (argsio_get_real(A, "beta", &(P->beta))) {
    if (P->beta < 0.0) return false;
  }

  if (argsio_get_real(A, "eta", &(P->eta))) {
    if (P->eta < 0.0) return false;
  }

  if (argsio_get_real(A, "epsbn", &(P->epsbn))) {
    if (P->epsbn < 0.0 || P->epsbn > 1.0) return false;
  }

  if (argsio_get_real(A, "epsbt", &(P->epsbt))) {
    if (P->epsbt < 0.0 || P->epsbt > 1.0) return false;
  }

  // FIXME: these ones could use a safety check (buffer length?)
  argsio_get_value(A, "final-file", P->final_file);
  argsio_get_value(A, "param-file", P->param_file);

  // reused parsing buffers
  char textvalue[1024];
  const int vmax = 8;
  double vdata[vmax];

  if (argsio_lookup(A, "gravity") != NULL) {
    argsio_get_value(A, "gravity", textvalue);
    const int len = argsio_util_parse_vector(textvalue, vmax, vdata);
    if (len == 2) {
      P->gx = vdata[0];
      P->gy = vdata[1];
    } else {
      printf("ignored gravity (expected 2-vector gx,gy)\n");
    }
    if (isfinite(P->gx) == 0 || isfinite(P->gy) == 0) {
      printf("error: invalid gravity vector\n");
      return false;
    }
  }

  const int num_wall_barriers = argsio_count_prefix_matches(A, "barrier-wall:");

  P->num_barriers = 0;
  for (int i = 0; i < num_wall_barriers; i++) {
    if (P->num_barriers == MAX_BARRIERS) {
      printf("skipping wall barrier no. %i\n", i);
      continue;
    }
    argsio_get_prefix_value(A, "barrier-wall:", i, textvalue);
    const int len = argsio_util_parse_vector(textvalue, vmax, vdata);
    if (len != 4) {
      printf("ignored (expected 4-vector x0,y0,x1,y1): \"%s\"\n", textvalue);
      continue;
    }
    P->barrier[P->num_barriers].type = barrier_wall_type;
    P->barrier[P->num_barriers].x0   = vdata[0];
    P->barrier[P->num_barriers].y0   = vdata[1];
    P->barrier[P->num_barriers].x1   = vdata[2];
    P->barrier[P->num_barriers].y1   = vdata[3];
    P->num_barriers++;
  }

  const int num_ball_barriers = argsio_count_prefix_matches(A, "barrier-ball:");

  for (int i = 0; i < num_ball_barriers; i++) {
    if (P->num_barriers == MAX_BARRIERS) {
      printf("skipping ball barrier no. %i\n", i);
      continue;
    }
    argsio_get_prefix_value(A, "barrier-ball:", i, textvalue);
    const int len = argsio_util_parse_vector(textvalue, vmax, vdata);
    if (len != 3) {
      printf("ignored (expected 3-vector cx,cy,r): \"%s\"\n", textvalue);
      continue;
    }
    P->barrier[P->num_barriers].type = barrier_ball_type;
    P->barrier[P->num_barriers].x0   = vdata[0];
    P->barrier[P->num_barriers].y0   = vdata[1];
    P->barrier[P->num_barriers].r    = vdata[2];
    P->num_barriers++;
  }

  //printf("[%s]: num barriers = %i\n", __func__, P->num_barriers);
  for (int i = 0; i < P->num_barriers; i++) {
    bool isok = (isfinite(P->barrier[i].x0) != 0);
    isok = isok && (isfinite(P->barrier[i].y0) != 0);
    isok = isok && (isfinite(P->barrier[i].x1) != 0);
    isok = isok && (isfinite(P->barrier[i].y1) != 0);
    isok = isok && (isfinite(P->barrier[i].r) != 0);
    isok = isok && (P->barrier[i].type == barrier_wall_type || 
                    P->barrier[i].type == barrier_ball_type);
    if (!isok) {
      printf("error: invalid barrier specification(s)\n");
      return false;
    }
    //printf("r[%i] = %e\n", i, P->barrier[i].r);
  }

  return true;
}

bool serialize_parameters(const tSimParameters* P,
                          const char* paramsfilename)
{
  if (P == NULL) return false;

  char buffer[1024];

  tArgsio A;
  argsio_init(&A, 1024);

  sprintf(buffer, "steps=%i", P->steps);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "verbosity=%i", P->verbosity);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "viscosity=%i", P->viscosity);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "dt=%.16e", P->dt);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "threads=%i", P->threads);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "trace-steps=%i", P->trace_steps);
  argsio_add_kv(&A, buffer);

  if (strlen(P->trace_file) > 0) {
    sprintf(buffer, "trace-file=%s", P->trace_file);
    argsio_add_kv(&A, buffer);
  }

  sprintf(buffer, "frame-steps=%i", P->frame_steps);
  argsio_add_kv(&A, buffer);

  if (strlen(P->frame_file) > 0) {
    sprintf(buffer, "frame-file=%s", P->frame_file);
    argsio_add_kv(&A, buffer);
  }

  if (strlen(P->final_file) > 0) {
    sprintf(buffer, "final-file=%s", P->final_file);
    argsio_add_kv(&A, buffer);
  }

  sprintf(buffer, "kernel-name=%s", P->kernel_name);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "kernel-h=%.16e", P->kernel_h);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "stepper-name=%s", P->stepper_name);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "gamma=%.16e", P->gamma);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "kernel-eta=%.16e", P->kernel_eta);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "alpha=%.16e", P->alpha);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "beta=%.16e", P->beta);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "eta=%.16e", P->eta);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "epsbn=%.16e", P->epsbn);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "epsbt=%.16e", P->epsbt);
  argsio_add_kv(&A, buffer);

  sprintf(buffer, "gravity=%.16e,%.16e", P->gx, P->gy);
  argsio_add_kv(&A, buffer);

  for (int i = 0; i < P->num_barriers; i++) {
    switch (P->barrier[i].type) {
      case barrier_wall_type:
        sprintf(buffer, 
                "barrier-wall:%i=%.16e,%.16e,%.16e,%.16e", 
                i, P->barrier[i].x0, P->barrier[i].y0, P->barrier[i].x1, P->barrier[i].y1);
        argsio_add_kv(&A, buffer);
        break;
      case barrier_ball_type:
        sprintf(buffer, 
                "barrier-ball:%i=%.16e,%.16e,%.16e", 
                i, P->barrier[i].x0, P->barrier[i].y0, P->barrier[i].r);
        argsio_add_kv(&A, buffer);
        break;
    }
  }

  if (!argsio_all_unique(&A)) {
    printf("warning: non-unique key(s) after serialization\n");
  }

  bool wroteit = false;
  if (paramsfilename != NULL) {
    const char header[] = "serialized parameters";
    wroteit = argsio_export_file(&A, paramsfilename, header);
  }

  argsio_uninit(&A);
  return wroteit;
}

/* ------------------------------------------------------------ */
/* Basic I/O                                                    */
/* ------------------------------------------------------------ */

bool load_particlefile(const char* inputfilename, 
                       tParticle** ptr_ptr_particle,
                       tDotState** ptr_ptr_dotstate,
                       int* ptr_num_particles)
{
  *ptr_ptr_particle = NULL;
  *ptr_ptr_dotstate = NULL;
  *ptr_num_particles = 0;

  const char sep_char = ',';
  const bool detectComments = true;
  double* read_data = NULL;
  int read_m = 0, read_n = 0;
  const int read_return_code = textio_read_double_matrix(inputfilename, 
                                                         &read_data, 
                                                         &read_m, 
                                                         &read_n, 
                                                         TEXTIO_COLUMNMAJOR, 
                                                         sep_char,
                                                         detectComments);

  bool read_success = (read_data != NULL && 
                       read_n == 6 && 
                       read_m > 0 &&
                       read_return_code == read_m * read_n);

  for (int e = 0; e < read_m * read_n; e++) {
    read_success = (read_success && (isfinite(read_data[e]) != 0));
  }

  // this is the REQUIRED column order: x, y, vx, vy, u, m
  if (read_success) {
    tParticle* P = malloc(read_m * sizeof(tParticle));
    memset(P, 0, read_m * sizeof(tParticle));
    for (int i = 0; i < read_m; i++) {
      P[i].x  = read_data[read_m * 0 + i];
      P[i].y  = read_data[read_m * 1 + i];
      P[i].vx = read_data[read_m * 2 + i];
      P[i].vy = read_data[read_m * 3 + i];
      P[i].u  = read_data[read_m * 4 + i];
      P[i].m  = read_data[read_m * 5 + i];
    }
    *ptr_ptr_particle = P;

    tDotState* Q = malloc(read_m * sizeof(tDotState));
    memset(Q, 0, read_m * sizeof(tDotState));
    *ptr_ptr_dotstate = Q;

    *ptr_num_particles = read_m;
  }

  if (read_data != NULL)
    free(read_data);

  if (read_success) {
    const tParticle* P = *ptr_ptr_particle;
    for (int i = 0; i < read_m; i++) {
      read_success = (read_success && (P[i].m > 0.0) && (P[i].u > 0.0));
    }
  }

  return (read_success && 
          *ptr_ptr_particle != NULL && 
          *ptr_ptr_dotstate != NULL);
}

bool offload_particlefile(const char* outputfilename,
                          int num_particles,
                          tParticle** ptr_ptr_particle,
                          tDotState** ptr_ptr_dotstate)
{
  if (outputfilename != NULL && *ptr_ptr_particle != NULL) {
    double* M = malloc(num_particles * 6 * sizeof(double));
    const tParticle* sp = *ptr_ptr_particle;

    for (int i = 0; i < num_particles; i++) {
      M[i + num_particles * 0] = sp[i].x;
      M[i + num_particles * 1] = sp[i].y;
      M[i + num_particles * 2] = sp[i].vx;
      M[i + num_particles * 3] = sp[i].vy;
      M[i + num_particles * 4] = sp[i].u;
      M[i + num_particles * 5] = sp[i].m;
    }

    const char sep_char = ',';
    const char head_comment[] = "x,y,vx,vy,u,m";
    const char foot_comment[] = "autogenerated final state";

    textio_write_double_matrix(outputfilename,
                               M,
												  		 num_particles,
															 6,
															 NULL,
															 sep_char,
															 head_comment,
															 foot_comment);
 
    free(M);
  }

  if (*ptr_ptr_particle != NULL) {
    free(*ptr_ptr_particle);
    *ptr_ptr_particle = NULL;
  }

  if (*ptr_ptr_dotstate != NULL) {
    free(*ptr_ptr_dotstate);
    *ptr_ptr_dotstate = NULL;
  }

  return true;
}
