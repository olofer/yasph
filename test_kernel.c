#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "kernels2d.h"

double trapz_xy(kernel_2d_func_ptr kxyh,
                double h,
                int nx, 
                double xa, 
                double xb,
                int ny, 
                double ya, 
                double yb)
{
  const double hx = (xb - xa) / (nx - 1);
  const double hy = (yb - ya) / (ny - 1);
  double sum = 0.0;
  double kij = 0.0, kx = 0.0, ky = 0.0;
  for (int i = 0; i < nx; i++) {
    const double xi = i * hx + xa;
    const double wi = hx * ((i == 0 || i == (nx - 1)) ? 0.5 : 1.0);
    for (int j = 0; j < ny; j++) {
      const double yj = j * hy + ya;
      const double wj = hy * ((j == 0 || j == (ny - 1)) ? 0.5 : 1.0);
      (*kxyh)(xi, yj, h, &kij, &kx, &ky);
      sum += wi * wj * kij;
    }
  }
  return sum;
}

// Straight line integral from outside support (px, py) into center of kernel (origin).
// The integrand is the dot product of the kernel gradient and the directional unit vector.
// So the true value of the integral must be the (center) kernel value at (0, 0).

double line_trapz(kernel_2d_func_ptr kxyh,
                  double h,
                  int ns,
                  double px,
                  double py,
                  double* K0)
{
  const double L = sqrt(px * px + py * py);
  const double tx = px / L;
  const double ty = py / L;
  const double hs = L / (ns - 1);
  double sum = 0.0;
  double ki = 0.0, kix = 0.0, kiy = 0.0;
  for (int i = 0; i < ns; i++) {
    const double si = i * hs;        // s = 0 .. L
    const double xi = px - si * tx;  // dx/ds = -tx
    const double yi = py - si * ty;  // dy/ds = -ty
    (*kxyh)(xi, yi, h, &ki, &kix, &kiy);
    const double wi = hs * ((i == 0 || i == (ns - 1)) ? 0.5 : 1.0);
    sum += wi * (kix * tx + kiy * ty) * (-1.0);
  }
  if (K0 != NULL) *K0 = ki;
  return sum;
}

int main(int argc, const char** argv)
{
  if (argc != 4) {
    printf("usage: %s kernelname nabscissas bandwidth\n", argv[0]);
    return 1;
  }

  const char* kernelname = argv[1];

  kernel_2d_func_ptr kxyh = NULL;
  double width = 0.0;

  for (int i = 0; i < NUM_AVAILABLE_KERNELS; i++) {
    if (strcmp(kernel_name_func_data[i].name, kernelname) == 0) {
      kxyh = kernel_name_func_data[i].fptr;
      width = kernel_name_func_data[i].width;
      break;
    }
  }

  if (kxyh == NULL) {
    printf("did not recognize kernelname = \"%s\"\n", kernelname);
    return 1;
  }

  const int nabsc = atoi(argv[2]);

  if (nabsc < 2) {
    printf("at least 2 abscissas required\n");
    return 1;
  }

  const double h = atof(argv[3]);

  if (h <= 0.0) {
    printf("bandwidth > 0 required\n");
    return 1;
  }

  const double R = h * width;  // actual support radius
  const double I0 = trapz_xy(kxyh, h, nabsc, -R, R, nabsc, -R, R); // nabsc^2 evals

  printf("I0[\"%s\" | h = %e] = %.16e\n", kernelname, h, I0);

  double I1_check = 0.0;
  const double theta = -1.45435324;
  const double I1 = line_trapz(kxyh, h, nabsc * nabsc, R * cos(theta), R * sin(theta), &I1_check); // nabsc^2 evals

  printf("I1[\"%s\" | h = %e] = %.16e; I1-check = %.16e\n", kernelname, h, I1, I1_check);

  const double test_tol_0 = 1.0e-8;
  const double test_tol_1 = 1.0e-8;

  const bool I0_ok = fabs(I0 - 1.0) < test_tol_0;
  const bool I1_ok = fabs(I1 - I1_check) / I1_check < test_tol_1;

  return ((I0_ok && I1_ok) ? 0 : -1);
}
