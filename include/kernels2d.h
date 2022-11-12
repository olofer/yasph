#ifndef __KERNELS2D_H__
#define __KERNELS2D_H__

/*
 * 2D space exclusively; contains:
 *
 *   spline kernels
 *   wendland kernels
 *   gaussian
 *   super gaussian
 *
 */

typedef void (*kernel_2d_func_ptr)(double, double, double, double*, double*, double*);

/* ------------------------------------------------------------ */
/* Uniform kernel (mostly for reference) support radius 1*h     */
/* ------------------------------------------------------------ */

void evaluate_kernel_uf(double x, 
                        double y, 
                        double h, 
                        double* w,
                        double* wx,
                        double* wy)
{
  const double _one_pi = 3.14159265358979323846;
  const double h2 = h * h;
  const double sigma = 1.0 / (h2 * _one_pi);
  const double q2 = (x * x + y * y) / h2;
  *w = (q2 <= 1.0 ? sigma : 0.0);
  *wx = 0.0;
  *wy = 0.0;
  return;
}

/* ------------------------------------------------------------ */
/* Gaussian 2D Interpolation Kernel (infinite support radius)   */
/* ------------------------------------------------------------ */

void evaluate_kernel_g(double x, 
                       double y, 
                       double h, 
                       double* w,
                       double* wx,
                       double* wy)
{
  const double _one_pi = 3.14159265358979323846;
  const double sigma = 1.0 / _one_pi;
  const double h2 = h * h;
  const double q2 = (x * x + y * y) / h2;
  const double N = sigma / h2;
  const double W = N * exp(-q2);
  *w = W;
  *wx = -2.0 * x * W / h2;
  *wy = -2.0 * y * W / h2;
  return;
}

/* ------------------------------------------------------------ */
/* "Super Gaussian" (inf. support, and goes negative)           */
/* ------------------------------------------------------------ */

void evaluate_kernel_sg(double x, 
                        double y, 
                        double h, 
                        double* w,
                        double* wx,
                        double* wy)
{
  const double _one_pi = 3.14159265358979323846;
  const double sigma = 1.0 / _one_pi;
  const double h2 = h * h;
  const double N = sigma / h2;
  const double q2 = (x * x + y * y) / h2;
  const double emq2 = exp(-q2);
  *w = N * emq2 * (2.0 - q2);
  const double Wprime = -2.0 * (N / h2) * emq2 * (3.0 - q2);
  *wx = Wprime * x;
  *wy = Wprime * y;
  return;
}

/* ------------------------------------------------------------ */
/* Wendland C2; support = 2*h                                   */
/* ------------------------------------------------------------ */

void evaluate_kernel_wc2(double x, 
                         double y, 
                         double h, 
                         double* w,
                         double* wx,
                         double* wy)
{
  const double _one_pi = 3.14159265358979323846;
  const double sigma = 7.0 / (64.0 * _one_pi);
  const double h2 = h * h;
  const double N = sigma / h2;
  const double q2 = (x * x + y * y) / h2;
  if (q2 > 4.0) {
    *w = 0.0;
    *wx = 0.0;
    *wy = 0.0;
    return;
  }
  const double q = sqrt(q2);
  //const double tmp = (2.0 - q) * (2.0 - q);
  const double tmp = 4.0 - 4.0 * q + q2;
  const double W = N * (1.0 + 2.0 * q) * tmp * tmp;
  *w = W;
  //if (q == 0) {
  //  *wx = 0.0;
  //  *wy = 0.0;
  //  return;
  //}
  //const double Wprime = (N / (q * h2)) * (-10.0 * q * (2.0 - q) * tmp);
  const double Wprime = -10.0 * (N / h2) * (2.0 - q) * tmp;
  *wx = Wprime * x;
  *wy = Wprime * y;
  return;
}

/* ------------------------------------------------------------ */
/* Wendland C4; support = 2*h                                   */
/* ------------------------------------------------------------ */

void evaluate_kernel_wc4(double x, 
                         double y, 
                         double h, 
                         double* w,
                         double* wx,
                         double* wy)
{
  const double _one_pi = 3.14159265358979323846;
  const double sigma = 9.0 / (4.0 * _one_pi);
  const double c0 = 35.0 / 12.0;
  const double c1 = 7.0 / 96.0;
  const double h2 = h * h;
  const double N = sigma / h2;
  const double q2 = (x * x + y * y) / h2;
  if (q2 > 4.0) {
    *w = 0.0;
    *wx = 0.0;
    *wy = 0.0;
    return;
  }
  const double q = sqrt(q2);
  const double tmp = 0.25 * q2 - q + 1.0; 
  const double W = N * (tmp * tmp * tmp) * (c0 * q2 + 3.0 * q + 1.0);
  *w = W;
  //const double Wprime = c1 * (N / h2) * (5.0 * q2 * q2 * q2 - 48.0 * q2 * q2 * q + 180.0 * q2 * q2 - 320.0 * q2 * q + 240.0 * q2 - 64.0);
  const double Wprime = c1 * (N / h2) * (q2 * (q2 * (5.0 * q2 - 48.0 * q + 180.0) - 320.0 * q + 240.0) - 64.0);
  *wx = Wprime * x;
  *wy = Wprime * y;
  return;
}

/* ------------------------------------------------------------ */
/* Wendland C6; support = 2*h                                   */
/* ------------------------------------------------------------ */

void evaluate_kernel_wc6(double x, 
                         double y, 
                         double h, 
                         double* w,
                         double* wx,
                         double* wy)
{
  const double _one_pi = 3.14159265358979323846;
  const double sigma = 39.0 / (28672.0 * _one_pi);
  const double h2 = h * h;
  const double N = sigma / h2;
  const double q2 = (x * x + y * y) / h2;
  if (q2 > 4.0) {
    *w = 0.0;
    *wx = 0.0;
    *wy = 0.0;
    return;
  }
  const double q = sqrt(q2);
  const double q3 = q * q2;
  const double tmp = 4.0 - 4.0 * q + q2; // (2 - q)^2
  const double tmp3 = tmp * tmp * tmp; 
  const double W = N * (tmp3 * tmp) * (8.0 + 32.0 * q + 50.0 * q2 + 32.0 * q3);
  *w = W;
  // -44 (2 - q)^7 q (2 + 7 q + 8 q^2)
  const double Wprime = -44.0 * (N / h2) * tmp3 * (2.0 - q) * (2.0 + 7.0 * q + 8.0 * q2);
  *wx = Wprime * x;
  *wy = Wprime * y;
  return;
}

/* ------------------------------------------------------------ */
/* Quintic 2D Interpolation Kernel (support radius = 3*h)       */
/* ------------------------------------------------------------ */

void evaluate_kernel_5(double x, 
                       double y, 
                       double h, 
                       double* w,
                       double* wx,
                       double* wy)
{
  const double _one_pi = 3.14159265358979323846;
  const double sigma = 7.0 / (478.0 * _one_pi);
  const double q = sqrt(x * x + y * y) / h;
  const double h2 = h * h;
  const double N = sigma / h2;
  const double Nq = N / h2;

  if (q >= 3.0) {
    *w = 0.0;
    *wx = 0.0;
    *wy = 0.0;
    return;
  }

  const double c1 = 3.0 - q;
  const double t14 = c1 * c1 * c1 * c1;
  const double t15 = t14 * c1; 

  if (q >= 2.0) {
    *w = N * t15;
    const double wprime = Nq * (-5.0 * t14) / q;
    *wx = wprime * x;
    *wy = wprime * y;
    return;
  }

  const double c2 = 2.0 - q;
  const double t24 = c2 * c2 * c2 * c2;
  const double t25 = t24 * c2;

  if (q >= 1.0) {
    *w = N * (t15 - 6.0 * t25);
    const double wprime = Nq * (-5.0 * t14 + 30.0 * t24) / q;
    *wx = wprime * x;
    *wy = wprime * y;
    return;
  }

  double c3 = 1.0 - q;
  double t34 = c3 * c3 * c3 * c3;
  double t35 = t34 * c3;

  // else q >= 0 && q < 1.0

  *w = N * (t15 - 6.0 * t25 + 15.0 * t35);

  if (q != 0.0) {
    const double wprime = Nq * (-5.0 * t14 + 30.0 * t24 - 75.0 * t34) / q;
    *wx = wprime * x;
    *wy = wprime * y;
    return;
  }

  *wx = 0.0;
  *wy = 0.0;
  return;
}

/* ------------------------------------------------------------ */
/* Quartic 2D Interpolation Kernel (support radius = 2.5*h)     */
/* ------------------------------------------------------------ */

void evaluate_kernel_4(double x, 
                       double y, 
                       double h, 
                       double* w,
                       double* wx,
                       double* wy)
{
  const double _one_pi = 3.14159265358979323846;
  const double sigma = 96.0 / (1199.0 * _one_pi);
  const double q = sqrt(x * x + y * y) / h;
  const double h2 = h * h;
  const double N = sigma / h2;
  const double Nq = N / h2;

  if (q >= 2.5) {
    *w = 0.0;
    *wx = 0.0;
    *wy = 0.0;
    return;
  }

  const double c1 = 2.5 - q;
  const double t13 = c1 * c1 * c1;
  const double t14 = t13 * c1; 

  if (q >= 1.5) {
    *w = N * t14;
    const double wprime = Nq * (-4.0 * t13) / q;
    *wx = wprime * x;
    *wy = wprime * y;
    return;
  }

  const double c2 = 1.5 - q;
  const double t23 = c2 * c2 * c2;
  const double t24 = t23 * c2;

  if (q >= 0.5) {
    *w = N * (t14 - 5.0 * t24);
    const double wprime = Nq * (-4.0 * t13 + 20.0 * t23) / q;
    *wx = wprime * x;
    *wy = wprime * y;
    return;
  }

  double c3 = 0.5 - q;
  double t33 = c3 * c3 * c3;
  double t34 = t33 * c3;

  // else q >= 0 && q < 0.5

  *w = N * (t14 - 5.0 * t24 + 10.0 * t34);

  if (q != 0.0) {
    const double wprime = Nq * (-4.0 * t13 + 20.0 * t23 - 40.0 * t33) / q;
    *wx = wprime * x;
    *wy = wprime * y;
    return;
  }

  *wx = 0.0;
  *wy = 0.0;
  return;
}

/* ------------------------------------------------------------ */
/* Cubic 2D Interpolation Kernel (support radius = 2*h)         */
/* ------------------------------------------------------------ */

void evaluate_kernel_3(double x, 
                       double y, 
                       double h, 
                       double* w,
                       double* wx,
                       double* wy)
{
  const double _one_pi = 3.14159265358979323846;
  const double sigma = 10.0 / (7.0 * _one_pi);
  const double q = sqrt(x * x + y * y) / h;
  const double h2 = h * h;
  const double N = sigma / h2;
  const double Nq = N / h2;

  if (q >= 2.0) {
    *w = 0.0;
    *wx = 0.0;
    *wy = 0.0;
    return;
  }

  const double c1 = 2.0 - q;
  const double t12 = c1 * c1;
  const double t13 = t12 * c1; 

  if (q >= 1.0) {
    *w = N * 0.25 * t13;
    const double wprime = Nq * (-0.75 * t12) / q;
    *wx = wprime * x;
    *wy = wprime * y;
    return;
  }

  // else q >= 0.0 and q < 1.0
  const double c2 = 1.0 - q;
  const double t22 = c2 * c2;
  const double t23 = t22 * c2;

  *w = N * (0.25 * t13 - 1.0 * t23);
  if (q != 0.0) {
    const double wprime = Nq * (-0.75 * t12 + 3.0 * t22) / q;
    *wx = wprime * x;
    *wy = wprime * y;
    return; 
  }

  *wx = 0.0;
  *wy = 0.0;
  return;
}

/* ------------------------------------------------------------ */

typedef struct tKernNamePtr {
  char name[32];
  kernel_2d_func_ptr fptr;
  double width;
} tKernNamePtr;

#define NUM_AVAILABLE_KERNELS 9

const tKernNamePtr kernel_name_func_data[NUM_AVAILABLE_KERNELS] = { 
  {"uniform",    &evaluate_kernel_uf,  1.0 },
  {"cubic",      &evaluate_kernel_3,   2.0 },
  {"quartic",    &evaluate_kernel_4,   2.5 },
  {"quintic",    &evaluate_kernel_5,   3.0 },
  {"wendlandc2", &evaluate_kernel_wc2, 2.0 },
  {"wendlandc4", &evaluate_kernel_wc4, 2.0 },
  {"wendlandc6", &evaluate_kernel_wc6, 2.0 },
  {"gaussian",   &evaluate_kernel_g,   6.0 },
  {"sgaussian",  &evaluate_kernel_sg,  6.0 } 
};

#endif
