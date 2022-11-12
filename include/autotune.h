#ifndef __AUTOTUNE_H__
#define __AUTOTUNE_H__

// TODO: refactor to remove these globals
double autotune_kernel_h;
kernel_2d_func_ptr autotune_kernel_func;

/* convenient copy of coordinates when running the bandwidth autotuner */
typedef struct tAutotunePointData {
  double x;
  double y;
  double w;
} tAutotunePointData;

void autotune_number_density_callback(int i, 
                                      int j, 
                                      void* aux) 
{
  tAutotunePointData* p = (tAutotunePointData*) aux;

  double w, wx, wy;
  const double xi = p[i].x;
  const double yi = p[i].y;

  const double xj = p[j].x;
  const double yj = p[j].y;

  (*autotune_kernel_func)(xi - xj, yi - yj, autotune_kernel_h, &w, &wx, &wy);

  p[i].w += w;
}

// Evaluate the sum of all kernels at each kernel center; for a given bandwidth h
bool autotune_number_density(double h,
                             kernel_2d_func_ptr K,
                             double width,
                             int numpts,
                             tHashIndex2D* hti,
                             tAutotunePointData* pts)
{
  const double support_radius = width * h;
  for (int i = 0; i < numpts; i++) {
    hti->key_[i].xi = (int32_t) floor(pts[i].x / support_radius);
    hti->key_[i].yi = (int32_t) floor(pts[i].y / support_radius);
  }

  create_HashIndex2D(hti, numpts);
  if (!test_HashIndex2D(hti))
    return false;

  autotune_kernel_h = h;
  autotune_kernel_func = K;

  for (int i = 0; i < numpts; i++) {
    pts[i].w = 0.0;
    singleInteract_HashIndex2D(hti, i, &autotune_number_density_callback, pts);
  }

  return true;
}

double autotune_average_type1(int numpts,
                              tAutotunePointData* pts)
{
  double sum_inv = 0.0;
  for (int i = 0; i < numpts; i++)
    sum_inv += 1.0 / sqrt(pts[i].w);
  const double avg_inv = sum_inv / numpts;
  return avg_inv;
}

double autotune_average_type2(int numpts,
                              tAutotunePointData* pts)
{
  double sum_sum = 0.0;
  for (int i = 0; i < numpts; i++)
    sum_sum += pts[i].w;
  const double avg_sum = sum_sum / numpts;
  return (1.0 / sqrt(avg_sum));
}

typedef double (*autotune_average_func_ptr)(int, tAutotunePointData*);

bool autotune_kernel_bandwidth(double eta,
                               int avgtype,
                               kernel_2d_func_ptr K,
                               double width,
                               int numpts,
                               const void* x,
                               int incx,
                               const void* y,
                               int incy,
                               double* auto_h,
                               bool debug_printf)
{
  if (K == NULL || x == NULL || y == NULL || width <= 0.0) return false;
  if (avgtype != 1 && avgtype != 2) return false;
  autotune_average_func_ptr autotune_avg = (avgtype == 1 ? autotune_average_type1 : autotune_average_type2);

  tAutotunePointData* auto_data = malloc(numpts * sizeof(tAutotunePointData));
  if (auto_data == NULL) return false;

  double xmin = *((const double*)((const char*)x + 0 * incx));
  double xmax = xmin;

  double ymin = *((const double*)((const char*)y + 0 * incy));
  double ymax = ymin;

  for (int i = 0; i < numpts; i++) {
    const double xi = *((const double*)((const char*)x + i * incx));
    if (xi < xmin) xmin = xi;
    if (xi > xmax) xmax = xi;
    auto_data[i].x = xi;
    const double yi = *((const double*)((const char*)y + i * incy));
    if (yi < ymin) ymin = yi;
    if (yi > ymax) ymax = yi;
    auto_data[i].y = yi;
  }

  if (debug_printf) {
    printf("xmin, xmax = %f, %f\n", xmin, xmax);
    printf("ymin, ymax = %f, %f\n", ymin, ymax);
  }

  const double xrange = xmax - xmin;
  const double yrange = ymax - ymin;
  const double hbar = (xrange + yrange) / 2.0 / 2.0;
  const double hzero = hbar / numpts;

  if (debug_printf)
    printf("hbar, hzero = %f, %f\n", hbar, hzero);

  tHashIndex2D hti;
  const double load_factor = 0.10;
  if (allocate_HashIndex2D(&hti, numpts, load_factor) != 0) {
    free(auto_data);
    return false;
  }

  const int itr_max = 50;

  int itr = 0;
  double h = hzero;

  // upper bound h

  while (++itr < itr_max) {
    if (!autotune_number_density(h,
                                 K,
                                 width,
                                 numpts,
                                 &hti,
                                 auto_data)) break;

    const double hdens = eta * (*autotune_avg)(numpts, auto_data);

    if (debug_printf)
      printf("[bound %i]: hdens = %e (h = %e)\n", itr, hdens, h);

    if (hdens < h) break;
    h *= 2.0;
  }

  if (itr == itr_max) {
    deallocate_HashIndex2D(&hti);
    free(auto_data);
    return false;
  }

  // bisect bracket [h/2, h]

  const double epstol = 1.0e-6;
  double ha = h / 2.0;
  double hb = h;
  itr = 0;

  while (itr++ < itr_max) {
    h = (ha + hb) / 2.0;
    if ((hb - ha) / h < epstol) break;

    if (!autotune_number_density(h,
                                 K,
                                 width,
                                 numpts,
                                 &hti,
                                 auto_data)) break;

    const double hdens = eta * (*autotune_avg)(numpts, auto_data);

    if (debug_printf)
      printf("[bisect %i]: hdens = %e (h = %e)\n", itr, hdens, h);

    if (hdens < h) hb = h;
      else ha = h;
  }

  deallocate_HashIndex2D(&hti);
  free(auto_data);

  if (auto_h != NULL)
    *auto_h = (itr < itr_max ? h : 0.0);

  if (debug_printf)
    printf("[%s]: exiting with auto_h = %f (success = %i)\n", __func__, h, itr < itr_max);

  return (itr < itr_max); 
}

#endif
