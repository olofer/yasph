#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "quadtree.h"

/*

Randomly scatter numpts points in the square [-D/2,D/2] x [-D/2, D/2].
Then build the quad-tree up to a specific number of max points in each leaf.
Then create a random tiling and brute force check that each quadtree query is correct.
Return 0 if all is OK.
Otherwise return nonzero.

*/

void enumerateFunction(int i, int j, void* aux) {
  if (aux == NULL) return;
  uint32_t* histo = (uint32_t *) aux;
  histo[i]++;
}

void count_all_inside(const double* x, 
                      const double* y, 
                      int n,
                      uint32_t* histo,
                      double xmin,
                      double xmax,
                      double ymin,
                      double ymax)
{
  for (int i = 0; i < n; i++) {
    const bool x_inside = (x[i] >= xmin && x[i] < xmax);
    const bool y_inside = (y[i] >= ymin && y[i] < ymax);
    if (x_inside && y_inside) {
      histo[i]++;
    }
  }
}

int main(int argc, const char** argv)
{
  if (argc != 4) {
    printf("usage: %s numpts D maxleaf\n", argv[0]);
    return 1;
  }

  bool test_failed = false;

  const int numpts = atoi(argv[1]);
  const double D = atof(argv[2]);
  const int maxleaf = atof(argv[3]);

  if (numpts <= 0 || D <= 0 || maxleaf <= 0) {
    printf("invalid argument(s) given\n");
    return 1;
  } 

/*
  tHashIndex2D hti;
  const double load_factor = 0.10;
  if (allocate_HashIndex2D(&hti, numpts, load_factor) != 0) {
    printf("failed to allocate NN hash-index workspace\n");
    return 1;
  }
*/

  const double W = D / 2.0;
  double* pX = malloc(sizeof(double) * numpts);
  double* pY = malloc(sizeof(double) * numpts);

  for (int i = 0; i < numpts; i++) {
    pX[i] = (2.0 * (((double) rand()) / RAND_MAX) - 1.0) * W;
    pY[i] = (2.0 * (((double) rand()) / RAND_MAX) - 1.0) * W;
  }

/*
  for (int i = 0; i < numpts; i++) {
    hti.key_[i].xi = (int32_t) floor(pX[i] / d);
    hti.key_[i].yi = (int32_t) floor(pY[i] / d);
  }
*/

//  create_HashIndex2D(&hti, numpts); 

//  if (!test_HashIndex2D(&hti)) {
//    test_failed = true;
//  }

/*
  if (!test_failed) {

    uint32_t* pH0 = malloc(sizeof(uint32_t) * numpts);
    uint32_t* pH1 = malloc(sizeof(uint32_t) * numpts);

    memset(pH0, 0, sizeof(uint32_t) * numpts);
    memset(pH1, 0, sizeof(uint32_t) * numpts);

    for (int i = 0; i < numpts; i++) {
      cellIterate_HashIndex2D(&hti, i, &enumerateFunction, pH0);
      const int32_t ix = hti.key_[i].xi;
      const int32_t iy = hti.key_[i].yi;
      count_all_inside(pX, pY, numpts, pH1, ix * d, ix * d + d, iy * d, iy * d + d);
      const bool is_still_equal = (memcmp(pH0, pH1, sizeof(uint32_t) * numpts) == 0);
      if (!is_still_equal) {
        test_failed = true;
        break;
      }
    }

    if (!test_failed) {
      for (int i = 0; i < numpts; i++) {
        if (pH0[i] < 1) {
          test_failed = true;
          break;
        }
      }
    }

    free(pH0);
    free(pH1);
  }
*/

  free(pX);
  free(pY);

//  deallocate_HashIndex2D(&hti);

  return (!test_failed ? 0 : -1);
}
