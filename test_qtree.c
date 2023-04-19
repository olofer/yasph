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
  histo[j]++;
}

void count_all_inside(const tQTPointPayload* pt, 
                      int n,
                      uint32_t* histo,
                      double xmin,
                      double xmax,
                      double ymin,
                      double ymax)
{
  for (int i = 0; i < n; i++) {
    const double xi = pt[i].x;
    const double yi = pt[i].y;
    const bool x_inside = (xi >= xmin && xi < xmax);
    const bool y_inside = (yi >= ymin && yi < ymax);
    if (x_inside && y_inside) {
      histo[i]++;
    }
  }
}

int main(int argc, const char** argv)
{
  if (argc != 5) {
    printf("usage: %s numpts D qhw maxleaf\n", argv[0]);
    return 1;
  }

  bool test_failed = false;

  const int numpts = atoi(argv[1]);
  const double D = atof(argv[2]);
  const double qhw = atof(argv[3]);
  const int maxleaf = atoi(argv[4]);

  if (numpts <= 0 || D <= 0.0 || maxleaf <= 0 || qhw <= 0.0) {
    printf("invalid argument(s) given\n");
    return 1;
  } 

  const int maxNumberOfNodes = numpts * 4; // probably only 2 needed here

  tQTPointPayload* pt = malloc(sizeof(tQTPointPayload) * numpts);
  tQTPointPayload* pt_scratch = malloc(sizeof(tQTPointPayload) * numpts);
  tQuadTree* nodes_store = malloc(sizeof(tQuadTree) * maxNumberOfNodes);

  if (pt == NULL || pt_scratch == NULL || nodes_store == NULL) {
    printf("memory allocation failed\n");
    return 1;
  }

  const double W = D / 2.0;

  printf("num pts, leaf max = %i, %i\n", numpts, maxleaf);
  printf("domain, query half-widths = %f, %f\n", W, qhw);

  double xmin = W;
  double xmax = -1.0 * W;
  double ymin = W;
  double ymax = -1.0 * W;

  for (int i = 0; i < numpts; i++) {
    pt[i].x = (2.0 * (((double) rand()) / RAND_MAX) - 1.0) * W;
    pt[i].y = (2.0 * (((double) rand()) / RAND_MAX) - 1.0) * W;
    pt[i].index = i;

    xmin = (pt[i].x < xmin ? pt[i].x : xmin);
    xmax = (pt[i].x > xmax ? pt[i].x : xmax);
    ymin = (pt[i].y < ymin ? pt[i].y : ymin);
    ymax = (pt[i].y > ymax ? pt[i].y : ymax);
  }

  const double hbwx = (xmax - xmin) / 2.0;
  const double hbwy = (ymax - ymin) / 2.0;

  tQuadTree rootNode;

  rootNode.cb.x = (xmin + xmax) / 2.0;
  rootNode.cb.y = (ymin + ymax) / 2.0;
  rootNode.hbw = (hbwx > hbwy ? hbwx : hbwy);
  rootNode.hbw *= (1.0 + 1.0e-10);

  printf("cx, cy     = %f, %f\n", rootNode.cb.x, rootNode.cb.y);
  printf("hbwx, hbwy = %f, %f\n", hbwx, hbwy);

  const int maxDepth = 50;
  const int maxInLeaf = maxleaf;
  const int rootLevel = 0;

  const int nno = build_quadtree(&rootNode, 
                                 maxInLeaf, 
                                 maxDepth, 
                                 rootLevel, 
                                 numpts, 
                                 pt, 
                                 pt_scratch, 
                                 maxNumberOfNodes, 
                                 nodes_store);

  const int nn_in_tree = count_quadtree_nodes(&rootNode);
  const int n_in_tree = count_quadtree_points(&rootNode);
  const int max_depth = count_maximum_depth(&rootNode, 0);
  const double avg_depth = count_average_depth(&rootNode, 0.0);

  printf("nno, nn_in_tree, npts, n_in_tree = %i, %i, %i, %i\n", nno, nn_in_tree, numpts, n_in_tree);
  printf("maxdepth, <depth> = %i, %f\n", max_depth, avg_depth);

  const int total_box_count = quadtree_box_query_count(&rootNode, &(rootNode.cb), rootNode.hbw);

  if (numpts != n_in_tree || nno + 1 != nn_in_tree || total_box_count != numpts) {
    printf("point and/or node counts failed\n");
    test_failed = true;
    goto free_and_exit;
  }

  // Restore original ordering of pt array into now-unused pt_scratch area
  for (int i = 0; i < numpts; i++) {
    pt_scratch[pt[i].index] = pt[i];
  }

  for (int i = 0; i < numpts; i++) {
    if (pt_scratch[i].index != i) {
      printf("original order recreation failure\n");
      test_failed = true;
      goto free_and_exit;
    }
  }

  // Setup a interaction callback test (check all points that should be touched, are actually touched, and no other)
  uint32_t* pH0 = malloc(sizeof(uint32_t) * numpts);
  uint32_t* pH1 = malloc(sizeof(uint32_t) * numpts);

  memset(pH0, 0, sizeof(uint32_t) * numpts);
  memset(pH1, 0, sizeof(uint32_t) * numpts);

  const double query_box_halfwidth = qhw;

  for (int i = 0; i < numpts; i++)
  {
    const tQTPoint query_pt_i = {pt_scratch[i].x, pt_scratch[i].y};
    const int index_query_i = i;

    quadtree_box_interact(&rootNode,
                          index_query_i, 
                          &query_pt_i, 
                          query_box_halfwidth,
                          &enumerateFunction,
                          (void *) pH0);

    count_all_inside(pt_scratch, 
                     numpts, 
                     pH1, 
                     query_pt_i.x - query_box_halfwidth, 
                     query_pt_i.x + query_box_halfwidth,
                     query_pt_i.y - query_box_halfwidth, 
                     query_pt_i.y + query_box_halfwidth);

    const bool is_still_equal = (memcmp(pH0, pH1, sizeof(uint32_t) * numpts) == 0);
    if (!is_still_equal) {
      printf("interaction failed @ i = %i\n", i);
      test_failed = true;
      break;
    }
  }

  if (!test_failed) {
    uint32_t hmin = 1000000000;
    uint32_t hmax = 0;
    for (int i = 0; i < numpts; i++) {
      if (pH0[i] < 1) {
        printf("never-self interaction for i = %i\n", i);
        test_failed = true;
        break;
      }
      if (pH0[i] < hmin) hmin = pH0[i];
      if (pH0[i] > hmax) hmax = pH0[i];
    }
    printf("histogram min,max = %i,%i\n", hmin, hmax);
  }  

  free(pH0);
  free(pH1);

free_and_exit:

  free(pt);
  free(pt_scratch);
  free(nodes_store);

  return (!test_failed ? 0 : -1);
}
