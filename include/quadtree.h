#ifndef __QUADTREE_H__
#define __QUADTREE_H__

/*

Basic quadtree.
Can be used to find nearest neighbors.
In principle more flexible than the hash-based direct index.
Could be extended to a tree code for long range interactions (for example).

See "test_qtree.c" for usage and test.

*/

// TODO: missing "qti" version of QT build call

typedef struct tQTPoint {
  double x;
  double y;
} tQTPoint;

typedef struct tQTPointPayload {
  double x;
  double y;
  int index;
} tQTPointPayload;

typedef struct tQuadTree {
  tQTPoint cb;          // box center = (cx, cy)
  double hbw;           // half box width, box is a square always

  int npts;             // if leaf node: num. points in leaf, otherwise sum-total of points in children
  tQTPointPayload* p;   // NULL for non-leaf nodes; otherwise actual point data storage

  struct tQuadTree* pp; // child-nodes, one per quadrant
  struct tQuadTree* pm; // NULL if empty
  struct tQuadTree* mm;
  struct tQuadTree* mp;
} tQuadTree;

// upper borders excluded from "inside"
bool point_inside_box(double px, 
                      double py, 
                      const tQTPoint* cb, 
                      double hbw)
{
  const double cx = cb->x;
  const double cy = cb->y;
  return (px >= cx - hbw && px < cx + hbw && py >= cy - hbw && py < cy + hbw);
}

// Does box1 (center c1, half-width hw1) and box2 (c2, hw2) overlap in any way?
bool box_overlap_box(const tQTPoint* c1, 
                     double hw1, 
                     const tQTPoint* c2, 
                     double hw2) 
{
  const double right1 = c1->x + hw1;
  const double left1 = c1->x - hw1;
  const double right2 = c2->x + hw2;
  const double left2 = c2->x - hw2;
  if (right2 < left1) return false;
  if (right1 < left2) return false;
  const double top1 = c1->y + hw1;
  const double down1 = c1->y - hw1;
  const double top2 = c2->y + hw2;
  const double down2 = c2->y - hw2;
  if (top2 < down1) return false;
  if (top1 < down2) return false;
  return true;
}

// 00 = 3rd quadrant (mm)
// 01 = 2nd quadrant (pm)
// 10 = 4th quadrant (mp)
// 11 = 1st quadrant (pp)
uint8_t quadrant_index(const tQTPointPayload* p, 
                       const tQTPoint* c) 
{
  uint8_t k = 0;
  if (p->x >= c->x) k |= 0x01;
  if (p->y >= c->y) k |= 0x02;
  return k;
}

// return the number of nodes needed during the construction process (does not include root)
// point is the set of points (it is rearranged internally & accessed through tree traversal later)
// scratch is a temporary sorting area
int build_quadtree(tQuadTree* root, 
                   int maxInLeaf,
                   int maxDepth,
                   int level,
                   int np,
                   tQTPointPayload* point, 
                   tQTPointPayload* scratch,
                   int nn, 
                   tQuadTree* node) 
{
  // On entry, assume that root->cb and root->hbw have been set already... and go from there

  root->pp = NULL;
  root->pm = NULL;
  root->mm = NULL;
  root->mp = NULL;

  if (np <= maxInLeaf || level == maxDepth) {
    root->npts = np;
    root->p = point;
    return 0;
  }

  // Too many points to store as a leaf; subdivide by recursion..
  // First apply counting sort on given point array and then 
  // call quadtree builder on each required (non-empty) quadrant point set.

  root->npts = 0;
  root->p = NULL;

  if (nn == 0)
    return 0; // cannot continue since there are no more nodes to grab from the heap

  int count[4] = {0, 0, 0, 0}; // mm, pm, mp, pp

  for (int i = 0; i < np; i++) {
    const uint8_t ikey = quadrant_index(&(point[i]), &(root->cb));
    count[ikey]++;
    scratch[i] = point[i];
  }

  for (int i = 1; i <= 3; i++) {
    count[i] += count[i - 1]; // will hold offsets in the end
  }

  for (int i = np - 1; i >= 0; i--) {
    const uint8_t ikey = quadrant_index(&(scratch[i]), &(root->cb));
    count[ikey]--;
    point[count[ikey]] = scratch[i];
  }

  // point array re-ordered!
  // index: count[0] .. count[1] - 1 : 3rd
  // index: count[1] .. count[2] - 1 : 2nd
  // index: count[2] .. count[3] - 1 : 4th
  // index: count[3] .. np - 1       : 1st

  int n_3rd = count[1] - count[0];
  int n_2nd = count[2] - count[1];
  int n_4th = count[3] - count[2];
  int n_1st = np - count[3];

  // OK now for each quadrant set which is non-empty; create a child node and
  // recurse; making sure the available node store is updated ..

  const double hhbw = root->hbw / 2.0;

  int nnodes = 0;

  if (n_3rd != 0) {
    nnodes++;
    root->mm = node++;
    nn--;
    root->mm->cb.x = root->cb.x - hhbw;
    root->mm->cb.y = root->cb.y - hhbw;
    root->mm->hbw = hhbw;
    int nn_3rd = build_quadtree(root->mm, 
                                maxInLeaf,
                                maxDepth,
                                level + 1,
                                n_3rd,
                                &(point[count[0]]), 
                                &(scratch[count[0]]),
                                nn, 
                                node);
    nnodes += nn_3rd;
    nn -= nn_3rd;
    node += nn_3rd;
  }

  if (n_2nd != 0) {
    nnodes++;
    root->pm = node++;
    nn--;
    root->pm->cb.x = root->cb.x + hhbw;
    root->pm->cb.y = root->cb.y - hhbw;
    root->pm->hbw = hhbw;
    int nn_2nd = build_quadtree(root->pm, 
                                maxInLeaf,
                                maxDepth,
                                level + 1,
                                n_2nd,
                                &(point[count[1]]), 
                                &(scratch[count[1]]),
                                nn, 
                                node);
    nnodes += nn_2nd;
    nn -= nn_2nd;
    node += nn_2nd;
  }

  if (n_4th != 0) {
    nnodes++;
    root->mp = node++;
    nn--;
    root->mp->cb.x = root->cb.x - hhbw;
    root->mp->cb.y = root->cb.y + hhbw;
    root->mp->hbw = hhbw;
    int nn_4th = build_quadtree(root->mp, 
                                maxInLeaf,
                                maxDepth,
                                level + 1,
                                n_4th,
                                &(point[count[2]]), 
                                &(scratch[count[2]]),
                                nn, 
                                node);
    nnodes += nn_4th;
    nn -= nn_4th;
    node += nn_4th;
  }

  if (n_1st != 0) {
    nnodes++;
    root->pp = node++;
    nn--;
    root->pp->cb.x = root->cb.x + hhbw;
    root->pp->cb.y = root->cb.y + hhbw;
    root->pp->hbw = hhbw;
    int nn_1st = build_quadtree(root->pp, 
                                maxInLeaf,
                                maxDepth,
                                level + 1,
                                n_1st,
                                &(point[count[3]]), 
                                &(scratch[count[3]]),
                                nn, 
                                node);
    nnodes += nn_1st;
    nn -= nn_1st;
    node += nn_1st;
  }

  return nnodes;
}

// Simple box-query that just returns the total number of points found
int quadtree_box_query_count(const tQuadTree* root, 
                             const tQTPoint* cq, 
                             double hwq)
{
  if (!box_overlap_box(&(root->cb), root->hbw, cq, hwq))
    return 0;

  if (root->p != NULL) {
    int nqpts = 0;
    for (int j = 0; j < root->npts; j++) {
      if (point_inside_box(root->p[j].x, root->p[j].y, cq, hwq)) {
        nqpts++;
      }
    }
    return nqpts;
  }

  const tQuadTree* child[4] = {root->pp, root->pm, root->mm, root->mp};
  int ncq = 0;

  for (int c = 0; c < 4; c++) {
    if (child[c] != NULL)
      ncq += quadtree_box_query_count(child[c], cq, hwq);
  }

  return ncq;
}

// Count the number of nodes in the tree (supposed to match the return value of the build code)
int count_quadtree_nodes(const tQuadTree* root) {
  if (root->p != NULL)
    return 1;  // this is a leaf node
  const tQuadTree* child[4] = {root->pp, root->pm, root->mm, root->mp};
  int nn = 1;  // count this node.. plus all the nodes within the children..
  for (int c = 0; c < 4; c++) {
    if (child[c] == NULL)
      continue;
    nn += count_quadtree_nodes(child[c]);
  }
  return nn;
}

// Tally up the total number of points stored in the leaves.
int count_quadtree_points(const tQuadTree* root) {
  if (root->p != NULL)
    return root->npts;
  const tQuadTree* child[4] = {root->pp, root->pm, root->mm, root->mp};
  int n = 0;
  for (int c = 0; c < 4; c++) {
    if (child[c] == NULL)
      continue;
    n += count_quadtree_points(child[c]);
  }
  return n;
}

// Return total number of leaf nodes
int count_quadtree_leaves(const tQuadTree* root) {
  if (root->p != NULL)
    return 1;
  const tQuadTree* child[4] = {root->pp, root->pm, root->mm, root->mp};
  int n = 0;
  for (int c = 0; c < 4; c++) {
    if (child[c] == NULL)
      continue;
    n += count_quadtree_leaves(child[c]);
  }
  return n;
}

// Tally up maximum depth of tree.
int count_maximum_depth(const tQuadTree* root, 
                        int ref) 
{
  if (root->p != NULL)
    return ref;
  int d = ref;
  const tQuadTree* child[4] = {root->pp, root->pm, root->mm, root->mp};
  for (int c = 0; c < 4; c++) {
    if (child[c] == NULL)
      continue;
    const int thisd = count_maximum_depth(child[c], ref + 1);
    if (thisd > d)
      d = thisd;
  }
  return d;
}

// tally up the average depth of the tree
double count_average_depth(const tQuadTree* root, 
                           double ref) 
{
  if (root->p != NULL)
    return ref;
  const tQuadTree* child[4] = {root->pp, root->pm, root->mm, root->mp};
  int nc = 0;
  double dsum = 0.0;
  for (int c = 0; c < 4; c++) {
    if (child[c] == NULL)
      continue;
    const double thisd = count_average_depth(child[c], ref + 1.0);
    dsum += thisd;
    nc++;
  }
  return (dsum / nc);
}

// Take callback function and apply it to indices j interacting with  (iq, j).
// It is assumed that the query box centre is set to the actual point with index iq.
// And the query box half-width should be equal to the kernel support radius.

typedef void (*quadtree_interact_func_ptr)(int, int, void*);

void quadtree_box_interact(const tQuadTree* root,
                           int iq, 
                           const tQTPoint* cq, 
                           double hwq,
                           quadtree_interact_func_ptr callb_ij,
                           void* auxptr)
{
  if (!box_overlap_box(&(root->cb), root->hbw, cq, hwq))
    return;

  if (root->p != NULL) {
    for (int j = 0; j < root->npts; j++) {
      if (point_inside_box(root->p[j].x, root->p[j].y, cq, hwq)) {
        (*callb_ij)(iq, root->p[j].index, auxptr);
      }
    }
    return;
  }

  const tQuadTree* child[4] = {root->pp, root->pm, root->mm, root->mp};

  for (int c = 0; c < 4; c++) {
    if (child[c] != NULL) {
      quadtree_box_interact(child[c], iq, cq, hwq, callb_ij, auxptr);
    }
  }

  return;
}

typedef struct tQuadTreeIndex {
  int numpts;
  int maxnodes;
  tQTPointPayload* pt;
  tQTPointPayload* pt_scratch;
  tQuadTree* nodes_store;
  tQuadTree root;
} tQuadTreeIndex;

bool allocateQuadTreeIndex(tQuadTreeIndex* qti, 
                           int size) 
{
  memset(qti, 0, sizeof(tQuadTreeIndex));

  qti->numpts = size;
  qti->maxnodes = 4 * size;  // maybe 2 always works?

  qti->pt = malloc(sizeof(tQTPointPayload) * qti->numpts);
  qti->pt_scratch = malloc(sizeof(tQTPointPayload) * qti->numpts);
  qti->nodes_store = malloc(sizeof(tQuadTree) * qti->maxnodes);

  return (qti->pt != NULL && 
          qti->pt_scratch != NULL && 
          qti->nodes_store != NULL);
}

void freeQuadTreeIndex(tQuadTreeIndex* qti) {
  if (qti->pt != NULL) free(qti->pt);
  if (qti->pt_scratch != NULL) free(qti->pt_scratch);
  if (qti->nodes_store != NULL) free(qti->nodes_store);

  memset(qti, 0, sizeof(tQuadTreeIndex));
}

// qti->pt assumed to be preloaded with point data
// qti->root assumed to be preloaded with initial & correct bounding box/square
int rebuildQuadTreeIndex(tQuadTreeIndex* qti,
                         int maxInLeaf,
                         int maxDepth)
{
  const int rootLevel = 0;
  const int nno = build_quadtree(&qti->root, 
                                 maxInLeaf, 
                                 maxDepth, 
                                 rootLevel, 
                                 qti->numpts, 
                                 qti->pt, 
                                 qti->pt_scratch, 
                                 qti->maxnodes, 
                                 qti->nodes_store);
  return nno;
}

#endif 
