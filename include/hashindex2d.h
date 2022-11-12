#ifndef __HASHINDEX2D_H__
#define __HASHINDEX2D_H__

#include "murmur3.h"

/* ------------------------------------------------------------ */
/* Basic hashtable functionality                                */
/* ------------------------------------------------------------ */ 

typedef struct __mmkey {
  int32_t xi;
  int32_t yi;
} __mmkey;

const __mmkey __mmkey_unused = {INT32_MIN, INT32_MIN};

typedef struct __mmvalue {
  uint32_t index;
  uint32_t count;
} __mmvalue;

typedef struct __mmkeyvalue {
  __mmkey key;
  __mmvalue value;
} __mmkeyvalue;

const int __mmkeysize_bytes = sizeof(__mmkey);
const int __mmvaluesize_bytes = sizeof(__mmvalue); // sizeof(__mmkeyvalue) - sizeof(__mmkey);

typedef struct __mmhashtable {
  int items_;
  int hashsize_;
  int nbuckets_; // num buckets = 2^hashsize_
  int maxprobe_; // should be less than or equal to nbuckets_
  __mmkeyvalue* buckets_;
  uint32_t mask_;
  uint32_t seed_;
} __mmhashtable;

int __mmht_hashsize_from_loadfactor(int n, float lf) {
  // smallest j such that: n / 2^j <= lf
  int j = 0;
  while ((float) n / __mumu3_hashsize(j) > lf)
    j++;
  return j;
}

void __mmht_empty(__mmhashtable* ht) {
  ht->items_ = 0;
  for (int i = 0; i < ht->nbuckets_; i++)
    ht->buckets_[i].key = __mmkey_unused; 
}

// Set up an empty hash table with max items = 2^hsz
void __mmht_allocate(__mmhashtable* ht, int hsz) {
  ht->seed_ = 0xdeadbeef;
  ht->hashsize_ = hsz;
  ht->mask_ = __mumu3_hashmask(ht->hashsize_);
  ht->nbuckets_ = __mumu3_hashsize(ht->hashsize_);
  ht->maxprobe_ = (ht->nbuckets_ >> 1);
  ht->buckets_ = (__mmkeyvalue*) malloc(sizeof(__mmkeyvalue) * (ht->nbuckets_));
  __mmht_empty(ht);
  return;
}

void __mmht_deallocate(__mmhashtable* ht) {
  free(ht->buckets_);
  ht->buckets_ = NULL;
  ht->nbuckets_ = 0;
  ht->maxprobe_ = 0;
  ht->hashsize_ = 0;
  ht->items_ = 0;
  return;
}

// If returns true; then ht->buckets_[b].value, can be used to access value
// either for writing or reading.
// If returns false, then b can be used as the first available address where
// the key should've been (but only if b != ht->nbuckets_) 
bool __mmht_find_key(const __mmhashtable* ht, 
                     const __mmkey* k, 
                     uint32_t* b)
{
  uint32_t adr = MurmurHash3_x86_32((const void *) k, __mmkeysize_bytes, ht->seed_);
  adr &= ht->mask_;
  const __mmkey kunused = __mmkey_unused;
  int nprobes = 0;
  bool key_hit = false;
  while (nprobes < ht->maxprobe_) {
    key_hit = memcmp(&(ht->buckets_[adr].key), k, __mmkeysize_bytes) == 0;
    if (key_hit)
      break;
    if (memcmp(&(ht->buckets_[adr].key), &kunused, __mmkeysize_bytes) == 0)
      break;
    adr++;
    adr &= ht->mask_;
    nprobes++;
  }
  *b = ( nprobes == ht->maxprobe_ ? (uint32_t) ht->nbuckets_ : adr );
  return key_hit;
}

bool __mmht_address_ok(const __mmhashtable* ht, uint32_t b) {
  return (b < (uint32_t) ht->nbuckets_);
}

// These functions are 'unprotected', the bucket address b need to be valid
// It is the responsibility of the application to make it so.

void __mmht_get_value_at(const __mmhashtable* ht,
                         uint32_t b,
                         __mmvalue* v)
{
  const __mmvalue* pv = &(ht->buckets_[b].value);
  memcpy(v, pv, __mmvaluesize_bytes);
  return;
}

void __mmht_set_value_at(__mmhashtable* ht,
                         uint32_t b,
                         const __mmvalue* v)
{
  __mmvalue* pv = &(ht->buckets_[b].value);
  memcpy(pv, v, __mmvaluesize_bytes);
  return;
}

void __mmht_insert_at(__mmhashtable* ht,
                      uint32_t b,
                      const __mmkey* k,
                      const __mmvalue* v)
{
  __mmkeyvalue* pkv = &(ht->buckets_[b]);
  __mmkey* pk = &(pkv->key);
  memcpy(pk, k, __mmkeysize_bytes);
  __mmvalue* pv = &(pkv->value);
  memcpy(pv, v, __mmvaluesize_bytes);
  ht->items_++;
}

// use with great care
void __mmht_delete_at(__mmhashtable* ht,
                      uint32_t b)
{
  ht->buckets_[b].key = __mmkey_unused;
  ht->items_--;
}

bool __mmht_unused_at(const __mmhashtable* ht,
                      uint32_t b)
{
  const __mmkey kunused = __mmkey_unused;
  return (memcmp(&(ht->buckets_[b].key), &kunused, __mmkeysize_bytes) == 0);
}

bool __mmht_used_at(const __mmhashtable* ht,
                    uint32_t b)
{
  return !(__mmht_unused_at(ht, b));
}

/* ------------------------------------------------------------ */
/* Hash-based NN nearest-neighbor index "API"                   */
/* ------------------------------------------------------------ */

typedef struct tHashIndex2D {
  __mmhashtable ht_;
  __mmkey*  key_;
  uint32_t* kcell_;
  uint32_t* count_;
  uint32_t* cell_;
  uint32_t* index_;
  uint32_t* addr_;
  int n_;
  int c_;
} tHashIndex2D;

int allocate_HashIndex2D(tHashIndex2D* hti, int n, float LF) {
  int hsz = __mmht_hashsize_from_loadfactor(n, LF);
  __mmht_allocate(&(hti->ht_), hsz); // this will simply crash if it fails
  hti->key_   = (__mmkey*)  malloc(sizeof(__mmkey) * n);
  hti->kcell_ = (uint32_t*) malloc(sizeof(uint32_t) * n);
  hti->count_ = (uint32_t*) malloc(sizeof(uint32_t) * n);
  hti->cell_  = (uint32_t*) malloc(sizeof(uint32_t) * n);
  hti->index_ = (uint32_t*) malloc(sizeof(uint32_t) * n);
  hti->addr_ = NULL; // useful for storing bpos - facilitating special block deletion / roll back
  hti->n_ = 0;
  hti->c_ = 0;
  return 0;
}

void deallocate_HashIndex2D(tHashIndex2D* hti) {
  free(hti->key_);
  free(hti->kcell_);
  free(hti->count_);
  free(hti->cell_);
  free(hti->index_);
  __mmht_deallocate(&(hti->ht_));
  return;
}

// assume ht is cleared and has space allocated for the problem at hand
// there will be at most n keys 0..n-1

// the hti->key_ array need to be pre-filled correctly {floor(x/r), floor(y/r)}
// for the elements 0...n-1
int create_HashIndex2D(tHashIndex2D* hti, int n)
{
  __mmhashtable* ht = &(hti->ht_);
  __mmht_empty(ht);

  __mmkey theKey;
  __mmvalue theValue;
  uint32_t bpos;
  int unq = 0;

  for (int i = 0; i < n; i++) {
    memcpy(&theKey, &(hti->key_[i]), sizeof(__mmkey));
    bool haskey = __mmht_find_key(ht, &theKey, &bpos);
    if (!haskey && bpos < (uint32_t) ht->nbuckets_) {
      theValue.index = unq;
      __mmht_insert_at(ht, bpos, &theKey, &theValue);
      hti->count_[unq++] = 1;
    } else {
      __mmht_get_value_at(ht, bpos, &theValue);
      hti->count_[theValue.index]++;
    }
    hti->kcell_[i] = theValue.index;
  }

  for (int i = 1; i < unq; i++)
    hti->count_[i] += hti->count_[i - 1];

  for (int i = n - 1; i >= 0; i--) {
    uint32_t j = hti->kcell_[i];
    hti->count_[j]--;
    hti->cell_[hti->count_[j]] = j;
    hti->index_[hti->count_[j]] = i;
  }

  hti->n_ = n;
  hti->c_ = unq;

  return unq;
}

bool test_HashIndex2D(const tHashIndex2D* hti) {
  for (int i = 0; i < hti->n_; i++)
    if (hti->cell_[i] != hti->kcell_[hti->index_[i]])
      return false;
  for (int i = 0; i < hti->c_; i++)
    if ((uint32_t) i != hti->cell_[hti->count_[i]])
      return false;
  return true;
}

typedef void (*HashIndex2D_func_ptr)(int, int, void*);

// Run callback with arguments (j, j, aux) for every point j 
// in the same cell as query point i; the callbacks include the query j = i
void cellIterate_HashIndex2D(const tHashIndex2D* hti,
                             int i,
                             HashIndex2D_func_ptr callb_ij,
                             void* auxptr)
{
  const int npts = hti->n_;
  const uint32_t* cell = hti->cell_;
  const uint32_t* count = hti->count_;
  const uint32_t* index = hti->index_;
  const __mmhashtable* ht = &(hti->ht_);

  __mmkey Kq;
  memcpy(&Kq, &(hti->key_[i]), sizeof(__mmkey));

  uint32_t adr = 0;
  if (!__mmht_find_key(ht, &Kq, &adr)) return;

  __mmvalue Vq;
  __mmht_get_value_at(ht, adr, &Vq); 
  const int c = Vq.index;
  int ipos = count[c];
  while (ipos < npts && cell[ipos] == (uint32_t)c) {
    const int j = index[ipos];
    (*callb_ij)(j, j, auxptr);
    ipos++;
  }

  return;
}

// Run callback on all NN candidates to index i (includes self)
void singleInteract_HashIndex2D(const tHashIndex2D* hti,
                                int i,
                                HashIndex2D_func_ptr callb_ij,
                                void* auxptr)
{
  const int npts = hti->n_;
  const uint32_t* cell = hti->cell_;
  const uint32_t* count = hti->count_;
  const uint32_t* index = hti->index_;
  const __mmhashtable* ht = &(hti->ht_);
  
  const int32_t key_perturbation[2 * 9] =
    { 0, 0,
      1, 0,
      1, 1,
      0, 1,
     -1, 1,
     -1, 0,
     -1,-1,
      0,-1,
      1,-1 };

  __mmkey Kq;
  memcpy(&Kq, &(hti->key_[i]), sizeof(__mmkey));

  for (int p = 0; p < 9; p++) {
    __mmkey Kp = {Kq.xi + key_perturbation[2 * p], 
                  Kq.yi + key_perturbation[2 * p + 1]};
    uint32_t adr = 0;
    if (!__mmht_find_key(ht, &Kp, &adr)) continue;
    __mmvalue Vq;
    __mmht_get_value_at(ht, adr, &Vq); 
    const int c = Vq.index;
    int ipos = count[c];
    while (ipos < npts && cell[ipos] == (uint32_t)c) {
      const int j = index[ipos];
      (*callb_ij)(i, j, auxptr);
      ipos++;
    }
  }

  return;
}

// Run callback on all indices (i, j) where i is all indices in cell c
// and j includes all NN candidates (including all indices in cell c)
void cellInteract_HashIndex2D(const tHashIndex2D* hti,
                              int c,
                              HashIndex2D_func_ptr callb_ij,
                              void* auxptr)
{
  const int npts = hti->n_;
  const uint32_t* cell = hti->cell_;
  const uint32_t* count = hti->count_;
  const uint32_t* index = hti->index_;
  const __mmhashtable* ht = &(hti->ht_);

  const int32_t key_perturbation[2 * 9] =
    { 0, 0,
      1, 0,
      1, 1,
      0, 1,
     -1, 1,
     -1, 0,
     -1,-1,
      0,-1,
      1,-1 };

  int cqvec[9];

  if (c >= hti->c_ || c < 0) return;

  // K <- key[index[count[c]]]
  // then perturb that "central" key

  __mmkey Kq;
  memcpy(&Kq, &(hti->key_[index[count[c]]]), sizeof(__mmkey));

  for (int p = 0; p < 9; p++) {
    const __mmkey Kp = { Kq.xi + key_perturbation[2 * p], 
                         Kq.yi + key_perturbation[2 * p + 1] };
    uint32_t adr = 0;
    if (!__mmht_find_key(ht, &Kp, &adr)) {
      cqvec[p] = -1;
    } else {
      __mmvalue Vq;
      __mmht_get_value_at(ht, adr, &Vq);
      cqvec[p] = Vq.index;
    }
  }

  // Done with querying the hash table; now use cqvec[] to step through the 
  // possible pair interactions

  const int c0 = cqvec[0];
  int ipos = count[c0];
  while (ipos < npts && cell[ipos] == (uint32_t)c0) {
    const int iidx = index[ipos];

    //const double xi = x[iidx];
    //const double yi = y[iidx];
    //uint32_t nirq = 0;
    //double sumirq = 0.0;

    // this point should now interact with the complete listing under cqvec[0..8]
    for (int j = 0; j < 9; j++) {
      const int cj = cqvec[j];
      int jpos = count[cj];
      while (jpos < npts && cell[jpos] == (uint32_t)cj) {
        const int jidx = index[jpos];

        /*
         * Here the index iidx (cell c) interacts with jidx (cell c or its NNs)
         */

        (*callb_ij)(iidx, jidx, auxptr);

        /*const double xj = x[jidx];
        const double yj = y[jidx];
        const double dxij = xj - xi;
        const double dyij = yj - yi;
        const double r2ij = dxij * dxij + dyij * dyij;
        if (r2ij <= R2) {
          nirq++;
          sumirq += sqrt(r2ij);
        }*/

        jpos++;
      }
    }

    // store results
    //nrq[iidx] = nirq;
    //sumrq[iidx] = sumirq;

    ipos++;
  }

  return;
}

#endif
