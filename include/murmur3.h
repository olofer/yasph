#ifndef __MURMUR3_H__
#define __MURMUR3_H__

/* ------------------------------------------------------------ */
/* MurmurHash3                                                  */
/* ------------------------------------------------------------ */

#define __mumu3_hashsize(n) ((uint32_t)1<<(n))
#define __mumu3_hashmask(n) (__mumu3_hashsize(n)-1)

static inline uint32_t __mumu3_rotl32(uint32_t x, int8_t r) {
  return (x << r) | (x >> (32 - r));
}

static inline uint32_t __mumu3_getblock32(const uint32_t * p, int i) {
  return p[i];
}

static inline uint32_t __mumu3_fmix32(uint32_t h) {
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}

static inline uint32_t MurmurHash3_x86_32(const void* key, int len, uint32_t seed) {
  const uint8_t* data = (const uint8_t*) key;
  const int nblocks = len / 4;
  uint32_t h1 = seed;
  const uint32_t c1 = 0xcc9e2d51;
  const uint32_t c2 = 0x1b873593;
  const uint32_t* blocks = (const uint32_t *)(data + nblocks * 4);
  for (int i = -nblocks; i; i++) {
    uint32_t k1 = __mumu3_getblock32(blocks, i);
    k1 *= c1;
    k1 = __mumu3_rotl32(k1, 15);
    k1 *= c2;
    h1 ^= k1;
    h1 = __mumu3_rotl32(h1, 13); 
    h1 = h1 * 5 + 0xe6546b64;
  }
  const uint8_t* tail = (const uint8_t*)(data + nblocks * 4);
  uint32_t k1 = 0;
  switch (len & 3) {
  case 3: k1 ^= tail[2] << 16;
  case 2: k1 ^= tail[1] << 8;
  case 1: k1 ^= tail[0];
          k1 *= c1; k1 = __mumu3_rotl32(k1, 15); k1 *= c2; h1 ^= k1;
  }
  h1 ^= len;
  h1 = __mumu3_fmix32(h1);
  return h1;
}

#endif
