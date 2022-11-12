#ifndef __ARGSIO_H__
#define __ARGSIO_H__

// Management of settings strings of the form: "key1=val1 key2=val2 ...".
// Supports initialization from (argc, argv), can set and get values, can read/write to file.
// Has utilities for getting/setting integers and reals.
// Not intended for large (key, val) sets (inefficient linear probing of plain string).

// TODO: "protected" get_value (buffer size limit)
// TODO: lookup by index, delete by index

#define ARGSIO_LOCAL_BUFFER 256
#define ARGSIO_READ_BUFFER 1024

typedef struct tArgsio {
  char* text;
  int textsize;
} tArgsio;

void argsio_clear(tArgsio* A)
{
  memset(A->text, 0, sizeof(char) * A->textsize);
}

void argsio_init(tArgsio* A, 
                 int size) 
{
  A->text = malloc(sizeof(char) * size);
  A->textsize = size;
  argsio_clear(A);
}

int argsio_length(const tArgsio* A)
{
  return strlen(A->text);
}

int argsio_bufferlength(const tArgsio* A)
{
  return A->textsize;
}

void argsio_uninit(tArgsio* A)
{
  if (A->text != NULL) free(A->text);
  A->text = NULL;
  A->textsize = 0;
}

void argsio_extend(tArgsio* A)
{
  int textsize_ = A->textsize * 2;
  char* text_ = malloc(sizeof(char) * textsize_);
  memcpy(text_, A->text, sizeof(char) * A->textsize);
  argsio_uninit(A);
  A->text = text_;
  A->textsize = textsize_;
}

void argsio_shrink_to_fit(tArgsio* A)
{
  int textsize_ = strlen(A->text) + 1;
  char* text_ = malloc(sizeof(char) * textsize_);
  memcpy(text_, A->text, sizeof(char) * textsize_);
  argsio_uninit(A);
  A->text = text_;
  A->textsize = textsize_;
}

int argsio_count(const tArgsio* A)
{
  int n = 0;
  const char* kv = A->text;
  while (kv != NULL) {
    const char* pch = strchr(kv, '=');
    if (pch == NULL) break;
    n++;
    kv = strchr(pch, ' ');
    if (kv != NULL) kv++;
  }
  return n;
}

const char* argsio_ith_kv(const tArgsio* A, 
                          int i)
{
  const int count = argsio_count(A);
  if (i < 0 || i >= count) return NULL;
  const char* kv = A->text;
  int n = -1;
  while (++n != i) {
    kv = strchr(kv, ' ');
    kv++;
  }
  return kv;
}

bool argsio_util_kv_split(const char* kv,
                          char* key_dest,
                          char* val_dest)
{
  const char* pch = strchr(kv, '=');  
  if (pch == NULL) return false;
  if (key_dest != NULL) {
    const int len = pch - kv;
    memcpy(key_dest, kv, sizeof(char) * len);
    key_dest[len] = '\0';
  }
  if (val_dest != NULL) {
    pch++;
    const int len = &kv[strlen(kv)] - pch;
    memcpy(val_dest, pch, sizeof(char) * len);
    val_dest[len] = '\0';
  }
  return true;
}

// string s has this format: "r1,r2,r3,..."
int argsio_util_parse_vector(const char* s,
                             int maxv, 
                             double* v)
{
  int nv = 0;
  if (s == NULL || v == NULL || maxv <= 0) return nv;
  const char* beg = s;
  const char* pch = strchr(beg, ',');
  for (;;) {
    const int togo = strlen(beg);
    if (pch == NULL && togo != 0) pch = &beg[togo];
    if (pch == NULL || nv == maxv) break;
    const int len = pch - beg;
    char textvalue[len + 1];
    memcpy(textvalue, beg, sizeof(char) * len);
    textvalue[len] = '\0';
    //printf("[%i]: \"%s\"\n", nv, textvalue);
    v[nv++] = atof(textvalue);
    if (*pch == '\0') break;
    beg = pch + 1;
    pch = strchr(beg, ',');
  }
  return nv;
}

bool argsio_get_key_from_index(const tArgsio* A,
                               int i,
                               char* key_dest)
{
  const char* kv = argsio_ith_kv(A, i);
  if (kv == NULL) return false;
  const char* pch = strchr(kv, '=');
  const int len = pch - kv;
  if (key_dest != NULL) {
    memcpy(key_dest, kv, sizeof(char) * len);
    key_dest[len] = '\0';
  }
  return true;
}

bool argsio_get_val_from_index(const tArgsio* A,
                               int i,
                               char* val_dest)
{
  const char* kv = argsio_ith_kv(A, i);
  if (kv == NULL) return false;
  const char* pch = strchr(kv, '=') + 1;
  const char* pend = strchr(kv, ' ');
  const int len = (pend != NULL ? (pend - pch) : (A->text + strlen(A->text) - pch));
  if (val_dest != NULL) {
    memcpy(val_dest, pch, sizeof(char) * len);
    val_dest[len] = '\0';
  }
  return true;
}

const char* argsio_lookup(const tArgsio* A, 
                          const char* key) 
{
  const int klen = strlen(key);
  const char* pval = NULL;
  const char* kv = A->text;
  while (true) {
    const char* pch = strchr(kv, '=');
    if (pch == NULL) break; // if A is empty
    const bool length_match = ((int)(pch - kv) == klen);
    if (length_match && memcmp(kv, key, sizeof(char) * klen) == 0) {
      pval = pch + 1;
      return pval;
    }
    kv = strchr(pch, ' ');
    if (kv == NULL) break;
    kv += 1;
  }
  return pval;
}

int argsio_count_matches(const tArgsio* A,
                         const char* key)
{
  int m = 0;
  const int klen = strlen(key);
  const char* kv = A->text;
  while (true) {
    const char* pch = strchr(kv, '=');
    if (pch == NULL) break;
    const bool length_match = ((int)(pch - kv) == klen);
    if (length_match && memcmp(kv, key, sizeof(char) * klen) == 0) m++;
    kv = strchr(pch, ' ');
    if (kv == NULL) break;
    kv += 1;
  }
  return m;
}

bool argsio_is_unique(const tArgsio* A,
                      const char* key)
{
  return (argsio_count_matches(A, key) == 1);
}

int argsio_count_prefix_matches(const tArgsio* A,
                                const char* key_prefix)
{
  int m = 0;
  const int klen_prefix = strlen(key_prefix);
  const char* kv = A->text;
  while (true) {
    const char* pch = strchr(kv, '=');
    if (pch == NULL) break;
    const bool length_surpass_or_equal = ((int)(pch - kv) >= klen_prefix);
    if (length_surpass_or_equal && 
        memcmp(kv, key_prefix, sizeof(char) * klen_prefix) == 0) m++;
    kv = strchr(pch, ' ');
    if (kv == NULL) break;
    kv += 1;
  }
  return m;
}

const char* argsio_prefix_lookup(const tArgsio* A,
                                 const char* key_prefix,
                                 int subcount)
{
  int m = 0;
  const char* pval = NULL;
  const int klen_prefix = strlen(key_prefix);
  const char* kv = A->text;
  while (true) {
    const char* pch = strchr(kv, '=');
    if (pch == NULL) break;
    const bool length_surpass_or_equal = ((int)(pch - kv) >= klen_prefix);
    if (length_surpass_or_equal && 
        memcmp(kv, key_prefix, sizeof(char) * klen_prefix) == 0)
    {
      if (m == subcount) {
        pval = pch + 1;
        break;
      }
      m++;
    }
    kv = strchr(pch, ' ');
    if (kv == NULL) break;
    kv += 1;
  }
  return pval;
}

bool argsio_all_unique(const tArgsio* A)
{
  const int c = argsio_count(A);
  for (int i = 0; i < c; i++) {
    char ikey[ARGSIO_LOCAL_BUFFER];
    argsio_get_key_from_index(A, i, ikey);
    if (!argsio_is_unique(A, ikey))
      return false;
  }
  return true;
}

bool argsio_get_value(const tArgsio* A,
                      const char* key,
                      char* val_dest)
{
  const char* pval = argsio_lookup(A, key);
  if (pval == NULL)
    return false;
  const char* pend = strchr(pval, ' ');
  if (pend == NULL)
    pend = A->text + strlen(A->text) + 1;
  const int val_len = (int) (pend - pval);
  memcpy(val_dest, pval, sizeof(char) * val_len);
  val_dest[val_len] = '\0';
  return true;
}

bool argsio_get_prefix_value(const tArgsio* A,
                             const char* key_prefix,
                             int subcount,
                             char* val_dest)
{
  const char* pval = argsio_prefix_lookup(A, key_prefix, subcount);
  if (pval == NULL)
    return false;
  const char* pend = strchr(pval, ' ');
  if (pend == NULL)
    pend = A->text + strlen(A->text) + 1;
  const int val_len = (int) (pend - pval);
  memcpy(val_dest, pval, sizeof(char) * val_len);
  val_dest[val_len] = '\0';
  return true;
}

bool argsio_get_int(const tArgsio* A, 
                    const char* key,
                    int* val)
{
  char strval[ARGSIO_LOCAL_BUFFER];
  const bool haskey = argsio_get_value(A, key, strval);
  if (haskey && val != NULL) {
    *val = atoi(strval);
  }
  return haskey;
}

bool argsio_get_real(const tArgsio* A, 
                     const char* key,
                     double* val)
{
  char strval[ARGSIO_LOCAL_BUFFER];
  const bool haskey = argsio_get_value(A, key, strval);
  if (haskey && val != NULL) {
    *val = atof(strval);
  }
  return haskey;
}

bool argsio_get_int_via_cast(const tArgsio* A, 
                             const char* key,
                             int* val)
{
  double tmpreal;
  const bool haskey = argsio_get_real(A, key, &tmpreal);
  if (haskey && val != NULL) {
    *val = (int) tmpreal;
  }
  return haskey;
}

// Add a string of the form "key=val" (assumed to be in that form)
bool argsio_add_kv(tArgsio* A,
                   const char* keyval) 
{
  const int len = strlen(A->text);
  const int kvlen = strlen(keyval);
  if (len + kvlen >= A->textsize)
    argsio_extend(A);
  if (len != 0)
    strcat(A->text, " ");
  const char* ret = strcat(A->text, keyval);
  return (ret == A->text);
}

bool argsio_add(tArgsio* A,
                const char* key, 
                const char* val) 
{
  if (strchr(key, ' ') != NULL || strchr(val, ' ') != NULL)
    return false;
  const int klen = strlen(key);
  const int vlen = strlen(val);
  char buf[klen + vlen + 2];
  buf[0] = '\0';
  strcat(buf, key);
  strcat(buf, "=");
  strcat(buf, val);
  return argsio_add_kv(A, buf);
}

bool argsio_remove_key(tArgsio* A,
                       const char* key)
{
  const char* pval = argsio_lookup(A, key);
  if (pval == NULL)
    return false;
  char* pstart = (char *) pval;
  while (pstart != A->text && *pstart != ' ') pstart--;
  char* pend = strchr(pval, ' ');
  if (pend == NULL) {
    *pstart = '\0';
    return true;
  }
  if (*pstart == ' ') pstart++;
  char* dest = pstart;
  char* src = pend + 1;
  while (true) {
    *dest = *src;
    if (*src == '\0')
      break;
    dest++;
    src++;
  }
  return true;
}

bool argsio_set_value(tArgsio* A,
                      const char* key,
                      const char* val_src,
                      bool add_if_missing,
                      bool overwrite_if_not_missing)
{
  const char* pval = argsio_lookup(A, key);
  const bool is_missing = (pval == NULL);
  if (is_missing && add_if_missing) {
    return argsio_add(A, key, val_src);
  }
  if (!is_missing && overwrite_if_not_missing) {
    if (!argsio_remove_key(A, key)) 
      return false;
    return argsio_add(A, key, val_src);
  }
  return true;
}

const char* argsio_is_keyval_string(const char* kv)
{
  const int kvlen = strlen(kv);
  if (kvlen == 0)
    return NULL;
  if (strchr(kv, ' ') != NULL)
    return NULL;
  const char* pch = strchr(kv, '=');
  if (pch == NULL)
    return NULL;
  const char* prch = strrchr(kv, '=');
  if (pch != prch)
    return NULL;
  if (pch == &kv[0] || pch == &kv[kvlen - 1])
    return NULL;
  return pch;
}

int argsio_add_cmdargs(tArgsio* A, 
                       int argc, 
                       const char** argv,
                       bool append)
{
  int nadded = 0;
  for (int i = 0; i < argc; i++) {
    const char* kv = argv[i];
    if (argsio_is_keyval_string(kv) == NULL)
      continue;
    if (append) {
      if (argsio_add_kv(A, kv))
        nadded++;
    } else {
      const int kvlen = strlen(kv);
      char keybuf[kvlen];
      char valbuf[kvlen];
      argsio_util_kv_split(kv, keybuf, valbuf);
      // printf("\"%s\":\"%s\"\n", keybuf, valbuf);
      if (argsio_set_value(A, keybuf, valbuf, true, true))
        nadded++;
    }
  }
  return nadded;
}

void argsio_printf(const tArgsio* A,
                   const char* header)
{
  if (header != NULL)
    printf("=== %s ===\n", header);
  const int c = argsio_count(A);
  for (int i = 0; i < c; i++) {
    char kbuf[ARGSIO_LOCAL_BUFFER];
    char vbuf[ARGSIO_LOCAL_BUFFER];
    argsio_get_key_from_index(A, i, kbuf);
    argsio_get_val_from_index(A, i, vbuf);
    printf("[%i]: \"%s\" = \"%s\"\n", i, kbuf, vbuf);
  }
}

bool argsio_export_file(const tArgsio* A,
                        const char* filename,
                        const char* header)
{
  FILE *pfile = NULL;
	pfile = fopen(filename, "w");
	if (!pfile) return false;
  if (header != NULL)
    fprintf(pfile, "=== %s ===\n", header);
  const int c = argsio_count(A);
  for (int i = 0; i < c; i++) {
    char kbuf[ARGSIO_LOCAL_BUFFER];
    char vbuf[ARGSIO_LOCAL_BUFFER];
    argsio_get_key_from_index(A, i, kbuf);
    argsio_get_val_from_index(A, i, vbuf);
    fprintf(pfile, "%s=%s\n", kbuf, vbuf);
  }
  fclose(pfile);
  return true;
}

int argsio_import_file(tArgsio* A,
                       const char* filename)
{
  char readbuffer[ARGSIO_READ_BUFFER];

  int added = 0;
  FILE* fp = fopen(filename, "r");
	if (!fp) return 0;

  while (true) {
    char* tmp = fgets(readbuffer, ARGSIO_READ_BUFFER, fp);
		if (tmp == NULL) break;
    const int len = strlen(tmp);
    if (tmp[len - 1] == '\n') tmp[len - 1] = '\0';
    if (argsio_is_keyval_string(tmp) == NULL) continue;
    if (argsio_add_kv(A, tmp)) added++;
  }

  fclose(fp);
  return added;
}

int argsio_insert(tArgsio* A, 
                  const tArgsio* B) 
{
  int inserted = 0;
  const int cB = argsio_count(B); 
  for (int i = 0; i < cB; i++) {
    char kbuf[ARGSIO_LOCAL_BUFFER];
    char vbuf[ARGSIO_LOCAL_BUFFER];
    argsio_get_key_from_index(B, i, kbuf);
    argsio_get_val_from_index(B, i, vbuf);
    if (argsio_add(A, kbuf, vbuf)) inserted++;
  }
  return inserted;
}

bool argsio_copy(tArgsio* A,
                 const tArgsio* B)
{
  argsio_clear(A);
  return (argsio_insert(A, B) == argsio_count(A));
}

int argsio_strcmp(const tArgsio* A,
                  const tArgsio* B)
{
  return strcmp(A->text, B->text);
}

#endif
