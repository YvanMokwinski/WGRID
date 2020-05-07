#ifndef WMESH_TYPES_H
#define WMESH_TYPES_H

#ifndef WMESH_ILP64
#error WMESH_ILP64 is undefined.
#endif

#if WMESH_ILP64
  typedef long long int wmesh_int_t;
  #define WMESH_INT_FORMAT "%Ld"
#else
  typedef int wmesh_int_t;
  #define WMESH_INT_FORMAT "%d"
#endif

typedef wmesh_int_t * __restrict__ 	wmesh_int_p;
typedef const wmesh_int_t * __restrict__ const_wmesh_int_p;

typedef char wmesh_str_t[256];

#endif

