#pragma once

#include "wmesh-types.h"
template<typename T>
static inline int wmesh_qsort_increasing_predicate(const void * a_,
						   const void * b_)
{
  const T*__restrict__ a = (const T*__restrict__ )a_;
  const T*__restrict__ b = (const T*__restrict__ )b_;
  if (*a < *b)
    {
      return -1;
    }
  else if (*a > *b)    
    {
      return 1;
    }
  else
    {
      return 0;
    }
}
