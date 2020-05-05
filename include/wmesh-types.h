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

typedef wmesh_int_t * __restrict__ wmesh_int_p;
typedef const wmesh_int_t * __restrict__ const_wmesh_int_p;

typedef char wmesh_str_t[256];

#include <stdio.h>
#ifdef __cplusplus
extern "C"
{
#endif
  struct wmesh_int_mat_t;
  struct wmesh_int_sparsemat_t;
  
  void wmesh_int_mat_def	(wmesh_int_mat_t * self_,
				 wmesh_int_t m,
				 wmesh_int_t n,
				 wmesh_int_p v,
				 wmesh_int_t ld);
  
  void wmesh_int_sparsemat_def	(wmesh_int_sparsemat_t*self_,
				 wmesh_int_t 	size_);
  
  void wmesh_int_sparsemat_get	(wmesh_int_sparsemat_t* self_,
				 wmesh_int_t 	idx_,
				 wmesh_int_mat_t* 	wint_mat_);
  
  void wmesh_int_sparsemat_set	(wmesh_int_sparsemat_t* 	self_,
				 const_wmesh_int_p  	m_,
				 const_wmesh_int_p  	n_,
				 const_wmesh_int_p  	ld_,
				 const_wmesh_int_p  ptr_,
				 wmesh_int_p  data_);
  
  void 		wmesh_int_sparsemat_free		(wmesh_int_sparsemat_t*self_);    
  void 		wmesh_int_sparsemat_info		(const wmesh_int_sparsemat_t* self_,
					 FILE * out_);
  
  wmesh_int_t 	wmesh_int_sparsemat_fprintf(const wmesh_int_sparsemat_t* 		self_,
					  FILE * out_);


#ifdef __cplusplus
}
#endif


#endif

