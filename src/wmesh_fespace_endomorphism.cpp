#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"
#include <chrono>
#include <iostream>
#include "bms.h"

#define WINT_SPARSE_MAT_PARAMS(_f)		\
  _f.m_ptr,					\
    _f.m_m,					\
    _f.m_n,					\
    _f.m_data,					\
    _f.m_ld

using namespace std::chrono;
extern "C"
{
  static int sort_predicate(const void * a_,
			    const void * b_)
  {
    const_wmesh_int_p a = (const_wmesh_int_p)a_;
    const_wmesh_int_p b = (const_wmesh_int_p)b_;
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
  
  wmesh_status_t wmesh_fespace_endomorphism	(const wmesh_t*__restrict__	mesh_,
						 wmesh_int_t 	degree_,
						 wmesh_int_p 	csr_size_,
						 wmesh_int_p*__restrict__ 	csr_ptr_,
						 wmesh_int_p*__restrict__ 	csr_ind_)
  {

    wmesh_status_t status;
    
    wmesh_int_t n2c_m 	= mesh_->m_num_nodes;
    wmesh_int_p n2c_ptr = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*(n2c_m+1)); 
    wmesh_int_p n2c_v 	= (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mesh_->m_c2n.m_ptr[mesh_->m_c2n.m_size]); 

    status = bms_n2c(mesh_->m_c2n.m_size,
		     WINT_SPARSE_MAT_PARAMS(mesh_->m_c2n),
		     n2c_ptr,
		     n2c_m,
		     n2c_v);
    
    WMESH_STATUS_CHECK(status);

    wmesh_int_p blank 	= (wmesh_int_p)calloc(n2c_m,sizeof(wmesh_int_t)); 
    wmesh_int_p select 	= (wmesh_int_p)calloc(n2c_m,sizeof(wmesh_int_t));
    csr_size_[0] = n2c_m;
    csr_ptr_[0] = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*(n2c_m + 1));
    csr_ptr_[0][0] = 0;
    for (wmesh_int_t idof=0;idof<n2c_m;++idof)
      {
	wmesh_int_t select_n = 0;
	//
	// Union 
	//
	for (wmesh_int_t s = n2c_ptr[idof];s<n2c_ptr[idof+1];++s)
	  {
	    wmesh_int_t c = n2c_v[s];
	    wmesh_int_t cindex, ctype;
	    
	    status = bms_n2c_cindex(c,
				    &cindex);
	    WMESH_STATUS_CHECK(status);
	    
	    status = bms_n2c_ctype	(c,
					 &ctype);
	    WMESH_STATUS_CHECK(status);
	    
	    for (wmesh_int_t i=0;i<mesh_->m_c2n.m_m[ctype];++i)
	      {
		wmesh_int_t jdof = mesh_->m_c2n.m_data[mesh_->m_c2n.m_ptr[ctype] + cindex * mesh_->m_c2n.m_ld[ctype] + i];
		if (0==blank[jdof])
		  {
		    select[select_n++] = jdof;
		    blank[jdof] = select_n;
		  }
	      }
	  }

	//
	// Reset.
	//
	for (wmesh_int_t s = 0;s<select_n;++s)
	  {
	    blank[select[s]] = 0;
	  }

	csr_ptr_[0][idof+1] = csr_ptr_[0][idof] + select_n;	
      }

    csr_ind_[0] = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*csr_ptr_[0][n2c_m]);
    for (wmesh_int_t idof=0;idof<n2c_m;++idof)
      {
	wmesh_int_t select_n = 0;
	//
	// Union 
	//
	for (wmesh_int_t s = n2c_ptr[idof];s<n2c_ptr[idof+1];++s)
	  {
	    wmesh_int_t c = n2c_v[s];
	    wmesh_int_t cindex, ctype;
	    
	    status = bms_n2c_cindex	(c,
					 &cindex);
	    WMESH_STATUS_CHECK(status);
	    
	    status = bms_n2c_ctype	(c,
					 &ctype);
	    WMESH_STATUS_CHECK(status);
	    
	    for (wmesh_int_t i=0;i<mesh_->m_c2n.m_m[ctype];++i)
	      {
		wmesh_int_t jdof = mesh_->m_c2n.m_data[mesh_->m_c2n.m_ptr[ctype] + cindex * mesh_->m_c2n.m_ld[ctype] + i];
		if (0==blank[jdof])
		  {
		    select[select_n++] = jdof;
		    blank[jdof] = select_n;
		  }
	      }
	  }

	//
	// Reset.
	//
	for (wmesh_int_t s = 0;s<select_n;++s)
	  {
	    blank[select[s]] = 0;
	  }

	//
	// Sort.
	//

	//
	// Copy back.
	//

	qsort(select,
	      select_n,
	      sizeof(wmesh_int_t),
	      sort_predicate);
	
	
	
	for (wmesh_int_t s = 0;s<select_n;++s)
	  {
	    csr_ind_[0][csr_ptr_[0][idof] + s] = select[s];
	  }
	
      }
    free(select);
    free(blank);
    free(n2c_ptr);
    free(n2c_v);
#if 0
    wmesh_status_t bms_n2c_cindex	(wmesh_int_t 		c_,
					 wmesh_int_p 		cindex_);
  
    wmesh_status_t bms_n2c_ctype	(wmesh_int_t 		c_,
					 wmesh_int_p 		ctype_);
#endif
    return WMESH_STATUS_SUCCESS;
  }
  
};
