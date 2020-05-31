
#include "bms.h"
#include "wmesh-utils.hpp"
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdarg.h>
#include "mmio.h"

template<typename T>
static void fprintf_real(FILE * out_,T v_);

template<>
 void fprintf_real<double>(FILE * out_,double v_)
{
  
  fprintf(out_, "%8.15e", v_);
}

template<>
 void fprintf_real<float>(FILE * out_,float v_)
{
  
  fprintf(out_, "%6.8e", v_);
}

template<typename T>
wmesh_status_t bms_template_matrix_market_dense_write(wmesh_int_t 			m_,
						      wmesh_int_t 			n_,
						      T * __restrict__  		v_,
						      wmesh_int_t 			ld_,
						      const char * __restrict__ 	filename_,
						       ...)
{
  WMESH_CHECK_POSITIVE(m_);
  WMESH_CHECK_POSITIVE(n_);
  WMESH_CHECK_POINTER(v_);
  WMESH_CHECK( m_ <= ld_ );
  WMESH_CHECK_POINTER( filename_ );
    
  wmesh_str_t filename;
    
  { va_list args;
    va_start(args,filename_);
    vsprintf(filename,filename_,args);
    va_end(args); }

  MM_typecode matcode;
  FILE * out = fopen(filename,"w");
  mm_initialize_typecode(&matcode);
  mm_set_matrix(&matcode);
  mm_set_array(&matcode);
  mm_set_real(&matcode);
  mm_write_banner(out, matcode);
  mm_write_mtx_array_size(out, m_, n_);  
  
  for (wmesh_int_t i=0;i < m_;++i)
    {
      for (wmesh_int_t j=0;j<n_;++j)
	{
	  fprintf(out," ");
	  fprintf_real(out,v_[ld_*j+i]);
	}
      fprintf(out,"\n");
    }
  fclose(out);  
    
  return WMESH_STATUS_SUCCESS;
}

template<typename T>
wmesh_status_t bms_template_matrix_market_csr_write(wmesh_int_t 	m_,
						    wmesh_int_t 	n_,
						    wmesh_int_t 	nnz_,
						    wmesh_int_p 	csr_ptr_,
						    wmesh_int_p		csr_ind_,
						    T*__restrict__	csr_val_,
						    const char * __restrict__ filename_,
						    ...)
{
  WMESH_CHECK_POSITIVE(m_);
  WMESH_CHECK_POSITIVE(n_);
  WMESH_CHECK_POSITIVE(nnz_);
  WMESH_CHECK_POINTER(csr_ptr_);
  WMESH_CHECK_POINTER(csr_val_);
  WMESH_CHECK_POINTER( filename_ );
  
  wmesh_str_t filename;
    
  { va_list args;
    va_start(args,filename_);
    vsprintf(filename,filename_,args);
    va_end(args); }

  MM_typecode matcode;

  FILE * out = fopen(filename,"w");
  mm_initialize_typecode(&matcode);
  mm_set_matrix(&matcode);
  mm_set_coordinate(&matcode);
  mm_set_real(&matcode);
  mm_write_banner(out, matcode);
  mm_write_mtx_crd_size(out, m_, n_, nnz_);  
  for (wmesh_int_t i=0;i<m_;++i)
    {
      for (wmesh_int_t k=csr_ptr_[i];k<csr_ptr_[i+1];++k)
	{
	  const wmesh_int_t j = csr_ind_[k]-1;
	  const T v = csr_val_[k];
	  fprintf(out, " " WMESH_INT_FORMAT " " WMESH_INT_FORMAT " ", i+1, j+1);	  
	  fprintf_real(out, v);
	  fprintf(out,"\n");
	}
    }  
  fclose(out);
  return WMESH_STATUS_SUCCESS;
}



extern "C"
{
  wmesh_status_t bms_matrix_market_dense_dwrite(wmesh_int_t 	m_,
						wmesh_int_t 	n_,
						double *   	v_,
						wmesh_int_t 	ld_,
						const char * 	filename_,
						...)
  {
    WMESH_CHECK_POINTER( filename_ );
    wmesh_str_t filename;
    
    { va_list args;
      va_start(args,filename_);
      vsprintf(filename,filename_,args);
      va_end(args); }

    return bms_template_matrix_market_dense_write(m_,
						  n_,
						  v_,
						  ld_,
						  filename);
  }
  

  wmesh_status_t bms_matrix_market_dense_fwrite(wmesh_int_t 	m_,
						wmesh_int_t 	n_,
						float *   	v_,
						wmesh_int_t 	ld_,
						const char * 	filename_,
						...)
  {
    WMESH_CHECK_POINTER( filename_ );
    wmesh_str_t filename;
    
    { va_list args;
      va_start(args,filename_);
      vsprintf(filename,filename_,args);
      va_end(args); }
    
    return bms_template_matrix_market_dense_write(m_,
						  n_,
						  v_,
						  ld_,
						  filename);
  }

  wmesh_status_t bms_matrix_market_csr_dwrite(wmesh_int_t 	m_,
					      wmesh_int_t 	n_,
					      wmesh_int_t 	nnz_,
					      wmesh_int_p 	csr_ptr_,
					      wmesh_int_p	csr_ind_,
					      double *__restrict__	csr_val_,
					      const char * 	filename_,
					      ...)
  {
    WMESH_CHECK_POINTER( filename_ );
    wmesh_str_t filename;
    
    { va_list args;
      va_start(args,filename_);
      vsprintf(filename,filename_,args);
      va_end(args); }
    
    return bms_template_matrix_market_csr_write(m_,
						n_,
						nnz_,
						csr_ptr_,
						csr_ind_,
						csr_val_,
						filename);
  }

  wmesh_status_t bms_matrix_market_csr_fwrite(wmesh_int_t 	m_,
					      wmesh_int_t 	n_,
					      wmesh_int_t 	nnz_,
					      wmesh_int_p 	csr_ptr_,
					      wmesh_int_p	csr_ind_,
					      float *__restrict__	csr_val_,
					      const char * 	filename_,
					      ...)
  {
    WMESH_CHECK_POINTER( filename_ );
    wmesh_str_t filename;
    
    { va_list args;
      va_start(args,filename_);
      vsprintf(filename,filename_,args);
      va_end(args); }
    
    return bms_template_matrix_market_csr_write(m_,
						n_,
						nnz_,
						csr_ptr_,
						csr_ind_,
						csr_val_,
						filename);
  }

};
