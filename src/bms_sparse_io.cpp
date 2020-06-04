
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

  inline int wmesh_sscanf(const char * str_,wmesh_int_t * i_,wmesh_int_t * j_)
  {
    return sscanf(str_,WMESH_INT_FORMAT " " WMESH_INT_FORMAT, i_, j_);
  };

  inline int wmesh_sscanf(const char * str_,wmesh_int_t * m_,wmesh_int_t * n_,wmesh_int_t * nnz_)
  {
    return sscanf(str_,WMESH_INT_FORMAT " " WMESH_INT_FORMAT " " WMESH_INT_FORMAT, m_, n_, nnz_);
  };

  template<typename real_t>
  inline int wmesh_sscanf(const char * str_,wmesh_int_t * i_,wmesh_int_t * j_,real_t*x_);

  template<typename real_t>
  inline int wmesh_sscanf(const char * str_,real_t*x_);

  template<>
  inline int wmesh_sscanf<double>(const char * str_,double*x_)
  {
    return sscanf(str_, "%le" , x_);
  }

  template<>
  inline int wmesh_sscanf<float>(const char * str_,float*x_)
  {
    return sscanf(str_, "%e" , x_);
  }

  template<>
  inline int wmesh_sscanf<double>(const char * str_,wmesh_int_t * i_,wmesh_int_t * j_,double*x_)
  {
    return sscanf(str_,WMESH_INT_FORMAT " " WMESH_INT_FORMAT " %le" , i_, j_, x_);
  }

  template<>
  inline int wmesh_sscanf<float>(const char * str_,wmesh_int_t * i_,wmesh_int_t * j_,float*x_)
  {
    return sscanf(str_,WMESH_INT_FORMAT " " WMESH_INT_FORMAT " %e" , i_, j_, x_);
  }


template<typename T>
wmesh_status_t bms_template_matrix_market_dense_read(wmesh_int_p 	m_,
						     wmesh_int_p 	n_,
						     T **   	v_,
						     wmesh_int_p 	ld_,
						     const char * 	filename_)
{
  FILE * f = fopen(filename_,"r");
  MM_typecode matcode;
  
  if (mm_read_banner(f, &matcode) != 0)
    {
      std::cerr << "cannot read Matrix Market banner" << std::endl;
      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
    }
  
  
  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */
  
  if (mm_is_matrix(matcode) && !mm_is_array(matcode) )
    {
      std::cerr << "no dense matrix found." << std::endl;
      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
    }
  
  /* find out size of sparse matrix .... */
  
  wmesh_int_t mat_nrows, mat_ncols;
  int ret_code = mm_read_mtx_array_size(f, &mat_nrows, &mat_ncols);
  if (ret_code != 0)
    {
      std::cerr << "cannot read dense matrix dimensions" << std::endl;
      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
    }
  
  T * x = (T*)malloc(sizeof(T)*mat_nrows*mat_ncols);
  if (!x)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
    }
  
  m_[0] = mat_nrows;
  n_[0] = mat_ncols;
  v_[0] = x;
  ld_[0] = mat_nrows;
  
  { size_t len;
    char * line = nullptr;
    for (wmesh_int_t irow=0;irow<mat_nrows;++irow)
      {
	ssize_t read = getline(&line, &len, f);
	if (read==-1)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	  }
	for (wmesh_int_t icol=0;icol<mat_ncols;++icol)
	  {
	    if (1 != wmesh_sscanf(line,&x[mat_nrows * icol + irow]))
	      {
		WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	      }
	  }  
      } if (line) free(line); }
  
  return WMESH_STATUS_SUCCESS;
};


extern "C"
{

  wmesh_status_t bms_matrix_market_dense_dread(wmesh_int_p 	m_,
					       wmesh_int_p 	n_,
					       double **   	v_,
					       wmesh_int_p 	ld_,
					       const char * 	filename_,
					       ...)
  {
    wmesh_str_t filename;
    
    { va_list args;
      va_start(args,filename_);
      vsprintf(filename,filename_,args);
      va_end(args); }

    return bms_template_matrix_market_dense_read(m_,
						 n_,
						 v_,
						 ld_,
						 filename);

  }

  wmesh_status_t bms_matrix_market_dense_fread(wmesh_int_p 	m_,
					       wmesh_int_p 	n_,
					       float **   	v_,
					       wmesh_int_p 	ld_,
					       const char * 	filename_,
					       ...)
  {
    WMESH_CHECK_POINTER( filename_ );
    wmesh_str_t filename;
    
    { va_list args;
      va_start(args,filename_);
      vsprintf(filename,filename_,args);
      va_end(args); }

    return bms_template_matrix_market_dense_read(m_,
						 n_,
						 v_,
						 ld_,
						 filename);
  }
  
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
