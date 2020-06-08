#pragma once
#include "wmesh-blas.hpp"

// #include "wls-sparse-matrix.hpp"
// #include "DiagonalBlockMatrix.hpp"
#if 0
int comp(const void * a,
	 const void * b)
{
  const_wmesh_int_p a_ = (const_wmesh_int_p)a;
  const_wmesh_int_p b_ = (const_wmesh_int_p)b;
  if (a_[0] < b_[0] )
    {
      return -1;
    }
  else if (a_[0] > b_[0] )
    {
      return 1;
    }
  else
    {
      return 0;
    }
}
#endif

struct DG_JACOBIAN
{
private:
  
  wmesh_int_t  m_nelm;
  wmesh_int_t  m_nfaceinelm;
  wmesh_int_t  m_size_block;
  wmesh_int_t  m_size_blockXsize_block;

public:

  //  WLS::sparse::matrix_t<double>* 	m_matrix;
  //  WLS::sparse::symbolic_t* 		m_sparsityPattern;
  
  wmesh_int_t  	m_n;
  wmesh_int_t  	m_nc;
  wmesh_int_p 	m_begin;
  wmesh_int_p 	m_index;
  
  double * __restrict__ 	m_values;
  double * __restrict__ 	m_original_values;
  
  double * __restrict__   	m_idiagonal {};
  double * __restrict__   	m_diagonal  {};
  
  bool 	has_inverse_diagonal{};
  
public:

  inline wmesh_int_t  nelm() const { return this->m_nelm; };
  inline wmesh_int_t  size_block() const { return m_size_block; };
  inline wmesh_int_t  size_blockXsize_block()const{return m_size_blockXsize_block;};
  
  void spy(const char * filename)
  {
    FILE * f = fopen(filename,"w");
    for (wmesh_int_t i = 0;i<m_n;++i)
      {
	for (wmesh_int_t at = m_begin[i];at<m_begin[i+1];++at)
	  {
	    fprintf(f," " WMESH_INT_FORMAT " " WMESH_INT_FORMAT " %8.15e\n",m_n - i, m_index[at] + 1,m_values[at]);
	    //  fprintf(f," " WMESH_INT_FORMAT " " WMESH_INT_FORMAT "\n",i, m_index[at]);
	  }
	fprintf(f,"\n");
      }
    fclose(f);
  };
  
  void gemv(const double * __restrict__ 	x,
	    double * __restrict__ 		y)
  {
    for (wmesh_int_t i=0;i<m_n;++i)
      {
	double s = 0.0;
	for (wmesh_int_t at = m_begin[i];at<m_begin[i+1];++at)
	  {
	    s+=x[m_index[at]] * m_original_values[at];
	  }
	y[i] = s;
      }
  };

  void gemv(const double * __restrict__ 	x,
	    wmesh_int_t 		xoff,
	    double * __restrict__ 		y,
	    wmesh_int_t 		yoff) const
  {
    for (wmesh_int_t i=0;i<m_n;++i)
      {
	wmesh_int_t row = i;
	wmesh_int_t ielm = row / m_size_block;
	wmesh_int_t iloc = row % m_size_block;
	double s = 0.0;
	for (wmesh_int_t at = m_begin[i];at<m_begin[i+1];++at)
	  {
	    wmesh_int_t col = m_index[at];
	    wmesh_int_t jelm = col / m_size_block;
	    wmesh_int_t jloc = col % m_size_block;
	    s+=x[jelm * xoff + jloc] * m_values[at];
	  }
	y[ielm*yoff + iloc] = s;
      }
  };

  

  void inverse_diagonal()
  {
    
    if (has_inverse_diagonal)
      {
	return;
      }
    const wmesh_int_t n = m_size_block;

    has_inverse_diagonal = true;
    if ( NULL == m_idiagonal)
      {
	{
	  m_original_values = (double * __restrict__)malloc(sizeof(double)*m_nc);
	  wmesh_int_t n1=1;
	  xcopy(&m_nc,m_values,&n1,m_original_values,&n1);
	}
	
	this->m_idiagonal = (double * __restrict__)malloc(sizeof(double)*m_nelm * n*n);
	this->m_diagonal  = (double * __restrict__)malloc(sizeof(double)*m_nelm * n*n);
      }
    
    //    spy("roger.txt");
    for (wmesh_int_t ielm=0;ielm<m_nelm;++ielm)
      {
	for (wmesh_int_t i = ielm*n;i<(ielm+1)*n;++i)
	  {
	    //	    printf("ielm " WMESH_INT_FORMAT " " WMESH_INT_FORMAT " " WMESH_INT_FORMAT "\n",ielm,m_begin[i],m_begin[i+1]);
	    for (wmesh_int_t at = m_begin[i];at<m_begin[i+1];++at)
	      {
		//		printf("j " WMESH_INT_FORMAT " ielm * n " WMESH_INT_FORMAT "\n",m_index[at],ielm*n);		
		if (m_index[at] == ielm * n)
		  {
		    for (wmesh_int_t j=0;j<n;++j)
		      {
			m_idiagonal[ielm * n*n + j * n + (i-ielm*n)] = m_values[at + j];
			m_diagonal[ielm * n*n + j * n + (i-ielm*n)] = m_values[at + j];
			m_values[at + j] = 0.0;
		      }
		    break;
		  }
	      }
	  }

      }
    
    wmesh_int_t lcperm[1024];
    double tmp[1024];
    for (wmesh_int_t ielm=0;ielm<m_nelm;++ielm)
      {
#if 0
	printf("BEFOdoubleE \n");
	for (wmesh_int_t j=0;j<n;++j)
	  {
	    for (wmesh_int_t i=0;i<n;++i)
	      {
		printf(" %8.15e",m_idiagonal[ielm*n*n+j*n+i]);
	      }
	    printf("\n");
	  }
	if (ielm==3)
	exit(1);
#endif
	for (wmesh_int_t j=0;j<n;++j)
	  {
	    for (wmesh_int_t i=0;i<n;++i)
	      {
		tmp[j*n+i] = m_idiagonal[ielm*this->m_size_blockXsize_block+j*this->m_size_block+i];
	      }
	  }
	for (wmesh_int_t j=0;j<n;++j)
	  {
	    for (wmesh_int_t i=0;i<n;++i)
	      {
		m_idiagonal[ielm*n*n+j*n+i] = 0.0;
	      }
	  }

	for (wmesh_int_t j=0;j<n;++j)
	  {
	    m_idiagonal[ielm*n*n+j*n+j] = 1.0;
	  }
	
	wmesh_int_t info_lapack;
	dgesv(&this->m_size_block,
	      &this->m_size_block,
	      tmp,
	      &this->m_size_block,
	      lcperm,
	      &m_idiagonal[ielm*this->m_size_blockXsize_block],
	      &this->m_size_block,
	      &info_lapack);
#if 0
	for (wmesh_int_t j=0;j<n;++j)
	  {
	    for (wmesh_int_t i=0;i<n;++i)
	      {
		printf(" %8.15e",m_idiagonal[ielm*n*n+j*n+i]);
	      }
	    printf("\n");
	  }
	printf("info_lapack %d\n",info_lapack);
#endif
	}		
  };

  void inverse_diagonal_gemv(double * __restrict__ y,const_wmesh_int_p yoff)
  {
    double tmp[128];
    wmesh_int_t n1=1;
    double r1=1.0;
    double r0=0.0;
    for (wmesh_int_t ielm=0;ielm<m_nelm;++ielm)
      {
	for (wmesh_int_t i=0;i<m_size_block;++i)
	  {
	    tmp[i] = y[yoff[0]*ielm+i];
	  }
	
	dgemv("N",
		   &this->m_size_block,
		   &this->m_size_block,
		   &r1,
		   &this->m_idiagonal[ielm*this->m_size_blockXsize_block],
		   &this->m_size_block,
		   tmp,
		   &n1,
		   &r0,
		   &y[yoff[0]*ielm],
		   &n1);
      }
  };

  void diagonal_gemv(double * __restrict__ y,const_wmesh_int_p yoff)
  {
    double tmp[128];
    wmesh_int_t n1=1;
    double r1=1.0;
    double r0=0.0;
    for (wmesh_int_t ielm=0;ielm<m_nelm;++ielm)
      {
	for (wmesh_int_t i=0;i<m_size_block;++i)
	  {
	    tmp[i] = y[yoff[0]*ielm+i];
	  }
	dgemv("N",
		   &this->m_size_block,
		   &this->m_size_block,
		   &r1,
		   &this->m_diagonal[ielm*this->m_size_blockXsize_block],
		   &this->m_size_block,
		   tmp,
		   &n1,
		   &r0,
		   &y[yoff[0]*ielm],
		   &n1);
      }
  };

  void update(wmesh_int_t s,wmesh_int_t n,double * __restrict__ x)
  {
    
    for (wmesh_int_t i = 0;i<n;++i)
      {
	//	double y=0.0;
	for (wmesh_int_t at = m_begin[s+i];at<m_begin[s+i+1];++at)
	  {
	    if (m_index[at] < s + i)
	      {
    //    printf("s+i " WMESH_INT_FORMAT " " WMESH_INT_FORMAT "\n",s+i,m_nelm*10);
#if 1
		x[s+i] -= m_values[at] * x[m_index[at]];
#endif
	      }
	  }
      }
  };

  void inverse_lower_gemv(double * __restrict__ y,const_wmesh_int_p yoff)
  {
    double tmp[128];
    wmesh_int_t n1=1;
    double r1=1.0;
    double r0=0.0;
    for (wmesh_int_t ielm=0;ielm<m_nelm;++ielm)
      {

	update(ielm*m_size_block,
	       m_size_block,
	       y);
	
	for (wmesh_int_t i=0;i<m_size_block;++i)
	  {
	    tmp[i] = y[yoff[0]*ielm+i];
	  }
	
	// -L1 * 
	// A1
	// L1 A2
	//	
	dgemv("N",
		   &this->m_size_block,
		   &this->m_size_block,
		   &r1,
		   &this->m_idiagonal[ielm*this->m_size_blockXsize_block],
		   &this->m_size_block,
		   tmp,
		   &n1,
		   &r0,
		   &y[yoff[0]*ielm],
		   &n1);
      }
  };



  void extra_gemv(double * __restrict__ x,double * __restrict__ tmp)
  {    
    for (wmesh_int_t i = 0;i<m_n;++i)
      {
	double y=0.0;
	for (wmesh_int_t at = m_begin[i];at<m_begin[i+1];++at)
	  {
	    y += m_values[at] * x[m_index[at]];
	  }
	tmp[i]=y;
      }
    for (wmesh_int_t i = 0;i<m_n;++i)
      {
	x[i] = -tmp[i];
      }    
  };

  void extra_gemv_upper(double * __restrict__ x,double * __restrict__ tmp)
  {    
    for (wmesh_int_t i = 0;i<m_n;++i)
      {
	double y=0.0;
	for (wmesh_int_t at = m_begin[i];at<m_begin[i+1];++at)
	  {
	    if (m_index[at] > i)
	      {
		y += m_values[at] * x[m_index[at]];
	      }
	  }
	tmp[i]=y;
      }
    for (wmesh_int_t i = 0;i<m_n;++i)
      {
	x[i] = -tmp[i];
      }    
  };


  DG_JACOBIAN(wmesh_int_t 	nelm_,
	      wmesh_int_t 	nfaceinelm_,
	      wmesh_int_t 	size_block_,
	      const_wmesh_int_p 	adj_,
	      const_wmesh_int_p 	adjoff_)
  {
    this->m_size_block 			= size_block_;
    this->m_size_blockXsize_block 	= size_block_*size_block_;
    this->m_nelm			= nelm_;
    this->m_nfaceinelm			= nfaceinelm_;
    this->m_n 				= this->m_size_block * this->m_nelm;

    // 
    // NC IS NOT EQUAL TO t/he NUMBER OF COEFF, YOU WILL NEED TO DEBUGx
    // 
    this->m_nc 				= this->m_size_block * this->m_size_block * (this->m_nfaceinelm+1) * this->m_nelm;
    this->m_begin 			= (wmesh_int_p)calloc(this->m_n+1,sizeof(wmesh_int_t));
    this->m_index 			= (wmesh_int_p)malloc(this->m_nc*sizeof(wmesh_int_t));
    this->m_values 			= (double * __restrict__)malloc(this->m_nc*sizeof(double));
    
    for (wmesh_int_t ielm=0;ielm<nelm_;++ielm)
      {
	//
	// Diagonal.
	//
	for (wmesh_int_t k = 0;k <  this->m_size_block;++k)
	  {
	    m_begin[ielm * m_size_block + k + 1] += m_size_block;
	  }
#if 1
	//
	// Extra diagonal.
	//
	for (wmesh_int_t localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    wmesh_int_t nei = adj_[adjoff_[0]*ielm+localFaceIndex];
	    if (nei)
	      {
		for (wmesh_int_t k = 0;k <  this->m_size_block;++k)
		  {
		    m_begin[ielm * m_size_block + k + 1] += m_size_block;
		  }
	      }
	  }
#endif
      }
    
    for (wmesh_int_t i=2;i<=this->m_n;++i)
      {
	m_begin[i]+=m_begin[i-1];
      }
#if 0
    fprintf(stdout," begin[" WMESH_INT_FORMAT "] = " WMESH_INT_FORMAT "\n",i,m_begin[i]);
    fprintf(stdout," " WMESH_INT_FORMAT " \n",m_begin[this->m_n]);
    exit(1);
#endif
    for (wmesh_int_t ielm=0;ielm<nelm_;++ielm)
      {
#if 1
	//
	// Before diagonal.
	//
	for (wmesh_int_t localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    wmesh_int_t nei = adj_[adjoff_[0]*ielm+localFaceIndex];
	    if (nei)
	      {
		if (nei-1 < ielm)
		  {
		    for (wmesh_int_t k = 0;k <  this->m_size_block;++k)
		      {
			for (wmesh_int_t j = 0;j <  this->m_size_block;++j)
			  {					    
			    m_index[m_begin[ielm * m_size_block + k]] = (nei-1) * m_size_block + j;
			    m_begin[ielm * m_size_block + k]+=1;
			  }
		      }
		  }
	      }
	  }	
#endif	
	//
	// Diagonal.
	//
	for (wmesh_int_t k = 0;k <  this->m_size_block;++k)
	  {
	    for (wmesh_int_t j = 0;j <  this->m_size_block;++j)
	      {		
		m_index[m_begin[ielm * m_size_block + k]] = ielm * m_size_block + j;
		m_begin[ielm * m_size_block + k]+=1;
	      }
	  }
#if 1
	//
	// After diagonal.
	//
	for (wmesh_int_t localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    wmesh_int_t nei = adj_[adjoff_[0]*ielm+localFaceIndex];
	    if (nei)
	      {
		if (nei - 1 > ielm)
		  {
		    for (wmesh_int_t k = 0;k <  this->m_size_block;++k)
		      {
			for (wmesh_int_t j = 0;j <  this->m_size_block;++j)
			  {					    
			    m_index[m_begin[ielm * m_size_block + k]] = (nei-1) * m_size_block + j;
			    m_begin[ielm * m_size_block + k]+=1;
			  }
		      }
		  }
	      }
	  }
#endif	
      }

    for (wmesh_int_t i=this->m_n;i>0;--i)
      {
	m_begin[i] = m_begin[i-1];
      }
    m_begin[0] = 0;

    
    for (wmesh_int_t i=0;i<this->m_n;++i)
      {
	qsort(&m_index[m_begin[i]],m_begin[i+1]-m_begin[i],sizeof(wmesh_int_t),comp);
      }
    //    this->spy("dg.txt");

#if 0
    {
      
      static constexpr bool use_fortran_indexing 	= false;
      static constexpr bool is_owner 			= true;
      this->m_sparsityPattern = new WLS::sparse::symbolic_t(use_fortran_indexing,
							    this->m_n,
							    this->m_n,
							    this->m_nc,
							    this->m_begin,
							    this->m_index,
							    is_owner);
    }
#endif
    
    //    this->m_matrix = new WLS::sparse::matrix_t<double>(this->m_sparsityPattern);
  };

#if 0
  temp_jacvar operator*(const DG_VAdouble&x_)
  {
    return {*this,x_};
  };
#endif
  
  void addelm(wmesh_int_t 	ielm_,
	      wmesh_int_t 	jelm_,
	      double         s_,
	      double * __restrict__ 	block_x_,
	      wmesh_int_t  	block_off_)
  {
    for (wmesh_int_t i = 0;i<this->m_size_block;++i)
      {
	wmesh_int_t bound = this->m_begin[ielm_*this->m_size_block + i + 1];
	for (wmesh_int_t at = this->m_begin[ielm_*this->m_size_block + i];at<bound;++at)
	  {
	    if (this->m_index[at] == jelm_ * this->m_size_block)
	      {
		for (wmesh_int_t k=0;k<this->m_size_block;++k)
		  {
		    this->m_values[at+k] = block_x_[block_off_*k+i] + s_ * this->m_values[at+k];
		  }
		break;
	      }
	  }
      }
  };

};
