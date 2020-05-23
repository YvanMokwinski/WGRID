#include "app.hpp"
#include <math.h>
#include <stdlib.h>
#include "bms.h"
#include "wmesh-blas.h"
#include "wmesh-math.hpp"
//#include "wmesh_nodes.h"
//#include "wmesh_nodes.h"

#include <iostream>

template<typename T>
wmesh_status_t bms_template_element_nodes(wmesh_int_t 		element_,
					  wmesh_int_t		c_storage_,
					  wmesh_int_t		c_m_,
					  wmesh_int_t		c_n_,
					  T*__restrict__ 	c_v_,
					  wmesh_int_t		c_ld_)
{

  T*__restrict__ r = c_v_ + (c_storage_ == WMESH_STORAGE_BLOCK) ? c_ld_ * 0 : 0;
  T*__restrict__ s = c_v_ + (c_storage_ == WMESH_STORAGE_BLOCK) ? c_ld_ * 1 : 1;
  T*__restrict__ t = c_v_ + (c_storage_ == WMESH_STORAGE_BLOCK) ? c_ld_ * 2 : 2;
  const wmesh_int_t inc = (c_storage_ == WMESH_STORAGE_BLOCK) ? 1 : c_ld_;
  
  switch(element_)
    {
    case WMESH_ELEMENT_NODE:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	break;
      }
  
#ifdef SET1
#error SET1 already defined
#endif
#define SET1(_i,_a)					\
  r[inc*(_i)] = (_a)

    case WMESH_ELEMENT_EDGE:
      {
	SET1(0, static_cast<T>(-1) );
	SET1(1, static_cast<T>(1)) ;
	return WMESH_STATUS_SUCCESS;
      }

#undef SET1

#ifdef SET2
#error SET2 already defined
#endif
#define SET2(_i,_a,_b)				\
      r[inc*(_i)] = (_a);			\
      s[inc*(_i)] = (_b)
      
    case WMESH_ELEMENT_TRIANGLE:
      {
	SET2(0, static_cast<T>(0), static_cast<T>(0) );
	SET2(1, static_cast<T>(1), static_cast<T>(0) );
	SET2(2, static_cast<T>(0), static_cast<T>(1) );
	return WMESH_STATUS_SUCCESS;
      }

    case WMESH_ELEMENT_QUADRILATERAL:
      {
	SET2(0, static_cast<T>(-1), static_cast<T>(-1) );
	SET2(1, static_cast<T>(1), static_cast<T>(-1) );
	SET2(2, static_cast<T>(1), static_cast<T>(1) );
	SET2(3, static_cast<T>(-1), static_cast<T>(1) );
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_TETRAHEDRON:
      {
	SET3(0, static_cast<T>(0), static_cast<T>(0), static_cast<T>(0) );
	SET3(1, static_cast<T>(1), static_cast<T>(0), static_cast<T>(0) );
	SET3(2, static_cast<T>(0), static_cast<T>(1), static_cast<T>(0) );
	SET3(3, static_cast<T>(0), static_cast<T>(0), static_cast<T>(1) );
	return WMESH_STATUS_SUCCESS;
      }

#undef SET2
      
#ifdef SET3
#error SET3 already defined
#endif
#define SET3(_i,_a,_b,_c)				\
      r[inc*(_i)] = (_a);				\
      s[inc*(_i)] = (_b);				\
      t[inc*(_i)] = (_c)
      
    case WMESH_ELEMENT_PYRAMID:
      {
	SET3(0, static_cast<T>(-1), static_cast<T>(-1), static_cast<T>(0) );
	SET3(1, static_cast<T>(1), static_cast<T>(-1), static_cast<T>(0) );
	SET3(2, static_cast<T>(1), static_cast<T>(1), static_cast<T>(0) );
	SET3(3, static_cast<T>(-1), static_cast<T>(1), static_cast<T>(0) );
	SET3(4, static_cast<T>(0), static_cast<T>(0), static_cast<T>(1) );
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_WEDGE:
      {
	SET3(0, static_cast<T>(0), static_cast<T>(0), static_cast<T>(0) );
	SET3(1, static_cast<T>(1), static_cast<T>(0), static_cast<T>(0) );
	SET3(2, static_cast<T>(0), static_cast<T>(1), static_cast<T>(0) );

	SET3(3, static_cast<T>(0), static_cast<T>(0), static_cast<T>(1) );
	SET3(4, static_cast<T>(1), static_cast<T>(0), static_cast<T>(1) );
	SET3(5, static_cast<T>(0), static_cast<T>(1), static_cast<T>(1) );
	return WMESH_STATUS_SUCCESS;
      }

    case WMESH_ELEMENT_HEXAHEDRON:
      {
	SET3(0, static_cast<T>(-1), static_cast<T>(-1), static_cast<T>(-1) );
	SET3(1, static_cast<T>(1), static_cast<T>(-1), static_cast<T>(-1) );
	SET3(2, static_cast<T>(1), static_cast<T>(1), static_cast<T>(-1) );
	SET3(3, static_cast<T>(-1), static_cast<T>(1), static_cast<T>(-1) );
	SET3(4, static_cast<T>(-1), static_cast<T>(-1), static_cast<T>(1) );
	SET3(5, static_cast<T>(1), static_cast<T>(-1), static_cast<T>(1) );
	SET3(6, static_cast<T>(1), static_cast<T>(1), static_cast<T>(1) );
	SET3(7, static_cast<T>(-1), static_cast<T>(1), static_cast<T>(1) );
	return WMESH_STATUS_SUCCESS;
      }

#undef SET3

    }

  return WMESH_STATUS_INVALID_ENUM;
}



wmesh_status_t bms_ordering_element_nodes(wmesh_int_t 	element_,
					  wmesh_int_t	c_storage_,
					  wmesh_int_t	c_m_,
					  wmesh_int_t	c_n_,
					  wmesh_int_p  	c_v_,
					  wmesh_int_t	c_ld_)
{  
  wmesh_int_p  r = c_v_ + ( (c_storage_ == WMESH_STORAGE_BLOCK) ? c_ld_ * 0 : 0);
  wmesh_int_p  s = c_v_ + ( (c_storage_ == WMESH_STORAGE_BLOCK) ? c_ld_ * 1 : 1);
  wmesh_int_p  t = c_v_ + ( (c_storage_ == WMESH_STORAGE_BLOCK) ? c_ld_ * 2 : 2);
  const wmesh_int_t inc = (c_storage_ == WMESH_STORAGE_BLOCK) ? 1 : c_ld_;
  static constexpr wmesh_int_t s_0 = 0;
  static constexpr wmesh_int_t s_1 = 1;
  switch(element_)
    {
    case WMESH_ELEMENT_NODE:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	break;
      }
  
#ifdef SET1
#error SET1 already defined
#endif
#define SET1(_i,_a)				\
      r[inc*(_i)] = (_a)
      
    case WMESH_ELEMENT_EDGE:
      {
	SET1(0, s_0 );
	SET1(1, s_1) ;
	return WMESH_STATUS_SUCCESS;
      }
      
#undef SET1
      
#ifdef SET2
#error SET2 already defined
#endif
#define SET2(_i,_a,_b)				\
      r[inc*(_i)] = (_a);			\
      s[inc*(_i)] = (_b)
      
    case WMESH_ELEMENT_TRIANGLE:
      {
	SET2(0, s_0, s_0 );
	SET2(1, s_1, s_0 );
	SET2(2, s_0, s_1 );
	return WMESH_STATUS_SUCCESS;
      }

    case WMESH_ELEMENT_QUADRILATERAL:
      {
	SET2(0, s_0, s_0 );
	SET2(1, s_1, s_0 );
	SET2(2, s_1, s_1 );
	SET2(3, s_0, s_1 );
	return WMESH_STATUS_SUCCESS;
      }
      

#undef SET2
      
#ifdef SET3
#error SET3 already defined
#endif
#define SET3(_i,_a,_b,_c)			\
      r[inc*(_i)] = (_a);			\
      s[inc*(_i)] = (_b);			\
      t[inc*(_i)] = (_c)
      
    case WMESH_ELEMENT_TETRAHEDRON:
      {
	SET3(0, s_0, s_0, s_0 );
	SET3(1, s_1, s_0, s_0 );
	SET3(2, s_0, s_1, s_0 );
	SET3(3, s_0, s_0, s_1 );
	return WMESH_STATUS_SUCCESS;
      }

    case WMESH_ELEMENT_PYRAMID:
      {
	SET3(0, s_0, s_0, s_0 );
	SET3(1, s_1, s_0, s_0 );
	SET3(2, s_1, s_1, s_0 );
	SET3(3, s_0, s_1, s_0 );
	SET3(4, s_0, s_0, s_1 );
	return WMESH_STATUS_SUCCESS;
      }

    case WMESH_ELEMENT_WEDGE:
      {
	SET3(0, s_0, s_0, s_0 );
	SET3(1, s_1, s_0, s_0 );
	SET3(2, s_0, s_1, s_0 );

	SET3(3, s_0, s_0, s_1 );
	SET3(4, s_1, s_0, s_1 );
	SET3(5, s_0, s_1, s_1 );
	return WMESH_STATUS_SUCCESS;
      }

    case WMESH_ELEMENT_HEXAHEDRON:
      {
	SET3(0, s_0, s_0, s_0 );
	SET3(1, s_1, s_0, s_0 );
	SET3(2, s_1, s_1, s_0 );
	SET3(3, s_0, s_1, s_0 );
	SET3(4, s_0, s_0, s_1 );
	SET3(5, s_1, s_0, s_1 );
	SET3(6, s_1, s_1, s_1 );
	SET3(7, s_0, s_1, s_1 );
	return WMESH_STATUS_SUCCESS;
      }
      
#undef SET3
    }

  return WMESH_STATUS_INVALID_ENUM;
}




			   
void usage(const char * appname_)
{
  fprintf(stderr,"//\n");
  fprintf(stderr,"// %s -s <element> -d <integer> -o <filename>\n",appname_);
  fprintf(stderr,"//\n");
  fprintf(stderr,"// Convert a mesh file.\n");
  fprintf(stderr,"// Example: %s example.mesh -o example.vtk\n",appname_);
  fprintf(stderr,"//\n");
  fprintf(stderr,"// Note: this is mainly used as a parrot for debugging.\n");
  fprintf(stderr,"//\n");
}

#if 0
wmesh_status_t bms_vandermonde_edge(wmesh_int_t 	trans_,
				    wmesh_int_t 	d_,
				    wmesh_int_t 	rs_m_,
				    wmesh_int_t 	rs_n_,
				    const double * 	rs_,
				    wmesh_int_t 	rs_ld_,
				    double* 		v_,
				    wmesh_int_t 	v_ld_,
				    wmesh_int_t 	work_n_,
				    double* 		work_)
{
  
  WMESH_CHECK(work_n_ >= 4 * rs_n_);
  
  double * x 	= work_ + 0*rs_n_;
  double * yi 	= work_ + 1*rs_n_;
  work_ 	= work_ + 2*rs_n_;
  work_n_ 	-=  2*rs_n_;
  wmesh_int_t idx = 0;
  for (wmesh_int_t i=0;i<=d_;++i)
    {      
      for (wmesh_int_t k = 0; k < rs_n_;++k)
	{
	  double r 	= rs_[rs_ld_*k+0];
	  double s 	= rs_[rs_ld_*k+1];
	  x[k] 		= (s < 1.0) ? (r * 2.0) / (1.0 - s) - 1.0 : 0.0;
	}
      
      bms_djacobip(0,
		   0,
		   i,		
		   rs_n_,
		   x,
		   1,
		   yi,
		   1,
		   work_n_,
		   work_);
      
      if (0 == trans_)
	{
	  for (wmesh_int_t k = 0; k < rs_n_;++k)
	    {
	      v_[v_ld_*idx+k] = yi[k];
	    }
	}
      else
	{
	  for (wmesh_int_t k = 0; k < rs_n_;++k)
	    {
	      v_[v_ld_*k+idx] = yi[k];
	    }
	}
      ++idx;    
    }
  return WMESH_STATUS_SUCCESS;
}
#endif



wmesh_status_t bms_vandermonde_triangle(wmesh_int_t 	d_,
					wmesh_int_t 	lcoo_m_,
					wmesh_int_t 	lcoo_n_,
					const double * 	lcoo_,
					wmesh_int_t 	lcoo_ld_,
					wmesh_int_t 	v_storage_,
					double* 	v_,
					wmesh_int_t 	v_ld_,
					wmesh_int_t 	work_n_,
					double* 	work_)
{
  WMESH_CHECK(work_n_ >= 5 * lcoo_n_);
  static constexpr double s1 = 1.0;
  static constexpr double s2 = 2.0;
  
  double * x 	= work_ + 0*lcoo_n_;
  double * yi 	= work_ + 1*lcoo_n_;
  double * yj 	= work_ + 2*lcoo_n_;

  work_ 	= work_ + 3*lcoo_n_;
  work_n_	-=  3*lcoo_n_;

  wmesh_int_t idx = 0;
  for (wmesh_int_t i=0;i<=d_;++i)
    {
      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
	{
	  double r 	= lcoo_[lcoo_ld_*k+0];
	  double s 	= lcoo_[lcoo_ld_*k+1];
	  x[k] 		= (s < s1) ? (r * s2) / (s1 - s) - s1 : 0.0;
	}

      bms_djacobip(0,
		   0,
		   i,		
		   lcoo_n_,
		   x,
		   1,
		   yi,
		   1,
		   work_n_,
		   work_);
      
      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
	{
	  double s 	= lcoo_[lcoo_ld_*k+1];
	  yi[k] 	*= wmesh_math<double>::xpow(s1 - s, i);
	}
      
      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
	{
	  double s 	= lcoo_[lcoo_ld_*k+1];
	  x[k] 		= s2*s - s1;
	}
      
      for (wmesh_int_t j=0;j<=d_-i;++j)
	{
	  bms_djacobip(2*i+1,
		       0,
		       j,		
		       lcoo_n_,
		       x,
		       1,
		       yj,
		       1,
		       work_n_,
		       work_);	  
	  //	  std::cout << "copy " << std::endl;
	  if (WMESH_STORAGE_BLOCK == v_storage_)
	    {
	      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
		{
		  v_[v_ld_*idx+k] = yi[k] * yj[k];
		}
	    }
	  else if (WMESH_STORAGE_INTERLEAVE == v_storage_)
	    {
	      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
		{
		  v_[v_ld_*k+idx] = yi[k] * yj[k];
		}
	    }
	  else
	    {
	      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
	    }
	  ++idx;
	}

      
    }
  return WMESH_STATUS_SUCCESS;
}

wmesh_status_t bms_vandermonde_tetrahedra(wmesh_int_t 		trans_,
					  wmesh_int_t 		d_,
					  wmesh_int_t 		rst_m_,
					  wmesh_int_t 		rst_n_,
					  const double * 	rst_,
					  wmesh_int_t 		rst_ld_,
					  double* 		v_,
					  wmesh_int_t 		v_ld_,
					  wmesh_int_t 		work_n_,
					  double* 		work_)
{  
  WMESH_CHECK(work_n_ >=  10 * rst_n_);
  
  double * xi 	= work_ + 0*rst_n_;
  double * xj 	= work_ + 1*rst_n_;
  double * xk 	= work_ + 2*rst_n_;

  double * si 	= work_ + 3*rst_n_;
  double * sj 	= work_ + 4*rst_n_;

  double * yi 	= work_ + 6*rst_n_;
  double * yj 	= work_ + 7*rst_n_;
  double * yk 	= work_ + 8*rst_n_;
  

  work_ 	= work_ + 8*rst_n_;
  work_n_      -= 10*rst_n_;
  
  for (wmesh_int_t l = 0; l < rst_n_;++l)
    {
      double r 	= rst_[rst_ld_*l+0];
      double s 	= rst_[rst_ld_*l+1];
      double t 	= rst_[rst_ld_*l+2];
      xi[l] 	= (s + t < 1.0) ? (r * 2.0) / (1.0-(s+t)) - 1.0 : 0.0;
    }
  
  for (wmesh_int_t l = 0; l < rst_n_;++l)
    {
      double s 	= rst_[rst_ld_*l+1];
      double t 	= rst_[rst_ld_*l+2];
      xj[l] 	= (t < 1.0) ? (s * 2.0) / (1.0 - t) - 1.0 : 0.0;
    }
  
  for (wmesh_int_t l = 0; l < rst_n_;++l)
    {
      double t 	= rst_[rst_ld_*l+2];
      xk[l] 	= t * 2.0 - 1.0;
    }
  
  wmesh_int_t idx = 0;
  for (wmesh_int_t i=0;i<=d_;++i)
    {      
      for (wmesh_int_t l = 0; l < rst_n_;++l)
	{
	  double s 	= rst_[rst_ld_*l+1];
	  double t 	= rst_[rst_ld_*l+2];
	  si[l] 	= wmesh_math<double>::xpow(1.0 - (s + t),i);
	}	      
      
      bms_djacobip(0,
		   0,
		   i,		
		   rst_n_,
		   xi,
		   1,
		   yi,
		   1,
		   work_n_,
		   work_);
      
      for (wmesh_int_t l = 0; l < rst_n_;++l)
	{
	  yi[l] *= si[l];
	}
      
      wmesh_int_t dj = 2*i+1;
      for (wmesh_int_t j=0;j<=d_-i;++j)
	{
	  for (wmesh_int_t l = 0; l < rst_n_;++l)
	    {	     
	      double t 	= rst_[rst_ld_*l+2];
	      sj[l] 	*= wmesh_math<double>::xpow(1.0 - t, j);
	    }
	  
	  bms_djacobip(dj,
		       0,
		       j,		
		       rst_n_,
		       xj,
		       1,
		       yj,
		       1,
		       work_n_,
		       work_);
	  
	  for (wmesh_int_t l = 0; l < rst_n_;++l)
	    {	     
	      yj[l] *= sj[l] * yi[l];
	    }

	  wmesh_int_t dk = 2*(i+j)+2;
	  for (wmesh_int_t k=0;k<=d_-i-j;++k)
	    {
	      bms_djacobip(dk,
			   0,
			   k,		
			   rst_n_,
			   xk,
			   1,
			   yk,
			   1,
			   work_n_,
			   work_);
	      
	      for (wmesh_int_t l = 0; l < rst_n_;++l)
		{	     
		  yk[l] *= yj[l];
		}
	      
	      if (0 == trans_)
		{
		  for (wmesh_int_t l = 0; l < rst_n_;++l)
		    {
		      v_[v_ld_*idx+l] = yk[l];
		    }
		}
	      else
		{
		  for (wmesh_int_t l = 0; l < rst_n_;++l)
		    {
		      v_[v_ld_*l+idx] = yk[l];
		    }
		}
	      ++idx;
	    }
	}
    }
  return WMESH_STATUS_SUCCESS;
}

#if 0
wmesh_status_t bms_vandermonde(wmesh_int_t 		element_,
			       wmesh_int_t 		d_,
			       wmesh_int_t 		ndofs_,
			       wmesh_int_t 		lcoo_storage_,
			       wmesh_int_t 		lcoo_m_,
			       wmesh_int_t 		lcoo_n_,
			       const double * 		lcoo_,
			       wmesh_int_t 		lcoo_ld_,
			       wmesh_int_t 		v_storage_,
			       double* 			v_,
			       wmesh_int_t 		v_ld_,
			       wmesh_int_t 		work_n_,
			       double* 			work_)
{
  wmesh_int_t d = 1;
  //  wmesh_int_t ndofs = ((d+1)*(d+2))/2;  
  wmesh_int_t ev_m = 2;
  wmesh_int_t ev_n = 3;
  double ev[] = {0.5,0.3,
		 0.23,0.45,
		 0.2,0.24};
  wmesh_int_t ev_ld = 2;
  wmesh_status_t status;
  double *vandermonde_dofs = (double*)malloc(sizeof(double)*ndofs * ndofs);
  double *vandermonde_eval = (double*)malloc(sizeof(double)*ndofs * ev_n);

  wmesh_int_t work_n = 0;
  switch(element_)
    {      
    case WMESH_ELEMENT_NODE:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
      }
      
    case WMESH_ELEMENT_EDGE:
      {
	work_n = 2*ndofs;
	break;
      }
      
    case WMESH_ELEMENT_TRIANGLE:
      {
	work_n 	= 5 * ndofs;
	break;
      }
      
    case WMESH_ELEMENT_QUADRILATERAL:
      {
	work_n = 10*ndofs;
	break;
      }
      
    case WMESH_ELEMENT_TETRAHEDRON:
      {
	work_n = 10*ndofs;
	break;
      }
    case WMESH_ELEMENT_PYRAMID:
      {
		work_n = 10*ndofs;
	break;
      }
    case WMESH_ELEMENT_WEDGE:
      {
		work_n = 10*ndofs;
	break;
      }
    case WMESH_ELEMENT_HEXAHEDRON:
      {
		work_n = 10*ndofs;
	break;
      }
      
    default:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
	break;
      }

    }
  
  double *    work   	= (double*)malloc(sizeof(double)*work_n);
  if (!work)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
    }

  
  wmesh_int_t rs_m   	= 2;
  wmesh_int_t rs_ld  	= 2;    
  wmesh_int_t work_n 	= 5 * ndofs;  
  switch(element_)
    {      
    case WMESH_ELEMENT_NODE:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
	break;
      }
      
    case WMESH_ELEMENT_EDGE:
      {
	break;
      }
      
    case WMESH_ELEMENT_TRIANGLE:
      {
	status = bms_vandermonde_triangle(degree_,
					  rs_m,
					  ndofs_,
					  rs,
					  rs_ld,
					  1,
					  vandermonde_dofs,
					  ndofs_,
					  work_n,
					  work);
	WMESH_STATUS_CHECK(status);
	break;
      }

    case WMESH_ELEMENT_QUADRILATERAL:
      {
	status = bms_vandermonde_quadrilateral(degree_,
					       rs_m,
					       ndofs_,
					       rs,
					       rs_ld,
					       1,
					       vandermonde_dofs,
					       ndofs_,
					       work_n,
					       work);
	WMESH_STATUS_CHECK(status);
	break;
      }
      
    case WMESH_ELEMENT_TETRAHEDRON:
      {
	status = bms_vandermonde_tetrahedron(degree_,
					     rs_m,
					     ndofs_,
					     rs,
					     rs_ld,
					     1,
					     vandermonde_dofs,
					     ndofs_,
					     work_n,
					     work);
	break;
      }
      
    case WMESH_ELEMENT_PYRAMID:
      {
	status = bms_vandermonde_pyramid(degree_,
					 rs_m,
					 ndofs_,
					 rs,
					 rs_ld,
					 1,
					 vandermonde_dofs,
					 ndofs_,
					 work_n,
					 work);
	break;
      }
    case WMESH_ELEMENT_WEDGE:
      {
	status = bms_vandermonde_wedge(degree_,
				       rs_m,
				       ndofs_,
				       rs,
				       rs_ld,
				       1,
				       vandermonde_dofs,
				       ndofs_,
				       work_n,
				       work);
	break;
      }
    case WMESH_ELEMENT_HEXAHEDRON:
      {
	status = bms_vandermonde_hexahedron(degree_,
					    rs_m,
					    ndofs_,
					    rs,
					    rs_ld,
					    1,
					    vandermonde_dofs,
					    ndofs_,
					    work_n,
					    work);
	break;
      }
      
    default:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
	break;
      }
    }

}

#endif
#if 0
  wmesh_int_t info_lapack;
  wmesh_int_p perm = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*ndofs);
  double* id = (double*)calloc(ndofs*ndofs,sizeof(double));
  for (wmesh_int_t i=0;i<ndofs;++i) id[i*ndofs+i]=1.0;
  LAPACK_dgesv((wmesh_int_p)&ndofs,
	       (wmesh_int_p)&ndofs,
	       vandermonde_dofs,
	       (wmesh_int_p)&ndofs,
	       perm,
	       id,
	       (wmesh_int_p)&ndofs,
	       (wmesh_int_p)&info_lapack);
  WMESH_CHECK( 0 == info_lapack );
#endif
#if 0
  wmesh_int_t info_lapack;
  wmesh_int_p perm 	= (wmesh_int_p)malloc(sizeof(wmesh_int_t)*ndofs);
  double* id 		= (double*)calloc(ndofs*ndofs,sizeof(double));
  for (wmesh_int_t i=0;i<ndofs;++i) id[i*ndofs+i] = 1.0;
  
  LAPACK_dgesv((wmesh_int_p)&ndofs,
	       (wmesh_int_p)&ev_n,
	       vandermonde_dofs,
	       (wmesh_int_p)&ndofs,
	       perm,
	       vandermonde_eval,
	       (wmesh_int_p)&ndofs,
	       (wmesh_int_p)&info_lapack);
  WMESH_CHECK( 0 == info_lapack );
  for (int i=0;i<ndofs;++i)
    {
      for (int j=0;j<ev_n;++j)
	{
	  std::cout << " " << vandermonde_eval[ndofs*j+i];	  
	}
      std::cout << std::endl;
    }
#endif

int main(int 		argc,
	 char** 	argv)
{

  
#if 1

#if 0
wmesh_nodes_t nodes;
  wmesh_nodes_def(&nodes,WMESH_ELEMENT_TETRAHEDRON,WMESH_NODES_FAMILY_LAGRANGE,3);
  
  std::cout <<  "d " << wmesh_nodes_degree(&nodes) << " 3" << std::endl;
  std::cout <<  "f " <<  wmesh_nodes_family(&nodes) << " " << WMESH_NODES_FAMILY_LAGRANGE << std::endl;
  std::cout <<  "e " << wmesh_nodes_element(&nodes) << " " << WMESH_ELEMENT_EDGE << std::endl;
  std::cout <<  "n " << wmesh_nodes_ndofs(&nodes) << std::endl;
#endif
#if 0
  {
    wmesh_int_t work_n = 1024;
    double c[512],work[1024];
    wmesh_int_t degree = 3;
    wmesh_int_t dim = 3;
    wmesh_int_t ndofs = 20;
    wmesh_int_t c_m = dim;
    wmesh_int_t c_n = ndofs;
    wmesh_int_t c_ld = dim;
    wmesh_status_t status = bms_dnodes(WMESH_ELEMENT_TETRAHEDRON,
				       WMESH_NODES_FAMILY_LAGRANGE,
				       degree,
				       WMESH_STORAGE_INTERLEAVE,
				       c_m,
				       c_n,
				       c,
				       c_ld,
				       work_n,
				       work);
    WMESH_STATUS_CHECK( status );
    for (wmesh_int_t j=0;j<c_n;++j)
      {
	for (wmesh_int_t i=0;i<c_m;++i)
	  {
	    std::cout << " " << c[c_ld*j+i];
	  }
	std::cout << std::endl;
      }
  }

  exit(1);
#endif  

  wmesh_int_t d = 1;
  wmesh_int_t ndofs = ((d+1)*(d+2))/2;
  
  wmesh_int_t ev_m = 2;
  wmesh_int_t ev_n = 3;
  double ev[] = {0.5,0.3,
		 0.23,0.45,
		 0.2,0.24};
  wmesh_int_t ev_ld = 2;
  
  double *vandermonde_dofs = (double*)malloc(sizeof(double)*ndofs * ndofs);
  double *vandermonde_eval = (double*)malloc(sizeof(double)*ndofs * ev_n);
  printf("yo1\n");  
  {
    // (VV^-1)^{-T}  VE^T
    double rs[]{0,0,1,0,0,1};
    wmesh_int_t rs_m  = 2;
    wmesh_int_t rs_ld = 2;    
    wmesh_int_t work_n 	= 20 * ndofs;
    double * 	work 	= (double*)malloc(sizeof(double)*work_n);
    bms_vandermonde_triangle(d,
			     rs_m,
			     ndofs,
			     rs,
			     rs_ld,
			     WMESH_STORAGE_INTERLEAVE,
			     vandermonde_dofs,
			     ndofs,
			     work_n,
			     work);
    free(work);
  }
  printf("yo2\n");  
  
  {
    wmesh_int_t work_n 	= 20 * ev_n;
    double * 	work 	= (double*)malloc(sizeof(double)*work_n);
    bms_vandermonde_triangle(d,
			     ev_m,
			     ev_n,
			     ev,
			     ev_ld,
			     WMESH_STORAGE_INTERLEAVE,
			     vandermonde_eval,
			     ndofs,
			     work_n,
			     work);
    free(work);
  }

#if 0
  wmesh_int_t info_lapack;
  wmesh_int_p perm = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*ndofs);
  double* id = (double*)calloc(ndofs*ndofs,sizeof(double));
  for (wmesh_int_t i=0;i<ndofs;++i) id[i*ndofs+i]=1.0;
  LAPACK_dgesv((wmesh_int_p)&ndofs,
	       (wmesh_int_p)&ndofs,
	       vandermonde_dofs,
	       (wmesh_int_p)&ndofs,
	       perm,
	       id,
	       (wmesh_int_p)&ndofs,
	       (wmesh_int_p)&info_lapack);
  WMESH_CHECK( 0 == info_lapack );
#endif

  wmesh_int_t info_lapack;
  wmesh_int_p perm 	= (wmesh_int_p)malloc(sizeof(wmesh_int_t)*ndofs);
  double* id 		= (double*)calloc(ndofs*ndofs,sizeof(double));
  for (wmesh_int_t i=0;i<ndofs;++i) id[i*ndofs+i] = 1.0;
  
  LAPACK_dgesv((wmesh_int_p)&ndofs,
	       (wmesh_int_p)&ev_n,
	       vandermonde_dofs,
	       (wmesh_int_p)&ndofs,
	       perm,
	       vandermonde_eval,
	       (wmesh_int_p)&ndofs,
	       (wmesh_int_p)&info_lapack);
  WMESH_CHECK( 0 == info_lapack );
  for (int i=0;i<ndofs;++i)
    {
      for (int j=0;j<ev_n;++j)
	{
	  std::cout << " " << vandermonde_eval[ndofs*j+i];	  
	}
      std::cout << std::endl;
    }


#if 0
  wmesh_int_t numPoints_=3,perm[32];
  wmesh_int_t nequal3=3,info_lapack;
  double B[]={1,0,0,0,1,0,0,0,1};
  for (int i=0;i<3;++i)
    {
      for (int j=0;j<3;++j)
	{
	  std::cout << " " << v[3*j+i];
	  
	}
      std::cout << std::endl;
    }
  std::cout << " dddddd " << std::endl;
  
  double V[32];
  for (wmesh_int_t i=0;i<3*3;++i)
    V[i]=v[i];
  LAPACK_dgesv((wmesh_int_p)&numPoints_,
	       (wmesh_int_p)&nequal3,
	       V,
	       (wmesh_int_p)&numPoints_,
	       perm,
	       B,
	       (wmesh_int_p)&numPoints_,
	       (wmesh_int_p)&info_lapack);
  std::cout << " " << info_lapack << std::endl;
  double r1=1.0;
  double r0=0.0;
  double C[32];
  BLAS_dgemm("N","N",&nequal3,&nequal3,&nequal3,&r1,B,&nequal3,v,&nequal3,&r0,C,&nequal3);
  for (int i=0;i<3;++i)
    {
      for (int j=0;j<3;++j)
	{
	  std::cout << " " << C[3*j+i];
	  
	}
      std::cout << std::endl;
    }
#endif
  exit(1);



#endif
  
  wmesh_status_t status;
  wmesh_int_t degree;
  wmesh_int_t element;
  wmesh_int_t nodes_family;  
  //
  // Command line.
  //
  WCOMMON::cmdline cmd(argc,
		       argv);
  
  if (cmd.get_nargs() == 1)
    {
      usage(argv[0]);
      return 1;
    }  
  else if ( cmd.option("-h") || cmd.option("--help") )
    {
      usage(argv[0]);
      return 0;
    }

  //
  // Get verbose.
  //
  bool verbose = cmd.option("-v");
  if (verbose)
    {
      cmd.disp_header(stdout);
    }     

  if (false == cmd.option("-d", &degree))
    {
      fprintf(stderr,"missing option, '-d <integer>'.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  else if (degree < 1)
    {
      fprintf(stderr,"wrong value from option, '-d <integer>', must be >=1.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;      
    }
  
  //
  // Get output filename.
  //
  WCOMMON::cmdline::str_t element_name;
  if (false == cmd.option("-e", element_name))
    {
      fprintf(stderr,"// wmesh-convert::error: missing element name, '-e' option.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  


  status = app_str2element(element_name,
			   &element);
  if (status != WMESH_STATUS_SUCCESS)
    {
      fprintf(stderr,"wrong element name '%s'\n",element_name);
      fprintf(stderr," available elements are:\n");
      for  (wmesh_int_t i=WMESH_ELEMENT_EDGE;i<WMESH_ELEMENT_ALL;++i)
	{
	  
	  status = app_element2str(i,
				     element_name);
	  WMESH_STATUS_CHECK(status);
	  fprintf(stderr," - '%s'\n",element_name);
	}      
      return WMESH_STATUS_INVALID_ARGUMENT;
    }

  WCOMMON::cmdline::str_t nodes_family_name;
  if (false == cmd.option("-n", nodes_family_name))
    {
      fprintf(stderr,"// bms_nodes.tests::error: missing nodes family name, '-n' option.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  
  status = app_str2nodesfamily(nodes_family_name,
			       &nodes_family);
  if (status != WMESH_STATUS_SUCCESS)
    {
      fprintf(stderr,"wrong nodes family name '%s'\n",nodes_family_name);
      fprintf(stderr," available elements are:\n");
      for  (wmesh_int_t i=0;i<WMESH_NODES_FAMILY_ALL;++i)
	{	  
	  status = app_nodesfamily2str(i,
				       nodes_family_name);
	  WMESH_STATUS_CHECK(status);
	  fprintf(stderr," - '%s'\n",nodes_family_name);
	}      
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  
  WCOMMON::cmdline::str_t ofilename;
  if (false == cmd.option("-o", ofilename))
    {
      fprintf(stderr,"// wmesh-convert::error: missing output file, '-o' option.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }

  //
  // Read the mesh.
  //
  wmesh_t* mesh;
  status = wmesh_def_rmacro(&mesh,
			    element,
			    nodes_family,
			    degree);
  WMESH_STATUS_CHECK(status);
  
  //
  // Write the mesh.
  //
  status = wmesh_write(mesh,
		       ofilename);
  WMESH_STATUS_CHECK(status);

  return 0;
}
