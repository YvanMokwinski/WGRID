#include "bms.hpp"
#include "wmesh-math.hpp"
#include "wmesh-blas.h"
#include <iostream>
#ifndef NDEBUG
#include <iostream>
#endif

template<typename T,typename F>  
static inline  wmesh_status_t bms_template_shape_eval(const_wmesh_int_p	diff_,
						      wmesh_int_t 	c_storage_,
						      
						      wmesh_int_t 	c_m_,
						      wmesh_int_t 	c_n_,
						      const T * 	c_,
						      wmesh_int_t 	c_ld_,

						      
						      wmesh_int_t 	b_storage_,
						      wmesh_int_t 	b_m_,
						      wmesh_int_t 	b_n_,
						      T* 		b_,
						      wmesh_int_t 	b_ld_,
						      
						      wmesh_int_t  	iw_n_,
						      wmesh_int_p  	iw_,
						      wmesh_int_t  	rw_n_,
						      T*		rw_,
						      F 		funcbasis)
  {
    
    const wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;
    
    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {

	      funcbasis(diff_,
			c_ + c_ld_*i,
			1,
			b_ + b_ld_*i,
			1);

	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {	      
	      funcbasis(diff_,
		    c_ + c_ld_*i,
		    1,
			
		    b_ + i,
			b_ld_);

	      
	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      funcbasis(diff_,
		    c_ + i,
		    c_ld_,
		    b_ + b_ld_ * i,
		    1);
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      funcbasis(diff_,
		    c_ + i,
		    c_ld_,
		    b_ + i,
		    b_ld_);
	    }
	  return WMESH_STATUS_SUCCESS;
	}








	
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;

    
  };


template<wmesh_int_t FAMILY_,typename T>
struct  bms_shape;


#include "bms_template_shape_jacobi.hpp"
#include "bms_template_shape_lagrange.hpp"
#include "bms_template_shape_legendre.hpp"
#include "bms_template_shape_orthogonal.hpp"




template<typename T>
struct 
bms_shape<WMESH_SHAPE_FAMILY_LAGRANGE,T>
{
  static inline  wmesh_status_t eval(wmesh_int_t 		element_,
				     wmesh_int_t 		degree_,
				     const_wmesh_int_p		diff_,				     
		      
				     wmesh_int_t 		c_storage_,					  
				     wmesh_int_t 		c_m_,
				     wmesh_int_t 		c_n_,
				     const T * 			c_,
				     wmesh_int_t 		c_ld_,
		      
				     wmesh_int_t 		b_storage_,
				     wmesh_int_t 		b_m_,
				     wmesh_int_t 		b_n_,
				     T* 			b_,
				     wmesh_int_t 		b_ld_,
		      
		      
				     wmesh_int_t		iw_n_,
				     wmesh_int_p		iw_,
				     wmesh_int_t		rw_n_,
				     T*		rw_)
  {

#define FORWARD					\
    degree_,					\
      diff_,					\
      c_storage_,				\
      c_m_,					\
      c_n_,					\
      c_,					\
      c_ld_,					\
      						\
      b_storage_,				\
      b_m_,					\
      b_n_,					\
      b_,					\
      b_ld_,					\
      						\
      iw_n_,					\
      iw_,					\
      rw_n_,					\
      rw_

    switch(element_)
      {
#define TREAT_CASE(_c) case _c: return bms_template_shape_lagrange<_c,T>::eval(FORWARD)

	TREAT_CASE(WMESH_ELEMENT_EDGE);
	TREAT_CASE(WMESH_ELEMENT_TRIANGLE);
	TREAT_CASE(WMESH_ELEMENT_QUADRILATERAL);
	TREAT_CASE(WMESH_ELEMENT_TETRAHEDRON);
	TREAT_CASE(WMESH_ELEMENT_PYRAMID);
	TREAT_CASE(WMESH_ELEMENT_WEDGE);
	TREAT_CASE(WMESH_ELEMENT_HEXAHEDRON);
#undef TREAT_CASE
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef FORWARD
  }
  
};


template<typename T>
struct 
bms_shape<WMESH_SHAPE_FAMILY_LEGENDRE,T>
{
  static inline      wmesh_status_t eval(wmesh_int_t 		element_,
					 wmesh_int_t 		degree_,
					 const_wmesh_int_p	diff_,				     
			
					 wmesh_int_t 		c_storage_,					  
					 wmesh_int_t 		c_m_,
					 wmesh_int_t 		c_n_,
					 const T * 		c_,
					 wmesh_int_t 		c_ld_,
			
					 wmesh_int_t 		b_storage_,
					 wmesh_int_t 		b_m_,
					 wmesh_int_t 		b_n_,
					 T* 			b_,
					 wmesh_int_t 		b_ld_,
			
					 wmesh_int_t		iw_n_,
					 wmesh_int_p		iw_,
					 wmesh_int_t		rw_n_,
					 T*		rw_)
  {
      
#define FORWARD					\
    degree_,					\
      diff_,					\
      c_storage_,				\
      c_m_,					\
      c_n_,					\
      c_,					\
      c_ld_,					\
      						\
      b_storage_,				\
      b_m_,					\
      b_n_,					\
      b_,					\
      b_ld_,					\
      						\
      iw_n_,					\
      iw_,					\
      rw_n_,					\
      rw_
      

    switch(element_)
      {
#define TREAT_CASE(_c) case _c: return bms_template_shape_legendre<_c,T>::eval(FORWARD)

	TREAT_CASE(WMESH_ELEMENT_EDGE);
	TREAT_CASE(WMESH_ELEMENT_TRIANGLE);
	TREAT_CASE(WMESH_ELEMENT_QUADRILATERAL);
	TREAT_CASE(WMESH_ELEMENT_TETRAHEDRON);
	TREAT_CASE(WMESH_ELEMENT_PYRAMID);
	TREAT_CASE(WMESH_ELEMENT_WEDGE);
	TREAT_CASE(WMESH_ELEMENT_HEXAHEDRON);
#undef TREAT_CASE
      }
    return WMESH_STATUS_INVALID_ENUM;      
#undef FORWARD
  }
};


template<typename T>
struct 
bms_shape<WMESH_SHAPE_FAMILY_ORTHOGONAL,T>
{
  static inline  wmesh_status_t eval(wmesh_int_t 		element_,
				     wmesh_int_t 		degree_,
				     const_wmesh_int_p		diff_,				     
			
				     wmesh_int_t 		c_storage_,					  
				     wmesh_int_t 		c_m_,
				     wmesh_int_t 		c_n_,
				     const T * 			c_,
				     wmesh_int_t 		c_ld_,
			
				     wmesh_int_t 		b_storage_,
				     wmesh_int_t 		b_m_,
				     wmesh_int_t 		b_n_,
				     T* 			b_,
				     wmesh_int_t 		b_ld_,
				     wmesh_int_t		iw_n_,
				     wmesh_int_p		iw_,
				     wmesh_int_t		rw_n_,
				     T*		rw_)
  {
      
#define FORWARD					\
    degree_,					\
      diff_,					\
      c_storage_,				\
      c_m_,					\
      c_n_,					\
      c_,					\
      c_ld_,					\
      						\
      b_storage_,				\
      b_m_,					\
      b_n_,					\
      b_,					\
      b_ld_,					\
      						\
      iw_n_,					\
      iw_,					\
      rw_n_,					\
      rw_
      

    switch(element_)
      {
#define TREAT_CASE(_c) case _c: return bms_template_shape_orthogonal<_c,T>::eval(FORWARD)
	  
	TREAT_CASE(WMESH_ELEMENT_EDGE);
	TREAT_CASE(WMESH_ELEMENT_TRIANGLE);
	TREAT_CASE(WMESH_ELEMENT_QUADRILATERAL);
	TREAT_CASE(WMESH_ELEMENT_TETRAHEDRON);
	TREAT_CASE(WMESH_ELEMENT_PYRAMID);
	TREAT_CASE(WMESH_ELEMENT_WEDGE);
	TREAT_CASE(WMESH_ELEMENT_HEXAHEDRON);
#undef TREAT_CASE
      }

    return WMESH_STATUS_INVALID_ENUM;      
#undef FORWARD
  }

    
};





template<typename T>
wmesh_status_t bms_template_shape(wmesh_int_t 		element_,
				  wmesh_int_t 		family_,
				  wmesh_int_t 		degree_,

				  const_wmesh_int_p	diff_,
				  
				  wmesh_int_t 		c_storage_,					  
				  wmesh_int_t 		c_m_,
				  wmesh_int_t 		c_n_,
				  const T * 		c_,
				  wmesh_int_t 		c_ld_,
				  
				  wmesh_int_t 		b_storage_,
				  wmesh_int_t 		b_m_,
				  wmesh_int_t 		b_n_,
				  T* 			b_,
				  wmesh_int_t 		b_ld_,
				      
				  wmesh_int_t		iw_n_,
				  wmesh_int_p		iw_,
				  wmesh_int_t		rw_n_,
				  T * rw_)
{
      
#define FORWARD					\
  element_,					\
    degree_,					\
    diff_,					\
						\
    c_storage_,					\
    c_m_,					\
    c_n_,					\
    c_,						\
    c_ld_,					\
						\
    b_storage_,					\
    b_m_,					\
    b_n_,					\
    b_,						\
    b_ld_,					\
      						\
    iw_n_,					\
    iw_,					\
    rw_n_,					\
    rw_

  
  return bms_shape<WMESH_SHAPE_FAMILY_LAGRANGE,T>::eval(FORWARD);
#if 0
  switch(element_)
    {
      
#define TREAT_CASE(_f) case _f: return bms_shape<_f,T>::eval(FORWARD)
      
      TREAT_CASE(WMESH_SHAPE_FAMILY_LAGRANGE);
      TREAT_CASE(WMESH_SHAPE_FAMILY_LEGENDRE);
      TREAT_CASE(WMESH_SHAPE_FAMILY_ORTHOGONAL);
      
#undef TREAT_CASE
      
    }
#endif
  return WMESH_STATUS_INVALID_ENUM;
  
#undef FORWARD

}

template
wmesh_status_t bms_template_shape<float>(wmesh_int_t 		element_,
				  wmesh_int_t 		family_,
				  wmesh_int_t 		degree_,

				  const_wmesh_int_p	diff_,
				  
				  wmesh_int_t 		c_storage_,					  
				  wmesh_int_t 		c_m_,
				  wmesh_int_t 		c_n_,
				  const float * 		c_,
				  wmesh_int_t 		c_ld_,
				  
				  wmesh_int_t 		b_storage_,
				  wmesh_int_t 		b_m_,
				  wmesh_int_t 		b_n_,
				  float* 			b_,
				  wmesh_int_t 		b_ld_,
				      
				  wmesh_int_t		iw_n_,
				  wmesh_int_p		iw_,
				  wmesh_int_t		rw_n_,
				  float * rw_);

template
wmesh_status_t bms_template_shape<double>(wmesh_int_t 		element_,
				  wmesh_int_t 		family_,
				  wmesh_int_t 		degree_,

				  const_wmesh_int_p	diff_,
				  
				  wmesh_int_t 		c_storage_,					  
				  wmesh_int_t 		c_m_,
				  wmesh_int_t 		c_n_,
				  const double * 		c_,
				  wmesh_int_t 		c_ld_,
				  
				  wmesh_int_t 		b_storage_,
				  wmesh_int_t 		b_m_,
				  wmesh_int_t 		b_n_,
				  double* 			b_,
				  wmesh_int_t 		b_ld_,
				      
				  wmesh_int_t		iw_n_,
				  wmesh_int_p		iw_,
				  wmesh_int_t		rw_n_,
				  double * rw_);


extern "C" wmesh_status_t bms_shape_buffer_size(wmesh_int_t 		element_,
						wmesh_int_t 		family_,
						wmesh_int_t 		degree_,
						wmesh_int_p 		iw_n_,
						wmesh_int_p 		rw_n_)
{
  iw_n_[0] = 0;
  rw_n_[0] = 0;
  return WMESH_STATUS_SUCCESS;
}

extern "C" wmesh_status_t bms_dshape(wmesh_int_t 		element_,
				     wmesh_int_t 		family_,
				     wmesh_int_t 		degree_,




				     
				     const_wmesh_int_p		diff_,
				     				     
				     wmesh_int_t 		c_storage_,					  
				     wmesh_int_t 		c_m_,
				     wmesh_int_t 		c_n_,
				     const double * 		c_,
				     wmesh_int_t 		c_ld_,
					    
				     wmesh_int_t 		b_storage_,
				     wmesh_int_t 		b_m_,
				     wmesh_int_t 		b_n_,
				     double* 			b_,
				     wmesh_int_t 		b_ld_,
				     
				     wmesh_int_t		iw_n_,
				     wmesh_int_p		iw_,
				     wmesh_int_t		rw_n_,
				     double*		rw_)
{    
  return bms_template_shape(element_,
			    family_,
			    degree_,




			    
			    diff_,


			    
			    c_storage_,
			    c_m_,
			    c_n_,
			    c_,
			    c_ld_,
				
			    b_storage_,
			    b_m_,
			    b_n_,
			    b_,
			    b_ld_,
				
			    iw_n_,
			    iw_,
			    rw_n_,
			    rw_);

}



  

  
