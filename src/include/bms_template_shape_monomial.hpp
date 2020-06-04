#pragma once

#include "bms_template_shape_monomial_splz.hpp"

template<wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_template_shape_monomial
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	degree_,
				     const_wmesh_int_p	diff_,				     
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
		      
				     wmesh_int_t   iw_n_,
				     wmesh_int_p   iw_,
				     wmesh_int_t   rw_n_,
				     T* rw_)
  {
#ifdef FORWARD
#error FORWARD already defined
#endif
#define FORWARD					\
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

    switch(degree_)
      {
      case 0:
	{
	  return bms_template_shape_eval(FORWARD,
					 bms_template_shape_monomial_splz<0,ELEMENT_,T>::basis);
	}
      case 1:
	{
	  return bms_template_shape_eval(FORWARD,
					 bms_template_shape_monomial_splz<1,ELEMENT_,T>::basis);
	}
	
      default:
	{
	  bms_template_shape_monomial_functor<ELEMENT_,T> f(degree_);
	  return bms_template_shape_eval(FORWARD,
					 f.basis);
	}
      }

    return WMESH_STATUS_INVALID_ARGUMENT;
#undef FORWARD
    
  }

};
