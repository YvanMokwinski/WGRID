#include "bms.hpp"
#include "wmesh-math.hpp"
#include "wmesh-blas.h"
#include <iostream>
#ifndef NDEBUG
#include <iostream>
#endif


template<wmesh_int_t FAMILY_,typename T>
struct  bms_shape;



template<wmesh_int_t 	degree_,
	 wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_shape_orthogonal_splz
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
			     wmesh_int_t 	c_m_,
			     wmesh_int_t 	c_n_,
			     const T * 	c_,
			     wmesh_int_t 	c_ld_,
		      
			     wmesh_int_t 	b_storage_,
			     wmesh_int_t 	b_m_,
			     wmesh_int_t 	b_n_,
			     T* 		b_,
			     wmesh_int_t 	b_ld_,
		      
			     wmesh_int_t   	iw_n_,
			     wmesh_int_p   	iw_,
			     wmesh_int_t   	rw_n_,
			     wmesh_int_p   	rw_);
};

template<wmesh_int_t ELEMENT_,typename T>
struct bms_shape_orthogonal_splz<0,ELEMENT_,T>
{
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
		      wmesh_int_p   rw_)
  {
    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;
    
    static constexpr T one = static_cast<T>(1);
    
#define l0       one
    
    switch(treat_case)
      {
	
      case 1:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      b_[b_ld_ * i + 0] = l0;	     
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      b_[b_ld_ * 0 + i] = l0;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      b_[b_ld_ * i + 0] = l0;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      b_[b_ld_ * 0 + i] = l0;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0

  }
};


template<typename T>
struct bms_shape_orthogonal_splz<1,WMESH_ELEMENT_EDGE,T>
{
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						  wmesh_int_p  	rw_)
  {
    
    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;

    static constexpr T one = static_cast<T>(1);
    static constexpr T two = static_cast<T>(2);

#define l0       ( one - r ) / two
#define l1       ( one + r ) / two

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1

  };
};

template<typename T>
struct bms_shape_orthogonal_splz<1,WMESH_ELEMENT_TRIANGLE,T>
{
  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						      wmesh_int_p  	rw_)
  {
    
    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;
    

    static constexpr T one = static_cast<T>(1);

#define l0       (one - (r+s))
#define l1       r
#define l2       s


    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2

  };
};

template<typename T>
struct bms_shape_orthogonal_splz<1,WMESH_ELEMENT_QUADRILATERAL,T>
{


  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
		      wmesh_int_p  	rw_)
  {

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;

        static constexpr T four = static_cast<T>(4);
    static constexpr T one = static_cast<T>(1);

#define l0       (one - r)*(one - s) / four
#define l1       (one + r)*(one + s) / four
#define l2       (one + r)*(one + s) / four
#define l3       (one - r)*(one + s) / four

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3

  };


};

template<typename T>
struct bms_shape_orthogonal_splz<1,WMESH_ELEMENT_TETRAHEDRON,T>
{


  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
							 wmesh_int_p  	rw_)
  {

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;
    
    static constexpr T one = static_cast<T>(1);
    
#define l0      ( one - (r + s + t) )
#define l1       r
#define l2       s
#define l3       t

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i],
		t = c_[c_ld_* 2  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i],
		t = c_[c_ld_* 2 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3

  };


};

template<typename T>
struct bms_shape_orthogonal_splz<1,WMESH_ELEMENT_WEDGE,T>
{

  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						     wmesh_int_p  	rw_)
  {
    static constexpr T one = static_cast<T>(1);

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;

#define l0     ( one - (r + s) ) * (one - t)
#define l1       r  * (one - t);
#define l2       s  * (one - t)
#define l3       ( one - (r + s) ) * t
#define l4       r * t
#define l5       s * t

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	      b_[b_ld_ * i + 5] = l5;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;
	      b_[b_ld_ * 5 + i] = l5;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i],
		t = c_[c_ld_* 2  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	      b_[b_ld_ * i + 5] = l5;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i],
		t = c_[c_ld_* 2 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;
	      b_[b_ld_ * 5 + i] = l5;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3
#undef l4
#undef l5

  };

};

template<typename T>
struct bms_shape_orthogonal_splz<1,WMESH_ELEMENT_HEXAHEDRON,T>
{

  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						     wmesh_int_p  	rw_)
  {
    static constexpr T one = static_cast<T>(1);
    static constexpr T eight = static_cast<T>(8);

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;


    	      
#define l0     (one - r)*(one - s)*(one - t) / eight
#define l1       (one + r)*(one - s)*(one - t) / eight
#define l2       (one + r)*(one + s)*(one - t) / eight
#define l3       (one - r)*(one + s)*(one - t) / eight
#define l4       (one - r)*(one - s)*(one + t) / eight
#define l5       (one + r)*(one - s)*(one + t) / eight
#define l6       (one + r)*(one + s)*(one + t) / eight
#define l7       (one - r)*(one + s)*(one + t) / eight

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	      b_[b_ld_ * i + 5] = l5;
	      b_[b_ld_ * i + 6] = l6;
	      b_[b_ld_ * i + 7] = l7;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;
	      b_[b_ld_ * 5 + i] = l5;
	      b_[b_ld_ * 6 + i] = l6;
	      b_[b_ld_ * 7 + i] = l7;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i],
		t = c_[c_ld_* 2  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	      b_[b_ld_ * i + 5] = l5;
	      b_[b_ld_ * i + 6] = l6;
	      b_[b_ld_ * i + 7] = l7;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i],
		t = c_[c_ld_* 2 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;
	      b_[b_ld_ * 5 + i] = l5;
	      b_[b_ld_ * 6 + i] = l6;
	      b_[b_ld_ * 7 + i] = l7;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3
#undef l4
#undef l5
#undef l6
#undef l7
  };
};

template<typename T>
struct bms_shape_orthogonal_splz<1,WMESH_ELEMENT_PYRAMID,T>
{
  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						     wmesh_int_p  	rw_)
  {
    static constexpr T one = static_cast<T>(1);
    static constexpr T four = static_cast<T>(4);

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;

#define l0 under ? ( (one - t) * (one - r)*(one - s) / four) : 0
#define l1 under ? ( (one - t) * (one + r)*(one - s) / four) : 0
#define l2 under ? ( (one - t) * (one + r)*(one + s) / four) : 0
#define l3 under ? ( (one - t) * (one - r)*(one + s) / four) : 0
#define l4 t

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      const bool under = (t < 1.0);

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      const bool under = (t < 1.0);

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i],
		t = c_[c_ld_* 2  + i];

	      const bool under = (t < 1.0);

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i],
		t = c_[c_ld_* 2 + i];


	      const bool under = (t < 1.0);

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3
#undef l4
  };
};


template<wmesh_int_t 	degree_,
	 wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_shape_lagrange_splz
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
		      wmesh_int_p   rw_);
};

template<wmesh_int_t ELEMENT_,typename T>
struct bms_shape_lagrange_splz<0,ELEMENT_,T>
{
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
		      wmesh_int_p   rw_)
  {
    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;
    
    static constexpr T one = static_cast<T>(1);
    
#define l0       one
    
    switch(treat_case)
      {
	
      case 1:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      b_[b_ld_ * i + 0] = l0;	     
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      b_[b_ld_ * 0 + i] = l0;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      b_[b_ld_ * i + 0] = l0;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      b_[b_ld_ * 0 + i] = l0;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0

  }
};


template<typename T>
struct bms_shape_lagrange_splz<1,WMESH_ELEMENT_EDGE,T>
{
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						  wmesh_int_p  	rw_)
  {
    
    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;

    static constexpr T one = static_cast<T>(1);
    static constexpr T two = static_cast<T>(2);

#define l0       ( one - r ) / two
#define l1       ( one + r ) / two

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1

  };
};

template<typename T>
struct bms_shape_lagrange_splz<1,WMESH_ELEMENT_TRIANGLE,T>
{
  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						      wmesh_int_p  	rw_)
  {
    
    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;
    

    static constexpr T one = static_cast<T>(1);

#define l0       (one - (r+s))
#define l1       r
#define l2       s


    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2

  };
};

template<typename T>
struct bms_shape_lagrange_splz<1,WMESH_ELEMENT_QUADRILATERAL,T>
{


  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
		      wmesh_int_p  	rw_)
  {

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;

        static constexpr T four = static_cast<T>(4);
    static constexpr T one = static_cast<T>(1);

#define l0       (one - r)*(one - s) / four
#define l1       (one + r)*(one + s) / four
#define l2       (one + r)*(one + s) / four
#define l3       (one - r)*(one + s) / four

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3

  };


};

template<typename T>
struct bms_shape_lagrange_splz<1,WMESH_ELEMENT_TETRAHEDRON,T>
{


  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
							 wmesh_int_p  	rw_)
  {

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;
    
    static constexpr T one = static_cast<T>(1);
    
#define l0      ( one - (r + s + t) )
#define l1       r
#define l2       s
#define l3       t

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i],
		t = c_[c_ld_* 2  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i],
		t = c_[c_ld_* 2 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3

  };


};

template<typename T>
struct bms_shape_lagrange_splz<1,WMESH_ELEMENT_WEDGE,T>
{

  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						     wmesh_int_p  	rw_)
  {
    static constexpr T one = static_cast<T>(1);

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;

#define l0     ( one - (r + s) ) * (one - t)
#define l1       r  * (one - t);
#define l2       s  * (one - t)
#define l3       ( one - (r + s) ) * t
#define l4       r * t
#define l5       s * t

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	      b_[b_ld_ * i + 5] = l5;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;
	      b_[b_ld_ * 5 + i] = l5;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i],
		t = c_[c_ld_* 2  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	      b_[b_ld_ * i + 5] = l5;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i],
		t = c_[c_ld_* 2 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;
	      b_[b_ld_ * 5 + i] = l5;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3
#undef l4
#undef l5

  };

};

template<typename T>
struct bms_shape_lagrange_splz<1,WMESH_ELEMENT_HEXAHEDRON,T>
{

  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						     wmesh_int_p  	rw_)
  {
    static constexpr T one = static_cast<T>(1);
    static constexpr T eight = static_cast<T>(8);

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;


    	      
#define l0     (one - r)*(one - s)*(one - t) / eight
#define l1       (one + r)*(one - s)*(one - t) / eight
#define l2       (one + r)*(one + s)*(one - t) / eight
#define l3       (one - r)*(one + s)*(one - t) / eight
#define l4       (one - r)*(one - s)*(one + t) / eight
#define l5       (one + r)*(one - s)*(one + t) / eight
#define l6       (one + r)*(one + s)*(one + t) / eight
#define l7       (one - r)*(one + s)*(one + t) / eight

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	      b_[b_ld_ * i + 5] = l5;
	      b_[b_ld_ * i + 6] = l6;
	      b_[b_ld_ * i + 7] = l7;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;
	      b_[b_ld_ * 5 + i] = l5;
	      b_[b_ld_ * 6 + i] = l6;
	      b_[b_ld_ * 7 + i] = l7;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i],
		t = c_[c_ld_* 2  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	      b_[b_ld_ * i + 5] = l5;
	      b_[b_ld_ * i + 6] = l6;
	      b_[b_ld_ * i + 7] = l7;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i],
		t = c_[c_ld_* 2 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;
	      b_[b_ld_ * 5 + i] = l5;
	      b_[b_ld_ * 6 + i] = l6;
	      b_[b_ld_ * 7 + i] = l7;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3
#undef l4
#undef l5
#undef l6
#undef l7
  };
};

template<typename T>
struct bms_shape_lagrange_splz<1,WMESH_ELEMENT_PYRAMID,T>
{
  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						     wmesh_int_p  	rw_)
  {
    static constexpr T one = static_cast<T>(1);
    static constexpr T four = static_cast<T>(4);

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;

#define l0 under ? ( (one - t) * (one - r)*(one - s) / four) : 0
#define l1 under ? ( (one - t) * (one + r)*(one - s) / four) : 0
#define l2 under ? ( (one - t) * (one + r)*(one + s) / four) : 0
#define l3 under ? ( (one - t) * (one - r)*(one + s) / four) : 0
#define l4 t

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      const bool under = (t < 1.0);

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      const bool under = (t < 1.0);

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i],
		t = c_[c_ld_* 2  + i];

	      const bool under = (t < 1.0);

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i],
		t = c_[c_ld_* 2 + i];


	      const bool under = (t < 1.0);

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3
#undef l4
  };
};


template<wmesh_int_t 	degree_,
	 wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_shape_legendre_splz
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
		      wmesh_int_p   rw_);
};

template<wmesh_int_t ELEMENT_,typename T>
struct bms_shape_legendre_splz<0,ELEMENT_,T>
{
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
		      wmesh_int_p   rw_)
  {
    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;
    
    static constexpr T one = static_cast<T>(1);
    
#define l0       one
    
    switch(treat_case)
      {
	
      case 1:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      b_[b_ld_ * i + 0] = l0;	     
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      b_[b_ld_ * 0 + i] = l0;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      b_[b_ld_ * i + 0] = l0;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      b_[b_ld_ * 0 + i] = l0;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0

  }
};


template<typename T>
struct bms_shape_legendre_splz<1,WMESH_ELEMENT_EDGE,T>
{
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						  wmesh_int_p  	rw_)
  {
    
    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;

    static constexpr T one = static_cast<T>(1);
    static constexpr T two = static_cast<T>(2);

#define l0       ( one - r ) / two
#define l1       ( one + r ) / two

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1

  };
};

template<typename T>
struct bms_shape_legendre_splz<1,WMESH_ELEMENT_TRIANGLE,T>
{
  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						      wmesh_int_p  	rw_)
  {
    
    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;
    

    static constexpr T one = static_cast<T>(1);

#define l0       (one - (r+s))
#define l1       r
#define l2       s


    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2

  };
};

template<typename T>
struct bms_shape_legendre_splz<1,WMESH_ELEMENT_QUADRILATERAL,T>
{


  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
		      wmesh_int_p  	rw_)
  {

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;

        static constexpr T four = static_cast<T>(4);
    static constexpr T one = static_cast<T>(1);

#define l0       (one - r)*(one - s) / four
#define l1       (one + r)*(one + s) / four
#define l2       (one + r)*(one + s) / four
#define l3       (one - r)*(one + s) / four

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3

  };


};

template<typename T>
struct bms_shape_legendre_splz<1,WMESH_ELEMENT_TETRAHEDRON,T>
{


  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
							 wmesh_int_p  	rw_)
  {

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;
    
    static constexpr T one = static_cast<T>(1);
    
#define l0      ( one - (r + s + t) )
#define l1       r
#define l2       s
#define l3       t

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i],
		t = c_[c_ld_* 2  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i],
		t = c_[c_ld_* 2 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3

  };


};

template<typename T>
struct bms_shape_legendre_splz<1,WMESH_ELEMENT_WEDGE,T>
{

  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						     wmesh_int_p  	rw_)
  {
    static constexpr T one = static_cast<T>(1);

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;

#define l0     ( one - (r + s) ) * (one - t)
#define l1       r  * (one - t);
#define l2       s  * (one - t)
#define l3       ( one - (r + s) ) * t
#define l4       r * t
#define l5       s * t

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	      b_[b_ld_ * i + 5] = l5;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;
	      b_[b_ld_ * 5 + i] = l5;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i],
		t = c_[c_ld_* 2  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	      b_[b_ld_ * i + 5] = l5;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i],
		t = c_[c_ld_* 2 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;
	      b_[b_ld_ * 5 + i] = l5;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3
#undef l4
#undef l5

  };

};

template<typename T>
struct bms_shape_legendre_splz<1,WMESH_ELEMENT_HEXAHEDRON,T>
{

  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						     wmesh_int_p  	rw_)
  {
    static constexpr T one = static_cast<T>(1);
    static constexpr T eight = static_cast<T>(8);

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;


    	      
#define l0     (one - r)*(one - s)*(one - t) / eight
#define l1       (one + r)*(one - s)*(one - t) / eight
#define l2       (one + r)*(one + s)*(one - t) / eight
#define l3       (one - r)*(one + s)*(one - t) / eight
#define l4       (one - r)*(one - s)*(one + t) / eight
#define l5       (one + r)*(one - s)*(one + t) / eight
#define l6       (one + r)*(one + s)*(one + t) / eight
#define l7       (one - r)*(one + s)*(one + t) / eight

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	      b_[b_ld_ * i + 5] = l5;
	      b_[b_ld_ * i + 6] = l6;
	      b_[b_ld_ * i + 7] = l7;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;
	      b_[b_ld_ * 5 + i] = l5;
	      b_[b_ld_ * 6 + i] = l6;
	      b_[b_ld_ * 7 + i] = l7;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i],
		t = c_[c_ld_* 2  + i];

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	      b_[b_ld_ * i + 5] = l5;
	      b_[b_ld_ * i + 6] = l6;
	      b_[b_ld_ * i + 7] = l7;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i],
		t = c_[c_ld_* 2 + i];


	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;
	      b_[b_ld_ * 5 + i] = l5;
	      b_[b_ld_ * 6 + i] = l6;
	      b_[b_ld_ * 7 + i] = l7;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3
#undef l4
#undef l5
#undef l6
#undef l7
  };
};

template<typename T>
struct bms_shape_legendre_splz<1,WMESH_ELEMENT_PYRAMID,T>
{
  
  static inline  wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
						     wmesh_int_p  	rw_)
  {
    static constexpr T one = static_cast<T>(1);
    static constexpr T four = static_cast<T>(4);

    wmesh_int_t treat_case
      = ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 1
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_INTERLEAVE == c_storage_)) ? 2
      : ((WMESH_STORAGE_INTERLEAVE == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 3
      : ((WMESH_STORAGE_BLOCK      == b_storage_)&&(WMESH_STORAGE_BLOCK == c_storage_)) ? 4
      : 0;

#define l0 under ? ( (one - t) * (one - r)*(one - s) / four) : 0
#define l1 under ? ( (one - t) * (one + r)*(one - s) / four) : 0
#define l2 under ? ( (one - t) * (one + r)*(one + s) / four) : 0
#define l3 under ? ( (one - t) * (one - r)*(one + s) / four) : 0
#define l4 t

    switch(treat_case)
      {
	
      case 1:
	{
	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      const bool under = (t < 1.0);

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 2:
	{	  
	  for (wmesh_int_t i=0;i<c_n_;++i)
	    {
	      const T
		r = c_[c_ld_*i + 0],
		s = c_[c_ld_*i + 1],
		t = c_[c_ld_*i + 2];

	      const bool under = (t < 1.0);

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;

	    }
	  return WMESH_STATUS_SUCCESS;
	}

      case 3:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0  + i],
		s = c_[c_ld_* 1  + i],
		t = c_[c_ld_* 2  + i];

	      const bool under = (t < 1.0);

	      b_[b_ld_ * i + 0] = l0;
	      b_[b_ld_ * i + 1] = l1;
	      b_[b_ld_ * i + 2] = l2;
	      b_[b_ld_ * i + 3] = l3;
	      b_[b_ld_ * i + 4] = l4;
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	
      case 4:
	{
	  for (wmesh_int_t i=0;i<c_m_;++i)
	    {
	      const T
		r = c_[c_ld_* 0 + i],
		s = c_[c_ld_* 1 + i],
		t = c_[c_ld_* 2 + i];


	      const bool under = (t < 1.0);

	      b_[b_ld_ * 0 + i] = l0;
	      b_[b_ld_ * 1 + i] = l1;
	      b_[b_ld_ * 2 + i] = l2;
	      b_[b_ld_ * 3 + i] = l3;
	      b_[b_ld_ * 4 + i] = l4;

	    }
	  return WMESH_STATUS_SUCCESS;
	}
	      
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      }

    return WMESH_STATUS_INVALID_ENUM;
#undef l0
#undef l1
#undef l2
#undef l3
#undef l4
  };
};








template<wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_shape_lagrange
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	degree_,
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
		      wmesh_int_p   rw_)
  {
    
#define FORWARD					\
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
    
    switch(degree_)
      {
      case 0:
	{
	  return bms_shape_lagrange_splz<0,ELEMENT_,T>::eval(FORWARD);
	}
      case 1:
	{
	  return bms_shape_lagrange_splz<1,ELEMENT_,T>::eval(FORWARD);
	}
      default:
	{
	  return WMESH_STATUS_SUCCESS;
	}
      }

    return WMESH_STATUS_INVALID_ARGUMENT;
#undef FORWARD
    
  }

};




  template<wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_shape_legendre
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	degree_,
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
		      wmesh_int_p   rw_)
  {
    
#define FORWARD					\
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
    
    switch(degree_)
      {
      case 0:
	{
	  return bms_shape_legendre_splz<0,ELEMENT_,T>::eval(FORWARD);
	}
      case 1:
	{
	  return bms_shape_legendre_splz<1,ELEMENT_,T>::eval(FORWARD);
	}
      default:
	{
	  return WMESH_STATUS_SUCCESS;
	}
      }

    return WMESH_STATUS_INVALID_ARGUMENT;
#undef FORWARD
    
  }

};





  



  template<wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_shape_orthogonal
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	degree_,
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
		      wmesh_int_p   rw_)
  {
    
#define FORWARD					\
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
    
    switch(degree_)
      {
      case 0:
	{
	  return bms_shape_orthogonal_splz<0,ELEMENT_,T>::eval(FORWARD);
	}
      case 1:
	{
	  return bms_shape_orthogonal_splz<1,ELEMENT_,T>::eval(FORWARD);
	}
      default:
	{
	  return WMESH_STATUS_SUCCESS;
	}
      }

    return WMESH_STATUS_INVALID_ARGUMENT;
#undef FORWARD
    
  }

};



  
template<typename T>
struct 
bms_shape<WMESH_SHAPE_FAMILY_LAGRANGE,T>
{
  static inline  wmesh_status_t eval(wmesh_int_t 		element_,
			     wmesh_int_t 		degree_,
		      
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
			     wmesh_int_p		rw_)
  {

#define FORWARD					\
    degree_,					\
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
#define TREAT_CASE(_c) case _c: return bms_shape_lagrange<_c,T>::eval(FORWARD)

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
				 wmesh_int_p		rw_)
  {
      
#define FORWARD					\
    degree_,					\
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
#define TREAT_CASE(_c) case _c: return bms_shape_legendre<_c,T>::eval(FORWARD)
	  
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
				     wmesh_int_p		rw_)
  {
      
#define FORWARD					\
    degree_,					\
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
#define TREAT_CASE(_c) case _c: return bms_shape_orthogonal<_c,T>::eval(FORWARD)
	  
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
				  wmesh_int_p		rw_)
{
      
#define FORWARD					\
  element_,					\
    degree_,					\
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
      
  switch(element_)
    {
      
#define TREAT_CASE(_f) case _f: return bms_shape<_f,T>::eval(FORWARD)
      
      TREAT_CASE(WMESH_SHAPE_FAMILY_LAGRANGE);
      TREAT_CASE(WMESH_SHAPE_FAMILY_LEGENDRE);
      TREAT_CASE(WMESH_SHAPE_FAMILY_ORTHOGONAL);
      
#undef TREAT_CASE
      
    }
  
  return WMESH_STATUS_INVALID_ENUM;
  
#undef FORWARD

}


extern "C" wmesh_status_t bms_dshape(wmesh_int_t 		element_,
				     wmesh_int_t 		family_,
				     wmesh_int_t 		degree_,
				     
				     wmesh_int_t 		c_storage_,					  
				     wmesh_int_t 		c_m_,
				     wmesh_int_t 		c_n_,
				     const double * 			c_,
				     wmesh_int_t 		c_ld_,
					    
				     wmesh_int_t 		b_storage_,
				     wmesh_int_t 		b_m_,
				     wmesh_int_t 		b_n_,
				     double* 				b_,
				     wmesh_int_t 		b_ld_,
				     wmesh_int_t			iw_n_,
				     wmesh_int_p			iw_,
				     wmesh_int_t			rw_n_,
				     wmesh_int_p			rw_)
{    
  return bms_template_shape(element_,
			    family_,
			    degree_,
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



  

  
#if 0
template<typename T>
wmesh_status_t bms_template_shape_vandermonde(wmesh_int_t 	element_,
					      wmesh_int_t 	family_,
					      wmesh_int_t 	degree_,
					      
					      wmesh_int_t 	c_storage_,					  
					      wmesh_int_t 	c_m_,
					      wmesh_int_t 	c_n_,
					      const T * 	c_,
					      wmesh_int_t 	c_ld_,
					      
					      wmesh_int_t 	v_storage_,
					      wmesh_int_t 	v_m_,
					      wmesh_int_t 	v_n_,
					      T* 		v_,
					      wmesh_int_t 	v_ld_)
{
  return WMESH_STATUS_SUCCESS;
}


template<typename T>
wmesh_status_t bms_template_shape_nodal(wmesh_int_t 		vn_storage_,
					wmesh_int_t 		vn_m_,
					wmesh_int_t 		vn_n_,
					T* 			vn_,
					wmesh_int_t 		vn_ld_,

					wmesh_int_t 		vc_storage_,
					wmesh_int_t 		vc_m_,
					wmesh_int_t 		vc_n_,
					const T* 		vc_,
					wmesh_int_t 		vc_ld_,
					
					wmesh_int_t 		e_storage_,
					wmesh_int_t 		e_m_,
					wmesh_int_t 		e_n_,
					const T* 		e_,
					wmesh_int_t 		e_ld_)
{

  if ( (vn_storage_ == WMESH_STORAGE_INTERLEAVE) && (vc_storage_ == WMESH_STORAGE_INTERLEAVE) && (e_storage_ == WMESH_STORAGE_INTERLEAVE) )
    {
      wmesh_int_t info_lapack;
      wmesh_int_p perm 	= (wmesh_int_p)malloc(sizeof(wmesh_int_t)*ndofs);
      double* id 		= (double*)calloc(ndofs*ndofs,sizeof(double));
      for (wmesh_int_t i=0;i<ndofs;++i) id[i*ndofs+i] = 1.0;
  
      LAPACK_dgesv((wmesh_int_p)&ndofs,
		   (wmesh_int_p)&ev_n,
		   vn_,
		   (wmesh_int_p)&vn_ld_,
		   perm,
		   vc_,
		   (wmesh_int_p)&vc_ld_,
		   (wmesh_int_p)&info_lapack);
      WMESH_CHECK( 0 == info_lapack );      
    }
  else if ( (vn_storage_ == WMESH_STORAGE_BLOCK) && (vc_storage_ == WMESH_STORAGE_BLOCK) && (e_storage_ == WMESH_STORAGE_BLOCK) )
    {
      
    }
  else
    {

    }
  
}

#endif
