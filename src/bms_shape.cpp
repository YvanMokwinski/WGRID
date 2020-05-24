
#include "bms.hpp"
#include "wmesh-math.hpp"
#include "wmesh-blas.h"
#include <iostream>
#ifndef NDEBUG
#include <iostream>
#endif





template<wmesh_int_t 	degree_,
	 wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_shape_orthogonal_splz
{  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
struct bms_shape_orthogonal_splz<0,ELEMENT_,T>
{
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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


  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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


  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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

  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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

  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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


  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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


  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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

  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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

  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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


  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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


  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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

  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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

  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
  
  static wmesh_status_t eval(wmesh_int_t 	c_storage_,					  
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
  static wmesh_status_t eval(wmesh_int_t 	degree_,
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
  static wmesh_status_t eval(wmesh_int_t 	degree_,
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
  static wmesh_status_t eval(wmesh_int_t 	degree_,
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



template<wmesh_int_t FAMILY_,typename T>
struct  bms_shape;

  
template<typename T>
struct 
bms_shape<WMESH_SHAPE_FAMILY_LAGRANGE,T>
{
   static wmesh_status_t eval(wmesh_int_t 		element_,
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
static     wmesh_status_t eval(wmesh_int_t 		element_,
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
    static wmesh_status_t eval(wmesh_int_t 		element_,
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
      element_,					\
	degree_,				\
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
#define TREAT_CASE(_f) case _f: return bms_shape<_f,T>::eval(FORWARD)	  
	  TREAT_CASE(WMESH_SHAPE_FAMILY_LAGRANGE);
	  TREAT_CASE(WMESH_SHAPE_FAMILY_LEGENDRE);
	  TREAT_CASE(WMESH_SHAPE_FAMILY_ORTHOGONAL);
#undef TREAT_CASE
	  
	}
return WMESH_STATUS_INVALID_ENUM;
#undef FORWARD


    }

    extern "C"    wmesh_status_t bms_dshape(wmesh_int_t 		element_,
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
      
      return bms_template_shape(
      element_,					\
	family_,\
	degree_,				\
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

template<typename T>
wmesh_status_t bms_template_shape(wmesh_int_t 	element_,
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

  
  status =  bms_template_shape_vandermonde(wmesh_int_t 	element_,
					   wmesh_int_t 	c_storage_,					  
					   wmesh_int_t 	c_m_,
					   wmesh_int_t 	c_n_,
					   const T * 	c_,
					   wmesh_int_t 	c_ld_,
					   
					   wmesh_int_t 	v_storage_,
					   wmesh_int_t 	v_m_,
					   wmesh_int_t 	v_n_,
					   T* 		v_,
					   wmesh_int_t 	v_ld_);
  

  status =  bms_template_shape_no(wmesh_int_t 	element_,
					   wmesh_int_t 	c_storage_,					  
					   wmesh_int_t 	c_m_,
					   wmesh_int_t 	c_n_,
					   const T * 	c_,
					   wmesh_int_t 	c_ld_,
					   
					   wmesh_int_t 	v_storage_,
					   wmesh_int_t 	v_m_,
					   wmesh_int_t 	v_n_,
					   T* 		v_,
					   wmesh_int_t 	v_ld_);
  

}

wmesh_status_t bms_vandermonde_monomial(wmesh_int_t 	dim_,
					wmesh_int_t 	degree_,

					wmesh_int_t 	c_storage_,					  
					wmesh_int_t 	c_m_,
					wmesh_int_t 	c_n_,
					const double * 	c_,
					wmesh_int_t 	c_ld_,
					
					wmesh_int_t 	v_storage_,
					wmesh_int_t 	v_m_,
					wmesh_int_t 	v_n_,
					double* 	v_,
					wmesh_int_t 	v_ld_)
{
  switch(dim_)
    {
    case 1:
      {
	wmesh_int_t c_inc = (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_ld_ : 1;
	switch(v_storage_)
	  {
	    
	  case WMESH_STORAGE_INTERLEAVE:
	    {	      
	      for (wmesh_int_t j=0;j<v_n_;++j)
		{
		  v_[v_ld_*j+0] = one;
		}
	      if (degree_ >=1)
		{
		  for (wmesh_int_t j=0;j<v_n_;++j)
		    {
		      v_[v_ld_*j+1] = c_[c_inc*j];
		    }
		  
		  for (wmesh_int_t d=2;d<=degree_;++d)
		    {
		      for (wmesh_int_t j=0;j<v_n_;++j)
			{
			  v_[v_ld_*j+d] = v_[v_ld_*j + (d-1)];
			}
		    } 	  
		}
	      return WMESH_STATUS_SUCCESS;
	    }
	    
	  case WMESH_STORAGE_BLOCK:
	    {	      
	      for (wmesh_int_t i=0;i<v_m_;++i)
		{
		  v_[v_ld_*0+i] = one;
		}
	      if (degree_ >=1)
		{
		  for (wmesh_int_t i=0;i<v_m_;++i)
		    {
		      v_[v_ld_*1+i] = c_[c_inc*i];
		    }		  
		  for (wmesh_int_t d=2;d<=degree_;++d)
		    {
		      for (wmesh_int_t i=0;i<v_m_;++i)
			{
			  v_[v_ld_*d+i] = v_[v_ld_*(d-1) + i];
			}
		    } 	  
		}
	      return WMESH_STATUS_SUCCESS;
	    }
	  }
	return WMESH_STATUS_INVALID_ENUM;
      }

    case 2:
      {
	wmesh_int_t c_inc = (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_ld_ : 1;
	switch(v_storage_)
	  {
	    
	  case WMESH_STORAGE_INTERLEAVE:
	    {	      
	      return WMESH_STATUS_SUCCESS;
	    }
	    
	  case WMESH_STORAGE_BLOCK:
	    {	      
	      for (wmesh_int_t i=0;i<v_m_;++i)
		{
		  v_[v_ld_*0+i] = one;
		}
	      if (degree_ >=1)
		{
		  for (wmesh_int_t i=0;i<v_m_;++i)
		    {
		      v_[v_ld_*1+i] = r[r_inc*i];
		    }
		  
		  for (wmesh_int_t i=0;i<v_m_;++i)
		    {
		      v_[v_ld_*2+i] = s[s_inc*i];
		    }
		  
		  for (wmesh_int_t d=2;d<=degree_;++d)
		    {
		      wmesh_int_t start1 = (d+1-1)*(d+2-1)/2;
		      wmesh_int_t start0 = (d+1-2)*(d+2-2)/2;
		      for (wmesh_int_t j=0;j<d;++j)
			{
			  for (wmesh_int_t i=0;i<v_m_;++i)
			    {
			      v_[v_ld_*(start1+j)+i] = r[rinc*i] * ( v_[v_ld_*(start0+j) + i] );
			    }
			}
		      
		      for (wmesh_int_t i=0;i<v_m_;++i)
			{
			  v_[v_ld_*(start1+d)+i] = s[sinc*i] * ( v_[v_ld_*(start0+d) + i] );
			}
		    }		  
		}
	      return WMESH_STATUS_SUCCESS;
	    }

	    
	  }
	return WMESH_STATUS_INVALID_ENUM;
      }


      
    case 3:
      {
	switch(v_storage_)
	  {
	    
	  case WMESH_STORAGE_INTERLEAVE:
	    {	      
	      return WMESH_STATUS_SUCCESS;
	    }
	    
	  case WMESH_STORAGE_BLOCK:
	    {	      
	      for (wmesh_int_t i=0;i<v_m_;++i)
		{
		  v_[v_ld_*0+i] = one;
		}
	      if (degree_ >=1)
		{
		  for (wmesh_int_t i=0;i<v_m_;++i)
		    {
		      v_[v_ld_*1+i] = r[r_inc*i];
		    }
		  
		  for (wmesh_int_t i=0;i<v_m_;++i)
		    {
		      v_[v_ld_*2+i] = s[s_inc*i];
		    }
		  
		  for (wmesh_int_t i=0;i<v_m_;++i)
		    {
		      v_[v_ld_*3+i] = t[t_inc*i];
		    }
		  
		  for (wmesh_int_t d=2;d<=degree_;++d)
		    {
		      wmesh_int_t start0 = ((d+1-2)*(d+2-2)*(d+3-2))/6;
		      wmesh_int_t start1 = ((d+1-1)*(d+2-1)*(d+3-1))/6;
		      // 1 r s t  ###### rr rs rt - ss st - tt
		      //           ##### rrr rrs rrt rss rst rtt -  sss sst stt - ttt -
		      //           ##### rrrr rrrs rrrt rrss rrst rrtt rsss rsst rstt rttt - ssss ssst sstt sttt - tttt
		      wmesh_int_t shift_s = ((d+1-2)*(d+2-2))/2;
		      wmesh_int_t nr      = ((d+1-1)*(d+2-1))/2;
		      wmesh_int_t ns      = d;
		      
		      for (wmesh_int_t j=0;j<start1-start0;++j)
			{
			  for (wmesh_int_t i=0;i<v_m_;++i)
			    {
			      v_[v_ld_*(start1+j)+i] = r[rinc*i] * v_[v_ld_*(start0+j) + i];
			    }
			}
		      
		      for (wmesh_int_t j=0;j<ns;++j)
			{
			  for (wmesh_int_t i=0;i<v_m_;++i)
			    {
			      v_[v_ld_*(start1+nr+j)+i] = s[sinc*i] * ( v_[v_ld_*(start0+shift_s+j) + i] );
			    }
			}

		      for (wmesh_int_t i=0;i<nt;++i)
			{
			  v_[v_ld_*(start1+nr+ns)+i] = t[tinc*i] * ( v_[v_ld_*(start1-1) + i] );
			}
		      
		    }		  
		}
	      return WMESH_STATUS_SUCCESS;
	    }

	    
	  }
	return WMESH_STATUS_INVALID_ENUM;
      }
      
    }
  return WMESH_STATUS_INVALID_ENUM;
}

wmesh_status_t bms_vandermonde_monomial(wmesh_int_t 	dim_,
					wmesh_int_t 	degree_,

					wmesh_int_t 	c_storage_,					  
					wmesh_int_t 	c_m_,
					wmesh_int_t 	c_n_,
					const double * 	c_,
					wmesh_int_t 	c_ld_,
					
					wmesh_int_t 	v_storage_,
					wmesh_int_t 	v_m_,
					wmesh_int_t 	v_n_,
					double* 	v_,
					wmesh_int_t 	v_ld_)
{
  switch(dim_)
    {
    case 1:
      {
	wmesh_int_t c_inc = (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_ld_ : 1;
	switch(v_storage_)
	  {
	    
	  case WMESH_STORAGE_INTERLEAVE:
	    {	      
	      for (wmesh_int_t j=0;j<v_n_;++j)
		{
		  v_[v_ld_*j+0] = one;
		}
	      if (degree_ >=1)
		{
		  for (wmesh_int_t j=0;j<v_n_;++j)
		    {
		      v_[v_ld_*j+1] = c_[c_inc*j];
		    }
		  
		  for (wmesh_int_t d=2;d<=degree_;++d)
		    {
		      for (wmesh_int_t j=0;j<v_n_;++j)
			{
			  v_[v_ld_*j+d] = v_[v_ld_*j + (d-1)];
			}
		    } 	  
		}
	      return WMESH_STATUS_SUCCESS;
	    }
	    
	  case WMESH_STORAGE_BLOCK:
	    {	      
	      for (wmesh_int_t i=0;i<v_m_;++i)
		{
		  v_[v_ld_*0+i] = one;
		}
	      if (degree_ >=1)
		{
		  for (wmesh_int_t i=0;i<v_m_;++i)
		    {
		      v_[v_ld_*1+i] = c_[c_inc*i];
		    }		  
		  for (wmesh_int_t d=2;d<=degree_;++d)
		    {
		      for (wmesh_int_t i=0;i<v_m_;++i)
			{
			  v_[v_ld_*d+i] = v_[v_ld_*(d-1) + i];
			}
		    } 	  
		}
	      return WMESH_STATUS_SUCCESS;
	    }
	  }
	return WMESH_STATUS_INVALID_ENUM;
      }

    case 2:
      {
	wmesh_int_t c_inc = (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_ld_ : 1;
	switch(v_storage_)
	  {
	    
	  case WMESH_STORAGE_INTERLEAVE:
	    {	      
	      return WMESH_STATUS_SUCCESS;
	    }
	    
	  case WMESH_STORAGE_BLOCK:
	    {	      
	      for (wmesh_int_t i=0;i<v_m_;++i)
		{
		  v_[v_ld_*0+i] = one;
		}
	      if (degree_ >=1)
		{
		  for (wmesh_int_t i=0;i<v_m_;++i)
		    {
		      v_[v_ld_*1+i] = r[r_inc*i];
		    }
		  
		  for (wmesh_int_t i=0;i<v_m_;++i)
		    {
		      v_[v_ld_*2+i] = s[s_inc*i];
		    }
		  
		  for (wmesh_int_t d=2;d<=degree_;++d)
		    {
		      wmesh_int_t start1 = (d+1-1)*(d+2-1)/2;
		      wmesh_int_t start0 = (d+1-2)*(d+2-2)/2;
		      for (wmesh_int_t j=0;j<d;++j)
			{
			  for (wmesh_int_t i=0;i<v_m_;++i)
			    {
			      v_[v_ld_*(start1+j)+i] = r[rinc*i] * ( v_[v_ld_*(start0+j) + i] );
			    }
			}
		      
		      for (wmesh_int_t i=0;i<v_m_;++i)
			{
			  v_[v_ld_*(start1+d)+i] = s[sinc*i] * ( v_[v_ld_*(start0+d) + i] );
			}
		    }		  
		}
	      return WMESH_STATUS_SUCCESS;
	    }

	    
	  }
	return WMESH_STATUS_INVALID_ENUM;
      }


      
    case 3:
      {
	switch(v_storage_)
	  {
	    
	  case WMESH_STORAGE_INTERLEAVE:
	    {	      
	      return WMESH_STATUS_SUCCESS;
	    }
	    
	  case WMESH_STORAGE_BLOCK:
	    {	      
	      for (wmesh_int_t i=0;i<v_m_;++i)
		{
		  v_[v_ld_*0+i] = one;
		}
	      if (degree_ >=1)
		{
		  for (wmesh_int_t i=0;i<v_m_;++i)
		    {
		      v_[v_ld_*1+i] = r[r_inc*i];
		    }
		  
		  for (wmesh_int_t i=0;i<v_m_;++i)
		    {
		      v_[v_ld_*2+i] = s[s_inc*i];
		    }
		  
		  for (wmesh_int_t i=0;i<v_m_;++i)
		    {
		      v_[v_ld_*3+i] = t[t_inc*i];
		    }
		  
		  for (wmesh_int_t d=2;d<=degree_;++d)
		    {
		      wmesh_int_t start0 = ((d+1-2)*(d+2-2)*(d+3-2))/6;
		      wmesh_int_t start1 = ((d+1-1)*(d+2-1)*(d+3-1))/6;
		      // 1 r s t  ###### rr rs rt - ss st - tt
		      //           ##### rrr rrs rrt rss rst rtt -  sss sst stt - ttt -
		      //           ##### rrrr rrrs rrrt rrss rrst rrtt rsss rsst rstt rttt - ssss ssst sstt sttt - tttt
		      wmesh_int_t shift_s = ((d+1-2)*(d+2-2))/2;
		      wmesh_int_t nr      = ((d+1-1)*(d+2-1))/2;
		      wmesh_int_t ns      = d;
		      
		      for (wmesh_int_t j=0;j<start1-start0;++j)
			{
			  for (wmesh_int_t i=0;i<v_m_;++i)
			    {
			      v_[v_ld_*(start1+j)+i] = r[rinc*i] * v_[v_ld_*(start0+j) + i];
			    }
			}
		      
		      for (wmesh_int_t j=0;j<ns;++j)
			{
			  for (wmesh_int_t i=0;i<v_m_;++i)
			    {
			      v_[v_ld_*(start1+nr+j)+i] = s[sinc*i] * ( v_[v_ld_*(start0+shift_s+j) + i] );
			    }
			}

		      for (wmesh_int_t i=0;i<nt;++i)
			{
			  v_[v_ld_*(start1+nr+ns)+i] = t[tinc*i] * ( v_[v_ld_*(start1-1) + i] );
			}
		      
		    }		  
		}
	      return WMESH_STATUS_SUCCESS;
	    }

	    
	  }
	return WMESH_STATUS_INVALID_ENUM;
      }
      
    }
  return WMESH_STATUS_INVALID_ENUM;
}

wmesh_status_t bms_vandermonde_tetrahedra(wmesh_int_t 		trans_,
					  wmesh_int_t 		d_,
					  
					  wmesh_int_t 		rst_storage_,					  
					  wmesh_int_t 		rst_m_,
					  wmesh_int_t 		rst_n_,
					  const double * 	rst_,
					  wmesh_int_t 		rst_ld_,
					  
					  double* 		v_,
					  wmesh_int_t 		v_ld_,
					  
					  wmesh_int_t 		work_n_,
					  double* 		work_);


#define solve_params(_type) wmesh_int_t*n,_type*jacM,_type*wr,_type*wi,_type*vl,_type*vr,_type*work,wmesh_int_t*work_n_
template<typename T>
void solve(solve_params(T));

template<>
void solve<double>(solve_params(double))
{
  static constexpr const  char trN[1] = {'N'};
  static constexpr const char trV[1] = {'V'};
  wmesh_int_t ldvl=1,info;
  LAPACK_dgeev(trN,
	trV,
	n,
	jacM,
	n, 
	wr, 
	wi, 
	vl, 
	&ldvl, 
	vr, 
	n, 
	work, 
	work_n_,
	&info);
}

template<>
void solve<float>(solve_params(float))
{
  static constexpr const  char trN[1] = {'N'};
  static constexpr const char trV[1] = {'V'};
  wmesh_int_t ldvl=1,info;
  LAPACK_sgeev(trN,
	trV,
	n,
	jacM,
	n, 
	wr, 
	wi, 
	vl, 
	&ldvl, 
	vr, 
	n, 
	work, 
	work_n_,
	&info);
}


    
static wmesh_status_t bms_cubature_legendre_buffer_size(wmesh_int_t nspl_,wmesh_int_p work_n_)
{
  work_n_[0] = (nspl_+2)*3 + 2*(nspl_+2)*(nspl_+2) + 4*(nspl_+1);
  return WMESH_STATUS_SUCCESS;
}

template<typename T>
wmesh_status_t bms_cubature_legendre(wmesh_int_t 		nspl_,
				     wmesh_int_t 		storagep_,
				     T*__restrict__ 		p_,
				     wmesh_int_t 		ldp_, 
				     T*__restrict__ 		w_,
				     wmesh_int_t 		incw_, 
				     wmesh_int_t 		work_n_,
				     T*__restrict__		work_)
{
  static constexpr T rhalf = static_cast<T>(0.5);
  static constexpr T r0 = static_cast<T>(0.0);
  static constexpr T r1 = static_cast<T>(1.0);
  static constexpr T r2 = static_cast<T>(2.0);
  static constexpr T r4 = static_cast<T>(4.0);
  
  //  if (WFE::storage_t::is_invalid(storagep_)) return WMESH_STATUS_INVALID_ENUM;
  
  if (nspl_ < 0)    	return WMESH_STATUS_INVALID_SIZE;
  if (!p_)      	return WMESH_STATUS_INVALID_POINTER;
  if (ldp_  < 1) 	return WMESH_STATUS_INVALID_SIZE;
  
  if (w_)
    {      
      if (incw_ < 1) 	return WMESH_STATUS_INVALID_SIZE;
    }
  wmesh_int_t required_work_n;
  wmesh_status_t status =  bms_cubature_legendre_buffer_size(nspl_,
							     &required_work_n);
  WMESH_STATUS_CHECK(status);
  if (required_work_n > work_n_)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
    }
  
  if (work_n_ < required_work_n ) return WMESH_STATUS_INVALID_SIZE;
  work_n_ -= (nspl_+2)*3 + 2*(nspl_+2)*(nspl_+2);
  
  bool p_is_block = (storagep_ == WMESH_STORAGE_BLOCK);
  if (p_is_block && ldp_ < nspl_) return WMESH_STATUS_INVALID_SIZE;  

  wmesh_int_t n    	= nspl_+1;
  wmesh_int_t nn 	= nspl_+2;
  T 
    * __restrict__ vl = NULL,
    * __restrict__ b    = &work_[0],
    * __restrict__ jacM = &work_[nn],
    * __restrict__ wr   = &work_[nn+nn*nn],
    * __restrict__ wi   = &work_[2*nn+nn*nn],
    * __restrict__ vr   = &work_[3*nn+nn*nn],
    * __restrict__ work  = &work_[3*nn+2*nn*nn];
  
#if 0
  T 
    * __restrict__ vl = NULL,
    * __restrict__ b    = &work_[nn],
    * __restrict__ jacM = &work_[2*nn],
    * __restrict__ wr   = &work_[2*nn+nn*nn],
    * __restrict__ wi   = &work_[2*nn+nn*nn+nn],
    * __restrict__ vr   = &work_[2*nn+nn*nn+2*nn],
    * __restrict__ work  = &work_[2*nn+2*nn*nn+2*nn];
#endif
  for (wmesh_int_t k=0;k<n;++k)
    {
      b[k]=(T)0.0;
    }
  
  b[0] = r2;
  b[1] = r0;
  for (wmesh_int_t k=2;k<n;++k)
    {
      b[k] = r1 / (r4 - r1/( ((T)(k-1))  * ((T)(k-1)) ) );
    }
  
  for (wmesh_int_t i=0;i<n;++i)
    {
      for (wmesh_int_t j=0;j<n;++j)
	{
	  jacM[i+j*n]=r0;
	}
    }
  
  for (wmesh_int_t i=0;i<n;++i)
    {
      jacM[i+i*n]=(T)0.0;
    }
  for (wmesh_int_t i=0;i<n-1;++i)
    {
      jacM[i+i*n+1]=wmesh_math<T>::xsqrt(b[1+i]);
    }
  for (wmesh_int_t i=1;i<n;++i)
    {
      jacM[i+i*n-1]=wmesh_math<T>::xsqrt(b[i]);
    }

  //
  // Eigen values / vectors decompositi
  //
  solve(&n,
	jacM,
	wr,
	wi,
	vl,
	vr,
	work,
	&work_n_);
 
  if (p_is_block)
    {
      for (wmesh_int_t i=0;i<n-1;++i) 
	{
	  p_[i] = wr[i];
	}
    }
  else
    {
      for (wmesh_int_t i=0;i<n-1;++i) 
	{
	  p_[i * ldp_] = wr[i];
	}	  
    }

  if (w_)
    {
      for (wmesh_int_t i=0;i<n-1;++i) 
	{
	  w_[i * incw_] = vr[1+i*n]*vr[1+i*n] * r2;
	}
    }

  for (int i = 0; i < (n-1)-1; i++)
    {
      int min_idx = i;
      for (int j = i+1; j < (n-1); j++)
	{
	  if (p_[j] < p_[min_idx])
	    {
	      min_idx = j;
	    }
	}
      
      { auto tmp = p_[min_idx]; p_[min_idx] = p_[i]; p_[i] = tmp; }

      if (w_)
	{
	  { auto tmp = w_[min_idx]; w_[min_idx] = w_[i]; w_[i] = tmp; }
	}
	  
    }
  
  //
  // Force symmetry.
  //
  if (nspl_ % 2 > 0)
    {
      p_[ (nspl_/2) *ldp_] = 0.0;
    }

  if (p_is_block)
    {
      for (int i=0;i<nspl_/2;++i)
	{
	  T tmp = (p_[i ] - p_[(nspl_-1-i) ])* rhalf;
	  p_[(nspl_-1-i)] = -tmp;
	  p_[i] = tmp;
	}
    }
  else
    {
      for (int i=0;i<nspl_/2;++i)
	{
	  T tmp = (p_[i *ldp_] - p_[(nspl_-1-i) *ldp_])*0.5;
	  p_[(nspl_-1-i) *ldp_] = -tmp;
	  p_[i *ldp_] = tmp;
	}
    }

  return WMESH_STATUS_SUCCESS;
}



wmesh_status_t
bms_cubature_num_nodes(wmesh_int_t	element_,
		       wmesh_int_t	family_,
		       wmesh_int_t	degree_,
		       wmesh_int_p	num_nodes_)
{
  wmesh_int_t num_nodes = (degree_%2 > 0) ? ((degree_ + 1) / 2) : ((degree_ + 2) / 2);
  num_nodes_[0] = 0;
  switch(family_)
    {
    case WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE:
      {
	switch(element_)
	  {
	  case WMESH_ELEMENT_NODE:
	    {
	      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	      break;
	    }
	    
	  case WMESH_ELEMENT_EDGE:
	    {
	      num_nodes_[0] = num_nodes;
	      return WMESH_STATUS_SUCCESS;
	    }
	    
	  case WMESH_ELEMENT_TRIANGLE:
	  case WMESH_ELEMENT_QUADRILATERAL:
	    {
	      num_nodes_[0] = num_nodes*num_nodes;
	      return WMESH_STATUS_SUCCESS;
	    }
	    
	  case WMESH_ELEMENT_TETRAHEDRON:
	  case WMESH_ELEMENT_PYRAMID:
	  case WMESH_ELEMENT_WEDGE:
	  case WMESH_ELEMENT_HEXAHEDRON:
	    {
	      num_nodes_[0] = num_nodes*num_nodes*num_nodes;
	      return WMESH_STATUS_SUCCESS;
	    }
	    
	  }
	return WMESH_STATUS_INVALID_ENUM;
      }
      
    case WMESH_CUBATURE_FAMILY_LOBATTO:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_NOT_IMPLEMENTED);	
      }

    }
  
  return WMESH_STATUS_INVALID_ENUM;  
}



template <typename T>
wmesh_status_t bms_template_cubature_transform(wmesh_int_t 	element_,
					       wmesh_int_t 	q_r_n_,
					       const T * 	q_r_v_,
					       wmesh_int_t 	q_r_inc_,
					       const T * 	q_w_v_,
					       wmesh_int_t 	q_w_inc_,
					       
					       wmesh_int_t	c_storage_,
					       wmesh_int_t	c_m_,
					       wmesh_int_t	c_n_,
					       T * 		c_v_,
					       wmesh_int_t	c_ld_,
					       
					       T * 		w_v_,
					       wmesh_int_t 	w_inc_)
#if 0
  wmesh_int_t		rst_storage_,
  wmesh_int_t		rst_m_,
  wmesh_int_t		rst_n_,
  const T*__restrict__ 	rst_v_,
  wmesh_int_t		rst_ld_;
#endif
#if 0  
T*  r_v = rst_v_ + (rst_storage_ == WMESH_STORAGE_BLOCK) ? rst_ld_ * 0 : 0;
T*  s_v = rst_v_ + (rst_storage_ == WMESH_STORAGE_BLOCK) ? rst_ld_ * 1 : 1;
T*  t_v = rst_v_ + (rst_storage_ == WMESH_STORAGE_BLOCK) ? rst_ld_ * 2 : 2;
#endif  

{

  WMESH_CHECK_POSITIVE(q_r_n_);
  WMESH_CHECK_POINTER(q_r_v_);
  WMESH_CHECK_POSITIVE(q_r_inc_);

  WMESH_CHECK_POINTER(q_w_v_);
  WMESH_CHECK_POSITIVE(q_w_inc_);

  WMESH_CHECK_POSITIVE(c_m_);
  WMESH_CHECK_POSITIVE(c_n_);
  WMESH_CHECK_POINTER(c_v_);
  WMESH_CHECK(c_ld_ >= c_m_);

  WMESH_CHECK_POINTER(w_v_);
  WMESH_CHECK_POSITIVE(w_inc_);
  
  static constexpr T s_one 	= static_cast<T>(1);
  static constexpr T s_two 	= static_cast<T>(2);
  static constexpr T s_three 	= static_cast<T>(3);
  static constexpr T s_four 	= static_cast<T>(4);

  const bool c_is_block = c_storage_ == WMESH_STORAGE_BLOCK;
  
  T*  r = c_v_ + ( c_is_block ? c_ld_ * 0 : 0);
  T*  s = c_v_ + ( c_is_block ? c_ld_ * 1 : 1);
  T*  t = c_v_ + ( c_is_block ? c_ld_ * 2 : 2);
  const wmesh_int_t inc = c_is_block  ? 1 : c_ld_;

  switch(element_)
    {
    case WMESH_ELEMENT_NODE:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	break;
      }
	    
    case WMESH_ELEMENT_EDGE:
      {		
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_TRIANGLE:
      {
	for (wmesh_int_t j=0;j<q_r_n_;++j)
	  {
	    for (wmesh_int_t i=0;i<q_r_n_;++i)
	      {
		r[ ( j * q_r_n_ + i ) * inc] 	= ( s_one + q_r_v_[i*q_r_inc_] ) * 0.5;
		s[ ( j * q_r_n_ + i ) * inc] 	= ( s_one - q_r_v_[i*q_r_inc_] ) * ( s_one - q_r_v_[j*q_r_inc_] ) * 0.25;		
		w_v_[ ( j * q_r_n_ + i ) * w_inc_] 	= ( s_one - q_r_v_[i*q_r_inc_] ) * q_w_v_[i*q_w_inc_] * q_w_v_[j*q_w_inc_] * 0.125;
	      }
	  }	
	return WMESH_STATUS_SUCCESS;
      }

    case WMESH_ELEMENT_QUADRILATERAL:
      {
	for (wmesh_int_t j=0;j<q_r_n_;++j)
	  {
	    for (wmesh_int_t i=0;i<q_r_n_;++i)
	      {
		r[ ( j * q_r_n_ + i ) * inc] 	= q_r_v_[i * q_r_inc_];
		s[ ( j * q_r_n_ + i ) * inc] 	= q_r_v_[j * q_r_inc_];
		w_v_[ ( j * q_r_n_ + i ) * w_inc_] 	= q_w_v_[i*q_w_inc_] * q_w_v_[j*q_w_inc_];
	      }
	  }
	
	return WMESH_STATUS_SUCCESS;
      }


    case WMESH_ELEMENT_TETRAHEDRON:
      {
	for (wmesh_int_t k=0;k<q_r_n_;++k)
	  {
	    for (wmesh_int_t j=0;j<q_r_n_;++j)
	      {
		for (wmesh_int_t i=0;i<q_r_n_;++i)
		  {
		    T u  = ( s_one + q_r_v_[i] ) / s_two;
		    T v  = ( s_one + q_r_v_[j] ) / s_two;
		    T w  = ( s_one + q_r_v_[k] ) / s_two;
		    
		    r[ ( k * q_r_n_ * q_r_n_ + j * q_r_n_ + i ) * inc] 	= u * v * w;
		    s[ ( k * q_r_n_ * q_r_n_ + j * q_r_n_ + i ) * inc] 	= u * v * (s_one - w);
		    t[ ( k * q_r_n_ * q_r_n_ + j * q_r_n_ + i ) * inc] 	= u * (s_one - u);
		    w_v_[ ( k * q_r_n_ * q_r_n_ + j * q_r_n_ + i ) * w_inc_] 	= (q_w_v_[i*q_w_inc_] *q_w_v_[i*q_w_inc_] * q_w_v_[j*q_w_inc_] * s_three) / s_four;
		  }
	      }
	  }		
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_HEXAHEDRON:
      {
	for (wmesh_int_t k=0;k<q_r_n_;++k)
	  {
	    for (wmesh_int_t j=0;j<q_r_n_;++j)
	      {
		for (wmesh_int_t i=0;i<q_r_n_;++i)
		  {
		    r[ ( q_r_n_ * q_r_n_ * k +  j * q_r_n_ + i ) * inc] 	= q_r_v_[i * q_r_inc_];
		    s[ ( q_r_n_ * q_r_n_ * k +  j * q_r_n_ + i ) * inc] 	= q_r_v_[j * q_r_inc_];
		    t[ ( q_r_n_ * q_r_n_ * k +  j * q_r_n_ + i ) * inc] 	= q_r_v_[k * q_r_inc_];
		    w_v_[ ( q_r_n_ * q_r_n_ * k +  j * q_r_n_ + i ) * w_inc_] 	= q_w_v_[i*q_w_inc_] * q_w_v_[j*w_inc_] * q_w_v_[k*q_w_inc_];
		  }
	      }
	  }
	return WMESH_STATUS_SUCCESS;
      }

    case WMESH_ELEMENT_WEDGE:
      {

	for (wmesh_int_t k=0;k<q_r_n_;++k)
	  {	    
	    for (wmesh_int_t j=0;j<q_r_n_;++j)
	      {
		for (wmesh_int_t i=0;i<q_r_n_;++i)
		  {
		    r[ ( k * q_r_n_* q_r_n_+j * q_r_n_ + i ) * inc] 	= ( s_one + q_r_v_[i*q_r_inc_] ) / s_two;
		    s[ ( k * q_r_n_* q_r_n_+j * q_r_n_ + i ) * inc] 	= ( s_one - q_r_v_[i*q_r_inc_] ) * ( s_one - q_r_v_[j*q_r_inc_] ) / s_four;
		    t[ ( k * q_r_n_* q_r_n_+j * q_r_n_ + i ) * inc] 	= ( s_one + q_r_v_[k*q_r_inc_] ) / s_two;		
		    w_v_[ ( k * q_r_n_* q_r_n_ + j * q_r_n_ + i ) * w_inc_] = ( s_one - q_r_v_[i*q_r_inc_] ) * q_w_v_[i*q_w_inc_] * q_w_v_[j*q_w_inc_] * 0.125 * q_w_v_[k *q_w_inc_];
		  }
	      }	    
	  }

	return WMESH_STATUS_SUCCESS;
      }

    case WMESH_ELEMENT_PYRAMID:
      {

	for (wmesh_int_t k=0;k<q_r_n_;++k)
	  {
	    T ww = q_r_v_[k*q_w_inc_];
	    T tt = (ww + s_one) / s_two;
	    
	    wmesh_int_t shiftk = k * q_r_n_* q_r_n_;	    
	    T scal = s_one - tt;
	    for (wmesh_int_t j=0;j<q_r_n_;++j)
	      {
		T vv = q_r_v_[j*q_w_inc_];
	    	wmesh_int_t shiftj = shiftk + j * q_r_n_;
		for (wmesh_int_t i=0;i<q_r_n_;++i)
		  {
		    T uu = q_r_v_[i*q_w_inc_];
		    wmesh_int_t at = shiftj + i;
		    
		    r[ at * inc] = uu * scal;
		    s[ at * inc] = vv * scal;
		    t[ at * inc] = tt;
		    
		    w_v_[ at * w_inc_] 	= q_w_v_[i*q_w_inc_] * q_w_v_[j*q_w_inc_] * q_w_v_[k *q_w_inc_] / s_two;
		  }
	      }	    
	  }
	return WMESH_STATUS_SUCCESS;
      }
      
    }
  
  return WMESH_STATUS_INVALID_ENUM;
}


template<typename T>
wmesh_status_t
bms_template_cubature(wmesh_int_t	element_,
		      wmesh_int_t	family_,
		      wmesh_int_t	n1d_,

		      wmesh_int_t	c_storage_,
		      wmesh_int_t	c_m_,
		      wmesh_int_t	c_n_,
		      T*__restrict__ 	c_v_,
		      wmesh_int_t	c_ld_,
	     
		      wmesh_int_t	w_n_,
		      T*__restrict__ 	w_v_,
		      wmesh_int_t	w_inc_,
	     
		      wmesh_int_t	rwork_n_,
		      T* __restrict__ 	rwork_)

{
  WMESH_CHECK_POSITIVE(c_m_);
  WMESH_CHECK_POSITIVE(c_n_);
  WMESH_CHECK_POINTER(c_v_);
  WMESH_CHECK(c_ld_ >= c_m_);

#ifndef NDEBUG
  std::cout << "// WMESH : bms_cubature : element = " << element_ << std::endl;
  std::cout << "// WMESH : bms_cubature : family  = " << family_ << std::endl;
  std::cout << "// WMESH : bms_cubature : n1d     = " << n1d_ << std::endl;
#endif  
  
  wmesh_status_t status;
  wmesh_int_t required_rwork_n;
  status = bms_cubature_buffer_size(element_,
				    family_,
				    n1d_,
				    &required_rwork_n);
  WMESH_STATUS_CHECK(status);

  if (required_rwork_n > rwork_n_)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
    }
  
  if (rwork_n_>0)
    {
      WMESH_CHECK_POINTER(rwork_);
    }
  
  switch(family_)
    {
    case WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE:
      {
	switch(element_)
	  {
	  case WMESH_ELEMENT_NODE:
	    {
	      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	      break;
	    }
	    
	  case WMESH_ELEMENT_EDGE:
	    {
	      status =  bms_cubature_legendre(n1d_,
					      c_storage_,
					      c_v_,
					      c_ld_, 
					      w_v_,					      
					      w_inc_,
					      rwork_n_,
					      rwork_);
	      WMESH_STATUS_CHECK(status);
	      return WMESH_STATUS_SUCCESS;
	    }
	    
	  case WMESH_ELEMENT_TRIANGLE:
	  case WMESH_ELEMENT_QUADRILATERAL:
	  case WMESH_ELEMENT_TETRAHEDRON:
	  case WMESH_ELEMENT_PYRAMID:
	  case WMESH_ELEMENT_WEDGE:
	  case WMESH_ELEMENT_HEXAHEDRON:
	    {
	      
	      T * p1d		=  rwork_;
	      T * w1d		=  rwork_ + n1d_;
	      rwork_n_ 	       -= n1d_*2;
	      rwork_ 	       += n1d_*2;
	      
	      status =  bms_cubature_legendre(n1d_,
					      WMESH_STORAGE_INTERLEAVE,
					      p1d,
					      1,
					      w1d,
					      1,
					      rwork_n_,
					      rwork_);
	      
	      WMESH_STATUS_CHECK(status);
	      status = bms_template_cubature_transform(element_,
						       n1d_,
						       p1d,
						       1,
						       w1d,
						       1,
						       c_storage_,
						       c_m_,
						       c_n_,
						       c_v_,
						       c_ld_,
						       w_v_,
						       w_inc_);
	      WMESH_STATUS_CHECK(status);
	      return WMESH_STATUS_SUCCESS;  
	    }
	  }	    
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_CUBATURE_FAMILY_LOBATTO:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_NOT_IMPLEMENTED);	
      }

    }
  
   return WMESH_STATUS_INVALID_ENUM;  
}

extern "C"
{

  
  wmesh_status_t bms_cubature_buffer_size(wmesh_int_t 	element_,
					  wmesh_int_t 	family_,
					  wmesh_int_t	n1d_,			
					  wmesh_int_p 	rwork_n_) 
  {
    WMESH_CHECK_POINTER(rwork_n_);

    rwork_n_[0] = 0;
    switch(family_)
      {
      case WMESH_NODES_FAMILY_LOBATTO:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_NOT_IMPLEMENTED);	
	}
      
      case WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE:
	{
	  wmesh_status_t status;	
	  status = bms_cubature_legendre_buffer_size(n1d_,
						     rwork_n_);
	  WMESH_STATUS_CHECK(status);
	
	  switch(element_)
	    {
	    
	    case WMESH_ELEMENT_NODE:
	      {
		WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	      }
	    
	    case WMESH_ELEMENT_EDGE:
	      {
		return WMESH_STATUS_SUCCESS;
	      }
	    
	    case WMESH_ELEMENT_TRIANGLE:
	    case WMESH_ELEMENT_QUADRILATERAL:
	      {
		rwork_n_[0] += n1d_*2;	      
		return WMESH_STATUS_SUCCESS;  
	      }
	    case WMESH_ELEMENT_TETRAHEDRON:
	    case WMESH_ELEMENT_PYRAMID:
	    case WMESH_ELEMENT_WEDGE:
	    case WMESH_ELEMENT_HEXAHEDRON:
	      {
		rwork_n_[0] += n1d_*2;	      
		return WMESH_STATUS_SUCCESS;  
	      }
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	

      }
  
    return WMESH_STATUS_INVALID_ENUM;  
  };

#if 1
  wmesh_status_t
  bms_cubature(wmesh_int_t 		element_,
	       wmesh_int_t 		family_,
	       wmesh_int_t		n1d_,
	       
	       wmesh_int_t		c_storage_,
	       wmesh_int_t		c_m_,
	       wmesh_int_t		c_n_,
	       double* 			c_v_,
	       wmesh_int_t		c_ld_,
	       
	       wmesh_int_t		w_n_,
	       double* 			w_v_,
	       wmesh_int_t		w_inc_,
	       
	       wmesh_int_t		rwork_n_,
	       double* __restrict__ 	rwork_)
  {
    return bms_template_cubature(element_,
				 family_,
				 n1d_,
				 c_storage_,
				 c_m_,
				 c_n_,
				 c_v_,
				 c_ld_,
				 w_n_,
				 w_v_,
				 w_inc_,
				 rwork_n_,
				 rwork_);    
  }
#endif
}

#if 0


//
//
// exact_degree = 2n-1;
// n mx_d
// 1 1
// 2 3
// 3 5
// 4 7
// 5 9
// 6 11

	
//
//
// Find min n to integrate d exactly
//
// if (d % 2 > 0) n = (d+1)/2 else (d+2)/2
//
//



  
#endif
#endif
