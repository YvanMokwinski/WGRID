#include "bms.h"
#include "wmesh-status.h"
#include "wmesh-enums.h"
#include "wmesh-math.hpp"
#include <iostream>

//
// This needs to conform with the face ordering, otherwise it is going to fail.
//
template<typename T>
wmesh_status_t bms_ordering_normals(wmesh_int_t 	element_,
				    wmesh_int_t 	normals_storage_,
				    wmesh_int_t 	normals_m_,
				    wmesh_int_t 	normals_n_,
				    T*__restrict__	normals_,
				    wmesh_int_t 	normals_ld_)
{
  WMESH_CHECK_POINTER(normals_);

  const wmesh_int_t normals_dimension  = (WMESH_STORAGE_INTERLEAVE == normals_storage_) ? normals_m_ : normals_n_;
  const wmesh_int_t normals_num_facets = (WMESH_STORAGE_INTERLEAVE == normals_storage_) ? normals_n_ : normals_m_;

  T* normals_r = normals_;
  T* normals_s = (WMESH_STORAGE_INTERLEAVE == normals_storage_) ? (normals_ + 1) : (normals_ + normals_ld_);
  T* normals_t = (WMESH_STORAGE_INTERLEAVE == normals_storage_) ? (normals_ + 2) : (normals_ + normals_ld_ * 2);

  const wmesh_int_t normals_inc = (WMESH_STORAGE_INTERLEAVE == normals_storage_) ? normals_ld_ : 1;
  
  static constexpr T zero 	= static_cast<T>(0);
  static constexpr T one 	= static_cast<T>(1);
  const T invsqrt2 = one / wmesh_math<T>::xsqrt(2);
  const T invsqrt3 = one / wmesh_math<T>::xsqrt(3);
  switch(element_)
    {
      
    case WMESH_ELEMENT_NODE:
      {
	normals_[0] = zero;
	return WMESH_STATUS_SUCCESS;
      }

    case WMESH_ELEMENT_EDGE:
      {
	WMESH_CHECK(1 == normals_dimension);
	WMESH_CHECK(2 == normals_num_facets);

	normals_r[normals_inc * 0] = -one;
	normals_r[normals_inc * 1] = one;
	return WMESH_STATUS_SUCCESS;
      }

    case WMESH_ELEMENT_TRIANGLE:
      {
	WMESH_CHECK(2 == normals_dimension);
	WMESH_CHECK(3 == normals_num_facets);
	
	normals_r[normals_inc * 0] = zero;
	normals_s[normals_inc * 0] = -one;
	
	normals_r[normals_inc * 1] = invsqrt2;
	normals_s[normals_inc * 1] = invsqrt2;
	
	normals_r[normals_inc * 2] = -one;
	normals_s[normals_inc * 2] = zero;
	
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_QUADRILATERAL:
      {
	WMESH_CHECK(2 == normals_dimension);
	WMESH_CHECK(4 == normals_num_facets);
	
	normals_r[normals_inc * 0] = zero;
	normals_s[normals_inc * 0] = -one;

	normals_r[normals_inc * 1] = one;
	normals_s[normals_inc * 1] = zero;

	normals_r[normals_inc * 2] = zero;
	normals_s[normals_inc * 2] = one;

	normals_r[normals_inc * 3] = -one;
	normals_s[normals_inc * 3] = zero;
	
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_TETRAHEDRON:
      {
	WMESH_CHECK(3 == normals_dimension);
	WMESH_CHECK(4 == normals_num_facets);
	
	normals_r[normals_inc * 0] = invsqrt3;
	normals_s[normals_inc * 0] = invsqrt3;
	normals_t[normals_inc * 0] = invsqrt3;
	
	normals_r[normals_inc * 1] = -one;
	normals_s[normals_inc * 1] = zero;
	normals_t[normals_inc * 1] = zero;

	normals_r[normals_inc * 2] = zero;
	normals_s[normals_inc * 2] = -one;
	normals_t[normals_inc * 2] = zero;

	normals_r[normals_inc * 3] = zero;
	normals_s[normals_inc * 3] = zero;
	normals_t[normals_inc * 3] = -one;

	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_PYRAMID:
      {
	WMESH_CHECK(3 == normals_dimension);
	WMESH_CHECK(5 == normals_num_facets);
	
	normals_r[normals_inc * 0] = zero;
	normals_s[normals_inc * 0] = -one;
	normals_t[normals_inc * 0] = one;

	normals_r[normals_inc * 1] = one;
	normals_s[normals_inc * 1] = zero;
	normals_t[normals_inc * 1] = one;

	normals_r[normals_inc * 2] = zero;
	normals_s[normals_inc * 2] = one;
	normals_t[normals_inc * 2] = one;
	
	normals_r[normals_inc * 3] = -one;
	normals_s[normals_inc * 3] = zero;
	normals_t[normals_inc * 3] = one;
	
	normals_r[normals_inc * 4] = zero;
	normals_s[normals_inc * 4] = zero;
	normals_t[normals_inc * 4] = -one;

	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_WEDGE:
      {
	WMESH_CHECK(3 == normals_dimension);
	WMESH_CHECK(5 == normals_num_facets);

	normals_r[normals_inc * 0] = zero;
	normals_s[normals_inc * 0] = zero;
	normals_t[normals_inc * 0] = -one;

	normals_r[normals_inc * 1] = zero;
	normals_s[normals_inc * 1] = zero;
	normals_t[normals_inc * 1] = one;

	normals_r[normals_inc * 2] = -one;
	normals_s[normals_inc * 2] = zero;
	normals_t[normals_inc * 2] = zero;

	normals_r[normals_inc * 3] = invsqrt2;
	normals_s[normals_inc * 3] = invsqrt2;
	normals_t[normals_inc * 3] = zero;
	
	normals_r[normals_inc * 4] = zero;
	normals_s[normals_inc * 4] = -one;
	normals_t[normals_inc * 4] = zero;

	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_HEXAHEDRON:
      {
	WMESH_CHECK(3 == normals_dimension);
	WMESH_CHECK(6 == normals_num_facets);

	normals_r[normals_inc * 0] = zero;
	normals_s[normals_inc * 0] = zero;
	normals_t[normals_inc * 0] = -one;

	normals_r[normals_inc * 1] = zero;
	normals_s[normals_inc * 1] = zero;
	normals_t[normals_inc * 1] = one;

	normals_r[normals_inc * 2] = zero;
	normals_s[normals_inc * 2] = -one;
	normals_t[normals_inc * 2] = zero;

	normals_r[normals_inc * 3] = one;
	normals_s[normals_inc * 3] = zero;
	normals_t[normals_inc * 3] = zero;
	
	normals_r[normals_inc * 4] = zero;
	normals_s[normals_inc * 4] = one;
	normals_t[normals_inc * 4] = zero;
	
	normals_r[normals_inc * 5] = -one;
	normals_s[normals_inc * 5] = zero;
	normals_t[normals_inc * 5] = zero;
	
	return WMESH_STATUS_SUCCESS;
      }
      
    }
  return WMESH_STATUS_INVALID_ENUM;
}


template
wmesh_status_t bms_ordering_normals<double>(wmesh_int_t 	element_,
					    wmesh_int_t 	normals_storage_,
					    wmesh_int_t 	normals_m_,
					    wmesh_int_t 	normals_n_,
					    double*__restrict__	normals_,
					    wmesh_int_t 	normals_ld_);


template
wmesh_status_t bms_ordering_normals<float>(wmesh_int_t 	element_,
					   wmesh_int_t 	normals_storage_,
					   wmesh_int_t 	normals_m_,
					   wmesh_int_t 	normals_n_,
					   float*__restrict__ 	normals_,
					   wmesh_int_t 	normals_ld_);
