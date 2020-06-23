
#include "wmesh.h"
#include "wmesh-status.h"
#include "bms.h"
#include "bms_traits.hpp"
#include "bms_templates.hpp"
#include <string.h>



template<typename T>
wmesh_status_t bms_mirrored_local_coordinates_edge(wmesh_int_t 			signed_rotation_,
						  wmesh_int_t 			c_storage_,
						  wmesh_int_t 			c_m_,
						  wmesh_int_t 			c_n_,
						  const T * __restrict__ 	c_v_,
						  wmesh_int_t 			c_ld_,
						  wmesh_int_t 			x_storage_,
						  wmesh_int_t 			x_m_,
						  wmesh_int_t 			x_n_,
						  T * __restrict__ 		x_v_,
						  wmesh_int_t 			x_ld_)
{
  
  const wmesh_int_t x_len 	= (x_storage_ == WMESH_STORAGE_INTERLEAVE) ? x_n_ : x_m_;
  const wmesh_int_t c_len 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_n_ : c_m_;
  WMESH_CHECK( x_len == c_len );
  
  T * __restrict__ u 		= x_v_;
  const wmesh_int_t x_inc 	= (x_storage_ == WMESH_STORAGE_INTERLEAVE) ? x_ld_ : 1;
  
  const T * __restrict__ r 	= c_v_;
  const wmesh_int_t c_inc 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_ld_ : 1;
  switch(signed_rotation_)
    {
    case 1:
      {
	xcopy(&x_len,r,&c_inc,u,&x_inc);
	return WMESH_STATUS_SUCCESS;
      }
    case 2:
      {
	for (wmesh_int_t i=0;i<x_len;++i)
	  {	    
	    u[x_inc*i] = r[c_inc*(x_len-1-i)];
	  }
	return WMESH_STATUS_SUCCESS;
      }
    case -1:
      {
	for (wmesh_int_t i=0;i<x_len;++i)
	  {	    
	    u[x_inc*i] = r[c_inc*(x_len-1-i)];
	  }
	return WMESH_STATUS_SUCCESS;
      }
    case -2:
      {
	xcopy(&x_len,r,&c_inc,u,&x_inc);
	return WMESH_STATUS_SUCCESS;
      }
    }
  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
};


template<typename T>
wmesh_status_t bms_mirrored_local_coordinates_triangle(wmesh_int_t 			signed_rotation_,
						      wmesh_int_t 			c_storage_,
						      wmesh_int_t 			c_m_,
						      wmesh_int_t 			c_n_,
						      const T * __restrict__ 	c_v_,
						      wmesh_int_t 			c_ld_,
						      wmesh_int_t 			x_storage_,
						      wmesh_int_t 			x_m_,
						      wmesh_int_t 			x_n_,
						      T * __restrict__ 		x_v_,
						      wmesh_int_t 			x_ld_)
{
  
  const wmesh_int_t x_len 	= (x_storage_ == WMESH_STORAGE_INTERLEAVE) ? x_n_ : x_m_;
  const wmesh_int_t c_len 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_n_ : c_m_;
  WMESH_CHECK( x_len == c_len );
  
  T * __restrict__ u 		= x_v_;
  T * __restrict__ v 		= x_v_ + ( (x_storage_ == WMESH_STORAGE_INTERLEAVE) ? 1 : x_ld_ );
  const wmesh_int_t x_inc 	= (x_storage_ == WMESH_STORAGE_INTERLEAVE) ? x_ld_ : 1;
  
  const T * __restrict__ r 	= c_v_;
  const T * __restrict__ s 	= c_v_ + ( (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? 1 : c_ld_ );
  const wmesh_int_t c_inc 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_ld_ : 1;
  
  switch(signed_rotation_)
    {
    case 1:
      {
	xcopy(&x_len,r,&c_inc,u,&x_inc);
	xcopy(&x_len,s,&c_inc,v,&x_inc);
	return WMESH_STATUS_SUCCESS;
      }
    case 2:
      {
	// (1-r-s) * (0,0) + r*(1,0) + s*(0,1)
	// (1-u-v) * (0,1) + u*(0,0) + v*(1,0)
	// u = 1 - s - r, v = r
	for (wmesh_int_t i=0;i<x_len;++i)
	  {	    
	    u[x_inc*i] = static_cast<T>(1) - (r[c_inc*i]  + s[c_inc*i]);
	    v[x_inc*i] = r[c_inc*i];
	  }
	return WMESH_STATUS_SUCCESS;
      }
    case 3:
      {
	// (1-r-s) * (0,0) + r*(1,0) + s*(0,1)
	// (1-u-v) * (1,0) + u*(0,1) + v*(0,0)
	// u = s, v = 1 - s -r
	for (wmesh_int_t i=0;i<x_len;++i)
	  {	    
	    u[x_inc*i] = s[c_inc*i];
	    v[x_inc*i] = static_cast<T>(1) - (r[c_inc*i]  + s[c_inc*i]);
	  }
	return WMESH_STATUS_SUCCESS;
      }
    case -1:
      {
	// (1-r-s) * (0,0) + r*(1,0) + s*(0,1)
	// (1-u-v) * (0,0) + u*(0,1) + v*(1,0)
	// u = s, v = r
	for (wmesh_int_t i=0;i<x_len;++i)
	  {	    
	    u[x_inc*i] = s[c_inc*i];
	    v[x_inc*i] = r[c_inc*i];
	  }
	return WMESH_STATUS_SUCCESS;
      }

    case -2:
      {
	// (1-r-s) * (0,0) + r*(1,0) + s*(0,1)
	// (1-u-v) * (1,0) + u*(0,0) + v*(0,1)
	// u = 1 - r - s, v = s
	for (wmesh_int_t i=0;i<x_len;++i)
	  {	    
	    u[x_inc*i] = static_cast<T>(1) - (r[c_inc*i]  + s[c_inc*i]);
	    v[x_inc*i] = s[c_inc*i];
	  }
	return WMESH_STATUS_SUCCESS;
      }

    case -3:
      {
	// (1-r-s) * (0,0) + r*(1,0) + s*(0,1)
	// (1-u-v) * (0,1) + u*(1,0) + v*(0,0)
	// u = r, v = 1 -  r - s
	for (wmesh_int_t i=0;i<x_len;++i)
	  {	    
	    u[x_inc*i] = r[c_inc*i];
	    v[x_inc*i] = static_cast<T>(1) - (r[c_inc*i]  + s[c_inc*i]);
	  }
	return WMESH_STATUS_SUCCESS;
      }
    }
  
  //
  // Face coordinates.
  //
  // 1
  // (1-r-s) * (0,0) + r*(1,0) + s*(0,1)
  // (1-u-v) * (0,0) + u*(1,0) + v*(0,1)
  // u = r, v = s
  // 2
  // (1-r-s) * (0,0) + r*(1,0) + s*(0,1)
  // (1-u-v) * (0,1) + u*(0,0) + v*(1,0)
  // u = 1 - s -r, v = r
  // 3
  // (1-r-s) * (0,0) + r*(1,0) + s*(0,1)
  // (1-u-v) * (1,0) + u*(0,1) + v*(0,0)
  // u = s, v = 1 - s -r
  //	  
  // -1
  // (1-r-s) * (0,0) + r*(1,0) + s*(0,1)
  // (1-u-v) * (0,0) + u*(0,1) + v*(1,0)
  // u = s, v = r
  // -2
  // (1-r-s) * (0,0) + r*(1,0) + s*(0,1)
  // (1-u-v) * (1,0) + u*(0,0) + v*(0,1)
  // u = 1 - r - s, v = s
  // -3
  // (1-r-s) * (0,0) + r*(1,0) + s*(0,1)
  // (1-u-v) * (0,1) + u*(1,0) + v*(0,0)
  // u = r, v = 1 -  r - s
  return WMESH_STATUS_INVALID_ARGUMENT;
};

template<typename T>
wmesh_status_t bms_mirrored_local_coordinates_quadrilateral(wmesh_int_t 			signed_rotation_,
							   wmesh_int_t 			c_storage_,
							   wmesh_int_t 			c_m_,
							   wmesh_int_t 			c_n_,
							   const T * __restrict__ 	c_v_,
							   wmesh_int_t 			c_ld_,
							   wmesh_int_t 			x_storage_,
							   wmesh_int_t 			x_m_,
							   wmesh_int_t 			x_n_,
							   T * __restrict__ 		x_v_,
							   wmesh_int_t 			x_ld_)
{
  
  const wmesh_int_t x_len 	= (x_storage_ == WMESH_STORAGE_INTERLEAVE) ? x_n_ : x_m_;
  const wmesh_int_t c_len 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_n_ : c_m_;
  WMESH_CHECK( x_len == c_len );
  
  T * __restrict__ u 		= x_v_;
  T * __restrict__ v 		= x_v_ + ( (x_storage_ == WMESH_STORAGE_INTERLEAVE) ? 1 : x_ld_ );
  const wmesh_int_t x_inc 	= (x_storage_ == WMESH_STORAGE_INTERLEAVE) ? x_ld_ : 1;
  
  const T * __restrict__ r 	= c_v_;
  const T * __restrict__ s 	= c_v_ + ( (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? 1 : c_ld_ );
  const wmesh_int_t c_inc 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_ld_ : 1;

  static constexpr T mr1 = static_cast<T>(-1);
  switch(signed_rotation_)
    {
    case 1:
      {
	xcopy(&x_len,r,&c_inc,u,&x_inc);
	xcopy(&x_len,s,&c_inc,v,&x_inc);
		return WMESH_STATUS_SUCCESS;;
      }
      
    case 2:
      {
	// (1-r)(1-s)/4 (-1,-1) + (1+r)(1-s)/4 (1,-1) + (1+r)(1+s)/4 (1,1) + (1-r)(1+s)/2 (-1,1)
	// (1-u)(1-v)/4 (-1,1) + (1+u)(1-v)/4 (-1,-1) + (1+u)(1+v)/4 (1,-1) + (1-u)(1+v)/2 (1,1)
	// u = -s,  v = r
	xcopy(&x_len,s,&c_inc,u,&x_inc);
	xscal(&x_len,&mr1,u,&x_inc);
	xcopy(&x_len,r,&c_inc,v,&x_inc);
		return WMESH_STATUS_SUCCESS;;
      }
      
    case 3:
      {
	// (1-r)(1-s)/4 (-1,-1) + (1+r)(1-s)/4 (1,-1) + (1+r)(1+s)/4 (1,1) + (1-r)(1+s)/2 (-1,1)
	// (1-u)(1-v)/4 (1,1) + (1+u)(1-v)/4 (-1,1) + (1+u)(1+v)/4 (-1,-1) + (1-u)(1+v)/2 (1,-1)
	// u = -r,  v = -s
	xcopy(&x_len,r,&c_inc,u,&x_inc);
	xscal(&x_len,&mr1,u,&x_inc);
	xcopy(&x_len,s,&c_inc,v,&x_inc);
	xscal(&x_len,&mr1,v,&x_inc);
		return WMESH_STATUS_SUCCESS;;
      }
      
    case 4:
      {
	// (1-r)(1-s)/4 (-1,-1) + (1+r)(1-s)/4 (1,-1) + (1+r)(1+s)/4 (1,1) + (1-r)(1+s)/2 (-1,1)
	// (1-u)(1-v)/4 (1,-1) + (1+u)(1-v)/4 (1,1) + (1+u)(1+v)/4 (-1,1) + (1-u)(1+v)/2 (-1,-1)
	// u = s,  v = -r
	xcopy(&x_len,s,&c_inc,u,&x_inc);
      
	xcopy(&x_len,r,&c_inc,v,&x_inc);
	xscal(&x_len,&mr1,v,&x_inc);
		return WMESH_STATUS_SUCCESS;;
      }
      
    case -1:
      {
	// (1-r)(1-s)/4 (-1,-1) + (1+r)(1-s)/4 (1,-1) + (1+r)(1+s)/4 (1,1) + (1-r)(1+s)/2 (-1,1)
	// (1-u)(1-v)/4 (-1,-1) + (1+u)(1-v)/4 (-1,1) + (1+u)(1+v)/4 (1,1) + (1-u)(1+v)/2 (1,-1)
	// u = s, v = r
	xcopy(&x_len,s,&c_inc,u,&x_inc);
	xcopy(&x_len,r,&c_inc,v,&x_inc);
		return WMESH_STATUS_SUCCESS;;
      }

    case -2:
      {
	// (1-r)(1-s)/4 (-1,-1) + (1+r)(1-s)/4 (1,-1) + (1+r)(1+s)/4 (1,1) + (1-r)(1+s)/2 (-1,1)
	// (1-u)(1-v)/4 (1,-1) + (1+u)(1-v)/4 (-1,-1) + (1+u)(1+v)/4 (-1,1) + (1-u)(1+v)/2 (1,1)
	// u = -r, v = s
	xcopy(&x_len,r,&c_inc,u,&x_inc);
	xscal(&x_len,&mr1,u,&x_inc);
	xcopy(&x_len,s,&c_inc,v,&x_inc);
		return WMESH_STATUS_SUCCESS;;
      }

    case -3:
      {
	// (1-r)(1-s)/4 (-1,-1) + (1+r)(1-s)/4 (1,-1) + (1+r)(1+s)/4 (1,1) + (1-r)(1+s)/2 (-1,1)
	// (1-u)(1-v)/4 (1,1) + (1+u)(1-v)/4 (1,-1) + (1+u)(1+v)/4 (-1,-1) + (1-u)(1+v)/2 (-1,1)
	// u =-s, v = -r
	xcopy(&x_len,s,&c_inc,u,&x_inc);
	xscal(&x_len,&mr1,u,&x_inc);
	xcopy(&x_len,r,&c_inc,v,&x_inc);
	xscal(&x_len,&mr1,v,&x_inc);
		return WMESH_STATUS_SUCCESS;;
      }
    case -4:
      {
	// (1-r)(1-s)/4 (-1,-1) + (1+r)(1-s)/4 (1,-1) + (1+r)(1+s)/4 (1,1) + (1-r)(1+s)/2 (-1,1)
	// (1-u)(1-v)/4 (-1,1) + (1+u)(1-v)/4 (1,1) + (1+u)(1+v)/4 (1,-1) + (1-u)(1+v)/2 (-1,-1)
	// u = r, v = -s
	xcopy(&x_len,r,&c_inc,u,&x_inc);
	xcopy(&x_len,s,&c_inc,v,&x_inc);
	xscal(&x_len,&mr1,v,&x_inc);
	return WMESH_STATUS_SUCCESS;
      }
    }
  return WMESH_STATUS_INVALID_ARGUMENT;    
};



template<typename T>
wmesh_status_t bms_mirrored_local_coordinates(wmesh_int_t 	       		element_,
					      wmesh_int_t 			signed_rotation_,
					      wmesh_int_t 			c_storage_,
					      wmesh_int_t 			c_m_,
					      wmesh_int_t 			c_n_,
					      const T * __restrict__ 		c_v_,
					      wmesh_int_t 			c_ld_,
					      wmesh_int_t 			x_storage_,
					      wmesh_int_t 			x_m_,
					      wmesh_int_t 			x_n_,
					      T * __restrict__ 			x_v_,
					      wmesh_int_t 			x_ld_)
{
  switch(element_)
    {
    case WMESH_ELEMENT_EDGE:
      {
	return  bms_mirrored_local_coordinates_edge(signed_rotation_,
					       c_storage_,
					       c_m_,
					       c_n_,
					       c_v_,
					       c_ld_,
					       x_storage_,
					       x_m_,
					       x_n_,
					       x_v_,
					       x_ld_);
	
      }
    case WMESH_ELEMENT_TRIANGLE:
      {
       	return  bms_mirrored_local_coordinates_triangle(signed_rotation_,
					       c_storage_,
					       c_m_,
					       c_n_,
					       c_v_,
					       c_ld_,
					       x_storage_,
					       x_m_,
					       x_n_,
					       x_v_,
					       x_ld_);

      }
    case WMESH_ELEMENT_QUADRILATERAL:
      {
       	return  bms_mirrored_local_coordinates_quadrilateral(signed_rotation_,
							     c_storage_,
							     c_m_,
							     c_n_,
							     c_v_,
							     c_ld_,
							     x_storage_,
							     x_m_,
							     x_n_,
							     x_v_,
							     x_ld_);
      }
    case WMESH_ELEMENT_NODE:
    case WMESH_ELEMENT_TETRAHEDRON:
    case WMESH_ELEMENT_PYRAMID:
    case WMESH_ELEMENT_WEDGE:
    case WMESH_ELEMENT_HEXAHEDRON:
      {
	return WMESH_STATUS_INVALID_ARGUMENT;    
      }
    }
  return WMESH_STATUS_INVALID_ENUM;    
}


template
wmesh_status_t bms_mirrored_local_coordinates<float>	(wmesh_int_t 	       		element_,
							 wmesh_int_t 			signed_rotation_,
							 wmesh_int_t 			c_storage_,
							 wmesh_int_t 			c_m_,
							 wmesh_int_t 			c_n_,
							 const float * __restrict__ 	c_v_,
							 wmesh_int_t 			c_ld_,
							 wmesh_int_t 			x_storage_,
							 wmesh_int_t 			x_m_,
							 wmesh_int_t 			x_n_,
							 float * __restrict__ 		x_v_,
							 wmesh_int_t 			x_ld_);


template
wmesh_status_t bms_mirrored_local_coordinates<double>	(wmesh_int_t 	       		element_,
							 wmesh_int_t 			signed_rotation_,
							 wmesh_int_t 			c_storage_,
							 wmesh_int_t 			c_m_,
							 wmesh_int_t 			c_n_,
							 const double * __restrict__ 	c_v_,
							 wmesh_int_t 			c_ld_,
							 wmesh_int_t 			x_storage_,
							 wmesh_int_t 			x_m_,
							 wmesh_int_t 			x_n_,
							 double * __restrict__ 		x_v_,
							 wmesh_int_t 			x_ld_);


extern "C"
{
  
  wmesh_status_t bms_ndofs(wmesh_int_t element_,
			   wmesh_int_t d_,
			   wmesh_int_p ndofs_)
  {
#ifdef TREAT_CASE
#error TREAT_CASE is already defined
#else
#define TREAT_CASE(_f) case _f: ndofs_[0] = bms_template_ndofs<_f>(d_); return WMESH_STATUS_SUCCESS
#endif
    switch(element_)
    {
    TREAT_CASE(WMESH_ELEMENT_NODE);
    TREAT_CASE(WMESH_ELEMENT_EDGE);
    TREAT_CASE(WMESH_ELEMENT_TRIANGLE);
    TREAT_CASE(WMESH_ELEMENT_QUADRILATERAL);
    TREAT_CASE(WMESH_ELEMENT_TETRAHEDRON);
    TREAT_CASE(WMESH_ELEMENT_PYRAMID);
    TREAT_CASE(WMESH_ELEMENT_WEDGE);
    TREAT_CASE(WMESH_ELEMENT_HEXAHEDRON);
    
    }
  return WMESH_STATUS_INVALID_ENUM;
#undef TREAT_CASE
}
  wmesh_status_t bms_topodim2elements(wmesh_int_t 	topodim_,
				      wmesh_int_p 	num_elements_,
				      wmesh_int_p 	elements_)
  {
    WMESH_CHECK_POINTER(elements_);
    WMESH_CHECK_POINTER(num_elements_);
    switch(topodim_)
      {
	
#ifdef TREAT_CASE
#error TREAT_CASE already defined
#endif
	

#define TREAT_CASE(_c) case _c: return bms_template_topodim2elements<_c>(num_elements_,elements_)
      
	TREAT_CASE(WMESH_TOPODIM_NODE);
	TREAT_CASE(WMESH_TOPODIM_EDGE);
	TREAT_CASE(WMESH_TOPODIM_FACE);
	TREAT_CASE(WMESH_TOPODIM_VOLUME);

#undef TREAT_CASE
      
      }
  
    return WMESH_STATUS_INVALID_ENUM;    
  }

  wmesh_status_t bms_topodim2numtypes(wmesh_int_t 	topodim_,
				      wmesh_int_p 	ntypes_)
  {
    WMESH_CHECK_POINTER(ntypes_);
    switch(topodim_)
      {
	
#ifdef TREAT_CASE
#error TREAT_CASE already defined
#endif
      
#define TREAT_CASE(_c) case _c: return bms_template_topodim2numtypes<_c>(ntypes_)
      
	TREAT_CASE(WMESH_TOPODIM_VOLUME);
	TREAT_CASE(WMESH_TOPODIM_FACE);
	TREAT_CASE(WMESH_TOPODIM_EDGE);
	TREAT_CASE(WMESH_TOPODIM_NODE);
      
#undef TREAT_CASE
      
      }
    return WMESH_STATUS_INVALID_ENUM;
  }

  
  wmesh_status_t bms_element2topodim(wmesh_int_t 	element_,
				     wmesh_int_p 	topodim_)
  {
    WMESH_CHECK_POINTER(topodim_);
    switch(element_)
      {

#ifdef TREAT_CASE
#error TREAT_CASE already defined
#endif
      
#define TREAT_CASE(_c) case _c: return bms_template_element2topodim<_c>(topodim_)

	TREAT_CASE(WMESH_ELEMENT_NODE);
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
  }


    

  wmesh_status_t
  bms_elements_num_nodes(wmesh_int_t 		num_elements_,
			 const_wmesh_int_p 	elements_,
			 wmesh_int_p 		num_nodes_)
  {
    return bms_template_elements_num_entities<WMESH_ELEMENT_NODE>(num_elements_,
								  elements_,
								  num_nodes_);
  }

  
  wmesh_status_t
  bms_elements_num_edges(wmesh_int_t 		num_elements_,
			   const_wmesh_int_p 	elements_,
			   wmesh_int_p 		num_edges_)
  {
    return bms_template_elements_num_entities<WMESH_ELEMENT_EDGE>(num_elements_,
								  elements_,
								  num_edges_);
  }


  
  wmesh_status_t
  bms_element_facets(wmesh_int_t 		element_,
			 wmesh_int_p 		num_facets_,
			 wmesh_int_p 		facets_)
  {
    switch(element_)
      {
	
#ifdef TREAT_CASE
#error TREAT_CASE already defined
#endif
	
#define TREAT_CASE(_c) case _c: return bms_template_element_facets<_c>(num_facets_,facets_)      
	TREAT_CASE(WMESH_ELEMENT_EDGE);
	TREAT_CASE(WMESH_ELEMENT_TRIANGLE);
	TREAT_CASE(WMESH_ELEMENT_QUADRILATERAL);
	TREAT_CASE(WMESH_ELEMENT_TETRAHEDRON);
	TREAT_CASE(WMESH_ELEMENT_PYRAMID);
	TREAT_CASE(WMESH_ELEMENT_WEDGE);
	TREAT_CASE(WMESH_ELEMENT_HEXAHEDRON);      
#undef TREAT_CASE
	
      case WMESH_ELEMENT_NODE:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	}
      
      }
    return WMESH_STATUS_INVALID_ENUM;
  }

  wmesh_status_t
  bms_elements_num_facets(wmesh_int_t 		topodim_,
			      wmesh_int_p 		num_facets_)
  {
    switch(topodim_)
      {
	
#ifdef TREAT_CASE
#error TREAT_CASE already defined
#endif
      
#define TREAT_CASE(_c) case _c: return bms_template_elements_num_facets<_c>(num_facets_)
      
	TREAT_CASE(WMESH_TOPODIM_VOLUME);
	TREAT_CASE(WMESH_TOPODIM_FACE);
	TREAT_CASE(WMESH_TOPODIM_EDGE);
	TREAT_CASE(WMESH_TOPODIM_NODE);
      
#undef TREAT_CASE
      
      }
    return WMESH_STATUS_INVALID_ENUM;
  }


  wmesh_status_t
  bms_elements_num_entities(wmesh_int_t 		num_elements_,
			    const_wmesh_int_p 		elements_,
			    wmesh_int_t 		entity_,
			    wmesh_int_p 		num_entities_)
  {
    switch(entity_)
      {
	
#ifdef TREAT_CASE
#error TREAT_CASE already defined
#endif
      
#define TREAT_CASE(_c) case _c: return bms_template_elements_num_entities<_c>(num_elements_,elements_,num_entities_)
	
	TREAT_CASE(WMESH_ELEMENT_NODE);
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
  }

  
};
