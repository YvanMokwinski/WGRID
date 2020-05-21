#include "wmesh_nodes.hpp"

#include "wmesh_utils.hpp"
#include "bms.h"


#if 0
template<typename T>
wmesh_status_t wmesh_nodes( WMESH_NODES_INTERFACE_PARAMS(T) );

template<>
wmesh_status_t wmesh_nodes<double>( WMESH_NODES_INTERFACE_PARAMS(double) )
{  
  wmesh_int_t element  = wmesh_nodes_element(self_);
  wmesh_int_t family = wmesh_nodes_family(self_);
  wmesh_int_t degree = wmesh_nodes_degree(self_);
  return  bms_dnodes(element,
		    family,
		    degree,
		    storagec_,
		    mc_,
		    nc_,
		    c_,
		    ldc_,
		    work_n_,
		    work_);  
}

template<>
wmesh_status_t wmesh_nodes<float>( WMESH_NODES_INTERFACE_PARAMS(float) )
{  
  wmesh_int_t family = wmesh_nodes_family(self_);
  wmesh_int_t degree = wmesh_nodes_degree(self_);
  wmesh_int_t element  = wmesh_nodes_element(self_);
  return  bms_snodes(element,
		    family,
		    degree,
		    storagec_,
		    mc_,
		    nc_,
		    c_,
		    ldc_,
		    work_n_,
		    work_);  
}
#endif

#if 0
template<typename T>
wmesh_status_t wmesh_nodes_triangle( WMESH_NODES_INTERFACE_PARAMS(T) )
{
  return WMESH_StATUS_SUCCESS;
}
template<typename T>
wmesh_status_t wmesh_nodes_quadrilateral( WMESH_NODES_INTERFACE_PARAMS(T) )
{
  return WMESH_StATUS_SUCCESS;
}
template<typename T>
wmesh_status_t wmesh_nodes_tetrahedron( WMESH_NODES_INTERFACE_PARAMS(T) )
{
  return WMESH_StATUS_SUCCESS;
}
template<typename T>
wmesh_status_t wmesh_nodes_pyramid( WMESH_NODES_INTERFACE_PARAMS(T) )
{
  return WMESH_StATUS_SUCCESS;
}

template<typename T>
wmesh_status_t wmesh_nodes_wedge( WMESH_NODES_INTERFACE_PARAMS(T) )
{
  return WMESH_StATUS_SUCCESS;
}

template<typename T>
wmesh_status_t wmesh_nodes_hexahedron( WMESH_NODES_INTERFACE_PARAMS(T) )
{
  return WMESH_StATUS_SUCCESS;
}

template<typename T>
wmesh_status_t wmesh_nodes( WMESH_NODES_INTERFACE_PARAMS(T) )
{
  wmesh_status_t status;
  wmesh_int_t element = wmesh_nodes_element(self_);
  switch(element)
    {
    case WMESH_ELEMENT_NODE:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);	
      }
    case WMESH_ELEMENT_EDGE:
      {
	return wmesh_nodes_edge(WMESH_NODES_INTERFACE_FORWARD);
      }
    case WMESH_ELEMENT_TRIANGLE:
      {
	return wmesh_nodes_triangle(WMESH_NODES_INTERFACE_FORWARD);
      }
    case WMESH_ELEMENT_QUADRILATERAL:
      {
	return wmesh_nodes_quadrilateral(WMESH_NODES_INTERFACE_FORWARD);
      }
    case WMESH_ELEMENT_TETRAHEDRON:
      {
	return wmesh_nodes_tetrahedron(WMESH_NODES_INTERFACE_FORWARD);
      }
    case WMESH_ELEMENT_PYRAMID:
      {
	return wmesh_nodes_pyramid(WMESH_NODES_INTERFACE_FORWARD);
      }
    case WMESH_ELEMENT_WEDGE:
      {
	return wmesh_nodes_wedge(WMESH_NODES_INTERFACE_FORWARD);
      }
    case WMESH_ELEMENT_HEXAHEDRON:
      {
	return wmesh_nodes_hexahedron(WMESH_NODES_INTERFACE_FORWARD);
      }
    }
  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
}
#endif

extern "C"
{
  static constexpr wmesh_nodes_t s_zero = 0;
  static constexpr wmesh_nodes_t s_nbits = sizeof(wmesh_nodes_t)*8;

  static constexpr wmesh_nodes_t s_nbits_element = 3;
  static constexpr wmesh_nodes_t s_nbits_family  = 3;

#if 0
  wmesh_int_t wmesh_nodes_degree(const wmesh_nodes_t *  self_)
  {
    static constexpr wmesh_nodes_t z = ((~s_zero) >> (s_nbits-3))<<(s_nbits_element + s_nbits_family);
    return (z & self_[0]) >> (s_nbits_element + s_nbits_family);
  };
  
  wmesh_int_t wmesh_nodes_family(const wmesh_nodes_t *  self_)
  {
    static constexpr wmesh_nodes_t z = ((~s_zero) >> (s_nbits-s_nbits_family)) << s_nbits_element;
    return (z & self_[0]) >> s_nbits_element;
  };
  
  wmesh_int_t wmesh_nodes_element(const wmesh_nodes_t *  self_)
  {
    static constexpr wmesh_nodes_t z = (~s_zero) >> (s_nbits - s_nbits_element);
    return z & self_[0]; 
  };
  
  wmesh_int_t wmesh_nodes_ndofs(const wmesh_nodes_t *  self_)
  {
    wmesh_int_t    ndofs;    
    wmesh_status_t status = wfe_ndofs(wmesh_nodes_element(self_),
				      wmesh_nodes_degree(self_),
				      &ndofs);
    WMESH_STATUS_CHECK(status);
    return ndofs;
  };
  
  wmesh_int_t wmesh_nodes_def(wmesh_nodes_t * 	self_,
			      wmesh_int_t 	element_,
			      wmesh_int_t 	family_,
			      wmesh_int_t 	degree_)
  {
    WMESH_CHECK_POINTER(self_);
    WMESH_CHECK(family_  < 8);
    WMESH_CHECK(element_ < 8);     
    self_[0] = element_ + (family_ << s_nbits_family) + (degree_ << (s_nbits_family + s_nbits_element));
    return WMESH_STATUS_SUCCESS;
  };
#endif
#if 0
  wmesh_status_t wmesh_snodes(WMESH_NODES_INTERFACE_PARAMS(float))
  {
    return wmesh_nodes(WMESH_NODES_INTERFACE_FORWARD);
  }
  
  wmesh_status_t wmesh_dnodes(WMESH_NODES_INTERFACE_PARAMS(double))
  {
    return wmesh_nodes(WMESH_NODES_INTERFACE_FORWARD);
  }
#endif
}

