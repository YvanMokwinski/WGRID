#pragma once
#include "bms_traits.hpp"

template<wmesh_int_t TOPODIM_>
inline wmesh_status_t bms_template_topodim2numtypes(wmesh_int_p 	ntypes_)
{
  WMESH_CHECK_POINTER(ntypes_);
  ntypes_[0] = bms_traits_topodim<TOPODIM_>::s_ntypes;
  return WMESH_STATUS_SUCCESS;
}

template<wmesh_int_t TOPODIM_>
wmesh_status_t bms_template_topodim2elements(wmesh_int_p 	num_elements_,
					     wmesh_int_p 	elements_)
{
  num_elements_[0] = bms_traits_topodim<TOPODIM_>::s_ntypes;
  for (wmesh_int_t i=0;i<bms_traits_topodim<TOPODIM_>::s_ntypes;++i)
    {
      elements_[i] = bms_traits_topodim<TOPODIM_>::s_ntypes + i;
    }
  return WMESH_STATUS_SUCCESS;
}

template<wmesh_int_t ELEMENT_>
inline wmesh_status_t bms_template_element2topodim(wmesh_int_p 	topodim_)
{
  WMESH_CHECK_POINTER(topodim_);
  topodim_[0] = bms_traits_element<ELEMENT_>::s_topodim;
  return WMESH_STATUS_SUCCESS;
}

template<wmesh_int_t TOPODIM_>
inline wmesh_status_t bms_template_elements_num_hyperfaces(wmesh_int_p 	num_hyperfaces_);

template<>
inline wmesh_status_t bms_template_elements_num_hyperfaces<WMESH_TOPODIM_VOLUME>(wmesh_int_p 	num_hyperfaces_)
{
  WMESH_CHECK_POINTER(num_hyperfaces_);
  num_hyperfaces_[0] = 4;
  num_hyperfaces_[1] = 5;
  num_hyperfaces_[2] = 5;
  num_hyperfaces_[3] = 6;
  return WMESH_STATUS_SUCCESS;
}

template<>
inline wmesh_status_t bms_template_elements_num_hyperfaces<WMESH_TOPODIM_FACE>(wmesh_int_p 	num_hyperfaces_)
{
  WMESH_CHECK_POINTER(num_hyperfaces_);
  num_hyperfaces_[0] = 3;
  num_hyperfaces_[1] = 4;
  return WMESH_STATUS_SUCCESS;
}


template<>
inline wmesh_status_t bms_template_elements_num_hyperfaces<WMESH_TOPODIM_EDGE>(wmesh_int_p 	num_hyperfaces_)
{
  WMESH_CHECK_POINTER(num_hyperfaces_);
  num_hyperfaces_[0] = 2;
  return WMESH_STATUS_SUCCESS;
}

template<>
inline wmesh_status_t bms_template_elements_num_hyperfaces<WMESH_TOPODIM_NODE>(wmesh_int_p 	num_hyperfaces_)
{
   WMESH_CHECK_POINTER(num_hyperfaces_);
  num_hyperfaces_[0] = 0;
  return WMESH_STATUS_SUCCESS;
}

template<wmesh_int_t ENTITY_>
inline wmesh_status_t bms_template_elements_num_entities(wmesh_int_t 		num_elements_,
							 const_wmesh_int_p 	elements_,
							 wmesh_int_p 		num_entities_)
{  
  static wmesh_int_t s_num_entities[8][8] = { {1,2,3,4,4,5,6,8},
					      {0,1,3,4,6,8,9,12},
					      {0,0,1,0,4,4,2,0},
					      {0,0,0,1,0,1,4,6},
					      {0,0,0,0,1,0,0,0},
					      {0,0,0,0,0,1,0,0},
					      {0,0,0,0,0,0,1,0},
					      {0,0,0,0,0,0,0,0} };
  for (wmesh_int_t i=0;i<num_elements_;++i)
    {
      num_entities_[i] = s_num_entities[ENTITY_][elements_[i]];
    }
  return WMESH_STATUS_SUCCESS;
}

template <wmesh_int_t ELEMENT> inline wmesh_int_t bms_template_ndofs(wmesh_int_t d_){ return -1; }

template <> inline wmesh_int_t bms_template_ndofs<WMESH_ELEMENT_EDGE>(wmesh_int_t d_) { return d_+1; }
template <> inline wmesh_int_t bms_template_ndofs<WMESH_ELEMENT_TRIANGLE>(wmesh_int_t d_) { return ((d_+1)*(d_+2)) / 2; }
template <> inline wmesh_int_t bms_template_ndofs<WMESH_ELEMENT_QUADRILATERAL>(wmesh_int_t d_) { return ((d_+1)*(d_+1)); }
template <> inline wmesh_int_t bms_template_ndofs<WMESH_ELEMENT_TETRAHEDRON>(wmesh_int_t d_) { return ((d_+1)*(d_+2)*(d_+3)) / 6; }
template <> inline wmesh_int_t bms_template_ndofs<WMESH_ELEMENT_PYRAMID>(wmesh_int_t d_) { return ((d_+2)*(d_+1)*(2*d_+3)) / 6; }
template <> inline wmesh_int_t bms_template_ndofs<WMESH_ELEMENT_WEDGE>(wmesh_int_t d_) { return ((d_+1)*(d_+1)*(d_+2)) / 2; }
template <> inline wmesh_int_t bms_template_ndofs<WMESH_ELEMENT_HEXAHEDRON>(wmesh_int_t d_) { return (d_+1)*(d_+1)*(d_+1); }


template <wmesh_int_t ELEMENT> inline wmesh_int_t bms_template_ndofs_interior(wmesh_int_t d_){ return -1;}

template <> inline wmesh_int_t bms_template_ndofs_interior<WMESH_ELEMENT_EDGE>(wmesh_int_t d_) { return (d_>0) ?d_-1 : 1; }
template <> inline wmesh_int_t bms_template_ndofs_interior<WMESH_ELEMENT_TRIANGLE>(wmesh_int_t d_) { return (d_>0) ?((d_-1)*(d_-2)) / 2 : 1; }
template <> inline wmesh_int_t bms_template_ndofs_interior<WMESH_ELEMENT_QUADRILATERAL>(wmesh_int_t d_) { return (d_>0) ?((d_-1)*(d_-1)) : 1; }
template <> inline wmesh_int_t bms_template_ndofs_interior<WMESH_ELEMENT_TETRAHEDRON>(wmesh_int_t d_) { return (d_>0) ?((d_-1)*(d_-2)*(d_-3)) / 6 : 1; }
template <> inline wmesh_int_t bms_template_ndofs_interior<WMESH_ELEMENT_PYRAMID>(wmesh_int_t d_) { return (d_>0) ?((d_-2)*(d_-1)*(2*d_-3)) / 6 : 1; }
template <> inline wmesh_int_t bms_template_ndofs_interior<WMESH_ELEMENT_WEDGE>(wmesh_int_t d_) { return (d_>0) ?((d_-1)*(d_-1)*(d_-2)) / 2 : 1; }
template <> inline wmesh_int_t bms_template_ndofs_interior<WMESH_ELEMENT_HEXAHEDRON>(wmesh_int_t d_) { return (d_>0) ? (d_-1)*(d_-1)*(d_-1) : 1; }






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
				  T * rw_);
