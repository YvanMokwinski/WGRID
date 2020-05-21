
#include "wmesh.h"
#include "wmesh-status.h"
#include "bms.h"
#include "bms_traits.hpp"
#include "bms_templates.hpp"
#include <string.h>
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
  bms_elements_num_hyperfaces(wmesh_int_t 		topodim_,
			      wmesh_int_p 		num_hyperfaces_)
  {
    switch(topodim_)
      {
	
#ifdef TREAT_CASE
#error TREAT_CASE already defined
#endif
      
#define TREAT_CASE(_c) case _c: return bms_template_elements_num_hyperfaces<_c>(num_hyperfaces_)
      
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
