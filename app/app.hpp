#pragma once 
#include "wmesh.h"

#include <string.h>

//!
//! @brief element names for applications.
//!
static const char * s_wmesh_nodes_family_names[WMESH_NODES_FAMILY_ALL]
= {"lagrange",
   "gausslobatto"};

//!
//! @brief Convert a string to an nodesfamily.
//!
inline  wmesh_status_t app_str2nodesfamily(wmesh_str_t nodesfamily_name_,
					   wmesh_int_p nodesfamily_)
{
  WMESH_CHECK_POINTER(nodesfamily_);
  for (wmesh_int_t i=0;i<WMESH_NODES_FAMILY_ALL;++i)
    {            
      if (!strcmp(s_wmesh_nodes_family_names[i],
		  nodesfamily_name_))
	{
	  nodesfamily_[0] = i;
	  return WMESH_STATUS_SUCCESS;  
	}
    }
  return WMESH_STATUS_INVALID_ARGUMENT;  
}

//!
//! @brief Convert nodesfamily to a string.
//!
inline  wmesh_status_t app_nodesfamily2str(wmesh_int_t nodesfamily_,
					   wmesh_str_t nodesfamily_name_)
{
  strcpy(nodesfamily_name_,s_wmesh_nodes_family_names[nodesfamily_]);
  return WMESH_STATUS_SUCCESS;  
}

//!
//! @brief element names for applications.
//!
static const char * s_wmesh_element_names[WMESH_ELEMENT_ALL]
= {"node",
   "edge",
   "triangle",
   "quadrilateral",
   "tetrahedron",
   "pyramid",
   "prism",
   "hexahedron"};

//!
//! @brief Convert a string to an element.
//!
inline  wmesh_status_t app_str2element(wmesh_str_t element_name_,
				       wmesh_int_p element_)
{
  WMESH_CHECK_POINTER(element_);
  for (wmesh_int_t i=0;i<WMESH_ELEMENT_ALL;++i)
    {            
      if (!strcmp(s_wmesh_element_names[i],
		  element_name_))
	{
	  element_[0] = i;
	  return WMESH_STATUS_SUCCESS;  
	}
    }
  return WMESH_STATUS_INVALID_ARGUMENT;  
}
  
//!
//! @brief Convert element to a string.
//!
inline  wmesh_status_t app_element2str(wmesh_int_t element_,
				       wmesh_str_t element_name_)
{
  strcpy(element_name_,s_wmesh_element_names[element_]);
  return WMESH_STATUS_SUCCESS;  
}


#include "cmdline.hpp"
