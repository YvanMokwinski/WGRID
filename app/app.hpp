#pragma once 
#include "wmesh.h"
#include <string.h>

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
