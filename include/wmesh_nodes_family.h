#ifndef WMESH_NODES_FAMILY_H
#define WMESH_NODES_FAMILY_H

#include "wmesh-types.h"

#define WMESH_NODES_FAMILY_BEZIER 	0
#define WMESH_NODES_FAMILY_LAGRANGE 	1
#define WMESH_NODES_FAMILY_GAUSSLOBATTO 2
#define WMESH_NODES_FAMILY_ALL 		3
 
#ifdef __cplusplus
extern "C"
{
#endif
  
  const char * wmesh_nodes_family_to_string(wmesh_int_t self_);
  
#ifdef __cplusplus
}
#endif

#endif
