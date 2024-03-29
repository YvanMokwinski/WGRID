#ifndef WMESH_ENUMS_H
#define WMESH_ENUMS_H

#define WMESH_TOPODIM_NODE 			0
#define WMESH_TOPODIM_EDGE 			1
#define WMESH_TOPODIM_FACE 			2
#define WMESH_TOPODIM_VOLUME 			3

#define WMESH_ELEMENT_NODE 			0
#define WMESH_ELEMENT_EDGE 			1
#define WMESH_ELEMENT_TRIANGLE 			2
#define WMESH_ELEMENT_QUADRILATERAL 		3
#define WMESH_ELEMENT_TETRAHEDRON 		4
#define WMESH_ELEMENT_PYRAMID 			5
#define WMESH_ELEMENT_WEDGE 			6
#define WMESH_ELEMENT_HEXAHEDRON		7
#define WMESH_ELEMENT_ALL			8

#define WMESH_STORAGE_BLOCK             	0
#define WMESH_STORAGE_INTERLEAVE        	1

#define WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE   	0
#define WMESH_CUBATURE_FAMILY_GAUSSLOBATTO 	1
#define WMESH_CUBATURE_FAMILY_ALL 		2

#define WMESH_NODES_FAMILY_LAGRANGE 		0
#define WMESH_NODES_FAMILY_GAUSSLOBATTO 	1
#define WMESH_NODES_FAMILY_ALL 			2

#define WMESH_SHAPE_FAMILY_LAGRANGE   		0
#define WMESH_SHAPE_FAMILY_LEGENDRE 		1
#define WMESH_SHAPE_FAMILY_ORTHOGONAL 		2
#define WMESH_SHAPE_FAMILY_ALL 			3

#endif
