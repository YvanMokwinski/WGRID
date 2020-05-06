#pragma once

#include "wmesh.h"
#include "wmesh-types.hpp"

#define WMESH_ELEMENT_NODE 		0
#define WMESH_ELEMENT_EDGE 		1
#define WMESH_ELEMENT_TRIANGLE 		2
#define WMESH_ELEMENT_QUADRILATERAL 	3
#define WMESH_ELEMENT_TETRAHEDRON 	4
#define WMESH_ELEMENT_PYRAMID 		5
#define WMESH_ELEMENT_WEDGE 		6
#define WMESH_ELEMENT_HEXAHEDRON	7

#define WMESH_ELEMENT_ALL		8


extern "C"
{
  
  struct wmesh_t
  {

    //!
    //! @brief Graph Cells-2-Nodes.
    //!
    wmesh_int_sparsemat_t 	m_c2n;

    //!
    //! @brief Graph Cells-2-Edges.
    //!
    wmesh_int_sparsemat_t 	m_c2e;    

    //!
    //! @brief Graph Cells-2-Faces.
    //!
    wmesh_int_sparsemat_t 	m_c2f;

    //!
    //! @brief Graph Cells-2-Faces-Triangles.
    //!
    wmesh_int_sparsemat_t 	m_c2f_t;
    
    //!
    //! @brief Graph Cells-2-Faces-Quadrilaterals.
    //!
    wmesh_int_sparsemat_t 	m_c2f_q;

    //!
    //! @brief Graph Cells-2-Cells (through faces).
    //!
    wmesh_int_sparsemat_t 	m_c2c;
    //!
    //! @brief Graph Cells-2-Cells-(Through)Triangles
    //!
    wmesh_int_sparsemat_t 	m_c2c_t;
    
    //!
    //! @brief Graph Cells-2-Cells-(Through)Quadrilaterals.
    //!
    wmesh_int_sparsemat_t 	m_c2c_q;

    //!
    //! @brief Graph BoundaryFaces-2-Nodes.
    //!
    wmesh_int_sparsemat_t 	m_bf2n;
    
    //!
    //! @brief Static Reference Graph Edges-2-Nodes.
    //!
    wmesh_int_sparsemat_t 	m_s_e2n;    

    //!
    //! @brief Static Reference Graph Triangles-2-Nodes.
    //!
    wmesh_int_sparsemat_t 	m_s_t2n;    

    //!
    //! @brief Static Reference Graph Quadrilaterals-2-Nodes.
    //!
    wmesh_int_sparsemat_t 	m_s_q2n;


    //
    // Cell codes.
    //
    wmesh_int_sparsemat_t	m_c_c;
    //
    // Boundary faces cdes.
    //
    wmesh_int_sparsemat_t	m_bf_c;
    
    //!
    //! @brief Coordinates.
    //!
    double * 			m_coo;    
    wmesh_int_t			m_coo_ld;
    wmesh_int_mat_t		m_n_c;
    
    wmesh_int_t 		m_num_entities[32];
    wmesh_int_t 		m_num_bfaces[2];
    wmesh_int_t 		m_num_cells;
    wmesh_int_t 		m_num_nodes;
    wmesh_int_t 		m_num_edges;
    wmesh_int_t 		m_num_faces;
    wmesh_int_t 		m_num_volumes;
    wmesh_int_t 		m_num_triangles;
    wmesh_int_t 		m_num_quadrilaterals;
    wmesh_int_t 		m_num_tetrahedra;
    wmesh_int_t 		m_num_pyramids;
    wmesh_int_t 		m_num_wedges;
    wmesh_int_t 		m_num_hexahedra;
    
#if 1
    wmesh_int_sparsemat_t 	m_c2d;
    wmesh_int_sparsemat_t 	m_c2d_n;
    wmesh_int_sparsemat_t 	m_c2d_e;
    wmesh_int_sparsemat_t 	m_c2d_t;
    wmesh_int_sparsemat_t 	m_c2d_q;
    wmesh_int_sparsemat_t 	m_c2d_i;
    double * 			m_coo_dofs;    
    wmesh_int_t 		m_ndofs;
#endif

  };

  struct wmeshspace_t
  {
    wmesh_t * 			m_mesh;
    wmesh_int_t 		m_degree;
    wmesh_t * 			m_patterns[4];

    wmesh_int_sparsemat_t 	m_c2d;
    wmesh_int_sparsemat_t 	m_c2d_n;
    wmesh_int_sparsemat_t 	m_c2d_e;
    wmesh_int_sparsemat_t 	m_c2d_t;
    wmesh_int_sparsemat_t 	m_c2d_q;
    wmesh_int_sparsemat_t 	m_c2d_i;
    wmesh_int_t 		m_ndofs;
    
    //    double* 		m_coo_dofs;
    //    wmesh_int_t		m_coo_dofs_ld;
  };

};
