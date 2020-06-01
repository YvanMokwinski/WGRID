#pragma once

#include "wmesh.h"
#include "wmesh-types.hpp"

#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif

  wmesh_status_t wmesh_init_c2c(const wmesh_int_t 		topodim_,
				const wmesh_int_sparsemat_t*	c2n_,
				wmesh_int_sparsemat_t*		c2c_);

  wmesh_status_t wmesh_init_bf2n(wmesh_int_t 			num_triangles_,
				 wmesh_int_t 			num_quadrilaterals_,
				 wmesh_int_sparsemat_t*		bf2n_);
  wmesh_status_t wmesh_init_c2f_q(wmesh_int_sparsemat_t*	c2f_,
				  wmesh_int_sparsemat_t*	c2f_q_);
  wmesh_status_t wmesh_init_c2f_t(wmesh_int_sparsemat_t*	c2f_,
				  wmesh_int_sparsemat_t*	c2f_t_);
  
  wmesh_status_t wmesh_init_c2e(const wmesh_int_sparsemat_t*	c2n_,
				wmesh_int_sparsemat_t*		c2e_,
				wmesh_int_t 			topology_dimension_);

     wmesh_status_t wmesh_init_c2f(const wmesh_int_sparsemat_t*	c2n_,
				   wmesh_int_sparsemat_t*		c2f_);

    wmesh_status_t wmesh_create_c2c(wmesh_t* 		self_,
				  wmesh_int_p 		num_faces_,
				    wmesh_int_p 		num_boundary_faces_);


  wmesh_status_t wmesh_build_s_e2n		(wmesh_int_sparsemat_t* 	self_,
						 wmesh_int_t 			dim_);
  
  wmesh_status_t wmesh_build_s_t2n		(wmesh_int_sparsemat_t* 	self_);
  
  wmesh_status_t wmesh_build_s_q2n		(wmesh_int_sparsemat_t* 	self_);

  struct wmesh_t
  {
    wmesh_int_t m_topology_dimension;
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
    wmesh_int_t			m_coo_m;
    wmesh_int_t			m_coo_n;
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


  inline double * wmesh_get_coo(wmesh_t*self_)
  {
    return self_->m_coo;
  }

  
static  inline void get_c2e(wmesh_int_t 	m_,
			    const_wmesh_int_p 	c2e_,
			    wmesh_int_t 	c2e_ld_,
			    wmesh_int_t 	cell_idx_,
			    wmesh_int_p 	cnc_)
{
  for (wmesh_int_t i=0;i<m_;++i)
    {
      cnc_[i] = c2e_[cell_idx_*c2e_ld_+i];      
    }
};

  static inline void get_c2n(wmesh_int_t 	numCellNodes_,
			     const_wmesh_int_p 	c2n_,
			     wmesh_int_t 	c2n_ld_,
			     wmesh_int_t 	idx_,
			     wmesh_int_p 	cnc_)
  {
    for (wmesh_int_t localNodeIndex=0;localNodeIndex<numCellNodes_;++localNodeIndex)
      {
	cnc_[localNodeIndex] = c2n_[idx_*c2n_ld_+localNodeIndex];      
      }
  }


static  inline void get_e2n(const_wmesh_int_p	c2n_,
			    const wmesh_int_t 	e_idx_,
			    wmesh_int_p		e2n_,
			    wmesh_int_t 	s_e2n_m_,
			    wmesh_int_t 	s_e2n_n_,
			    const_wmesh_int_p 	s_e2n_v_,
			    wmesh_int_t 	s_e2n_ld_)		    
{
  e2n_[0] = c2n_[s_e2n_v_[s_e2n_ld_ * e_idx_+ 0]];
  e2n_[1] = c2n_[s_e2n_v_[s_e2n_ld_ * e_idx_+ 1]];
}

static  inline void get_q2n(const_wmesh_int_p		c2n_,
			    const wmesh_int_t 		t_lidx_,
			    wmesh_int_p			q2n_,
			    wmesh_int_t 		s_q2n_m_,
			    wmesh_int_t 		s_q2n_n_,
			    const_wmesh_int_p 		s_q2n_v_,
			    wmesh_int_t 		s_q2n_ld_)		    
{
  q2n_[0] = c2n_[s_q2n_v_[s_q2n_ld_ * t_lidx_ + 0]];
  q2n_[1] = c2n_[s_q2n_v_[s_q2n_ld_ * t_lidx_ + 1]];
  q2n_[2] = c2n_[s_q2n_v_[s_q2n_ld_ * t_lidx_ + 2]];
  q2n_[3] = c2n_[s_q2n_v_[s_q2n_ld_ * t_lidx_ + 3]];
};

  static  inline void get_t2n(const_wmesh_int_p		c2n_,
			      const wmesh_int_t 	t_lidx_,
			      wmesh_int_p		t2n_,
			      wmesh_int_t 		s_t2n_m_,
			      wmesh_int_t 		s_t2n_n_,
			      const_wmesh_int_p 	s_t2n_v_,
			      wmesh_int_t 		s_t2n_ld_)		    
  {
    t2n_[0] = c2n_[s_t2n_v_[s_t2n_ld_ * t_lidx_ + 0]];
    t2n_[1] = c2n_[s_t2n_v_[s_t2n_ld_ * t_lidx_ + 1]];
    t2n_[2] = c2n_[s_t2n_v_[s_t2n_ld_ * t_lidx_ + 2]];
  };
  
  static inline void get_x2n(const_wmesh_int_p	c2n_,
			     const wmesh_int_t 	x_lidx_,
			     wmesh_int_p	x2n_,
			     wmesh_int_t 	s_x2n_m_,
			     wmesh_int_t 	s_x2n_n_,
			     const_wmesh_int_p 	s_x2n_v_,
			     wmesh_int_t 	s_x2n_ld_)		    
  {    
    for (wmesh_int_t i=0;i<s_x2n_m_;++i)
      {
	x2n_[i] = c2n_[s_x2n_v_[s_x2n_ld_ * x_lidx_ + i]];
      }
  };


  wmesh_int_t 		wmesh_read_medit	(wmesh_t**self_,
						 bool is_binary_,
						 const char * filename,
						 ...);
  
  wmesh_status_t 	wmesh_write_medit	(const wmesh_t* 	self_,
						 bool  			is_binary_,
						 const char * 		filename_,
						 ...);

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
    wmesh_int_p			m_dof_codes;

    
    //    double* 		m_coo_dofs;
    //    wmesh_int_t		m_coo_dofs_ld;
  };

#ifdef __cplusplus
}
#endif
