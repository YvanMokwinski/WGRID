#ifndef BMS_H
#define BMS_H
#include <stdlib.h>
#include "wmesh-status.h"

#ifdef __cplusplus
extern "C"
{
#endif




  wmesh_status_t bms_n2c	(wmesh_int_t 		c2n_size_,
				 const_wmesh_int_p 	c2n_ptr_,
				 const_wmesh_int_p 	c2n_m_,
				 const_wmesh_int_p 	c2n_n_,
				 const_wmesh_int_p 	c2n_v_,
				 const_wmesh_int_p 	c2n_ld_,				 
				 wmesh_int_p 		n2c_ptr_,
				 wmesh_int_t 		n2c_m_,
				 wmesh_int_p 		n2c_v_);
  
  wmesh_status_t bms_n2c_cindex	(wmesh_int_t 		c_,
				 wmesh_int_p 		cindex_);
  
  wmesh_status_t bms_n2c_ctype	(wmesh_int_t 		c_,
				 wmesh_int_p 		ctype_);
  
  //!
  //! @brief Get the size of the needed buffer.
  //!
  wmesh_status_t
  wbms_c2c_t_calculate_buffer_size
  (wmesh_int_t		c2n_size_,
   const_wmesh_int_p 	c2n_n_,
   wmesh_int_p 		work_n_);
    
  //!
  //! @brief Build the cell-to-cell graph through triangle faces.
  //!
  wmesh_status_t
  wbms_c2c_t_calculate
  (wmesh_int_t		c2n_size_,   
   const_wmesh_int_p 	c2n_ptr_,
   const_wmesh_int_p 	c2n_m_,
   const_wmesh_int_p 	c2n_n_,
   const_wmesh_int_p 	c2n_v_,
   const_wmesh_int_p 	c2n_ld_,
		       
   const_wmesh_int_p 	c2c_t_ptr_,
   const_wmesh_int_p 	c2c_t_m_,
   const_wmesh_int_p 	c2c_t_n_,
   wmesh_int_p 		c2c_t_v_,
   const_wmesh_int_p 	c2c_t_ld_,
		       
   const_wmesh_int_p 	s_t2n_ptr_,
   const_wmesh_int_p	s_t2n_m_,
   const_wmesh_int_p	s_t2n_n_,
   const_wmesh_int_p 	s_t2n_v_,
   const_wmesh_int_p	s_t2n_ld_,
		       
   wmesh_int_p		num_boundary_triangles_,
   wmesh_int_t		work_n_,
   wmesh_int_p		work_);
  
  //!
  //! @brief Get the size of the needed buffer.
  //!
  wmesh_status_t
  wbms_c2c_q_calculate_buffer_size
  (wmesh_int_t		c2n_size_,
   const_wmesh_int_p 	c2n_n_,
   wmesh_int_p 		work_n_);
  
  //!
  //! @brief Build the cell-to-cell graph through quadrilateral faces.
  //!
  wmesh_status_t  wbms_c2c_q_calculate
  (wmesh_int_t		c2n_size_,
   const_wmesh_int_p 	c2n_ptr_,
   const_wmesh_int_p 	c2n_m_,
   const_wmesh_int_p 	c2n_n_,
   const_wmesh_int_p 	c2n_v_,
   const_wmesh_int_p 	c2n_ld_,
   
   const_wmesh_int_p 	c2c_q_ptr_,
   const_wmesh_int_p 	c2c_q_m_,
   const_wmesh_int_p 	c2c_q_n_,
   wmesh_int_p 		c2c_q_v_,
   const_wmesh_int_p 	c2c_q_ld_,
   
   const_wmesh_int_p 	s_q2n_ptr_,
   const_wmesh_int_p	s_q2n_m_,
   const_wmesh_int_p	s_q2n_n_,
   const_wmesh_int_p 	s_q2n_v_,
   const_wmesh_int_p	s_q2n_ld_,
   
   wmesh_int_p		num_boundary_quadrilaterals_,
   wmesh_int_t		work_n_,
   wmesh_int_p		work_);
  
  //!
  //! @brief Build the cell-to-cell graph.
  //!
  wmesh_status_t
  wbms_c2c_calculate_buffer_size
  (wmesh_int_t		c2n_size_,
   const_wmesh_int_p 	c2n_n_,
   wmesh_int_p 		work_n_);
  
  //!
  //! @brief Build the cell-to-cell graph.
  //!
  wmesh_status_t
  wbms_c2c_calculate
  (wmesh_int_t		c2n_size_,
   const_wmesh_int_p 	c2n_ptr_,
   const_wmesh_int_p 	c2n_m_,
   const_wmesh_int_p 	c2n_n_,
   const_wmesh_int_p 	c2n_v_,
   const_wmesh_int_p 	c2n_ld_,
   
   wmesh_int_t		c2c_t_size_,
   const_wmesh_int_p 	c2c_t_ptr_,
   const_wmesh_int_p 	c2c_t_m_,
   const_wmesh_int_p 	c2c_t_n_,
   wmesh_int_p		c2c_t_v_,
   const_wmesh_int_p 	c2c_t_ld_,
						 
   wmesh_int_t		c2c_q_size_,
   const_wmesh_int_p 	c2c_q_ptr_,
   const_wmesh_int_p 	c2c_q_m_,
   const_wmesh_int_p 	c2c_q_n_,
   wmesh_int_p		c2c_q_v_,
   const_wmesh_int_p 	c2c_q_ld_,
						 
   wmesh_int_t 		s_t2n_size_,
   const_wmesh_int_p 	s_t2n_ptr_,
   const_wmesh_int_p 	s_t2n_m_,
   const_wmesh_int_p 	s_t2n_n_,
   const_wmesh_int_p 	s_t2n_v_,
   const_wmesh_int_p 	s_t2n_ld_,
						 
   wmesh_int_t 		s_q2n_size_,
   const_wmesh_int_p 	s_q2n_ptr_,
   const_wmesh_int_p 	s_q2n_m_,
   const_wmesh_int_p 	s_q2n_n_,
   const_wmesh_int_p 	s_q2n_v_,
   const_wmesh_int_p 	s_q2n_ld_,
						 
   wmesh_int_t 		work_n_,
   wmesh_int_p 		work_,
   wmesh_int_p 		num_faces_,
   wmesh_int_p 		num_bfaces_);


  
  wmesh_status_t bms_read_medit_stat(int64_t inm_,
				     wmesh_int_p num_nodes_,
				     wmesh_int_p num_edges_,
				     wmesh_int_p num_triangles_,
				     wmesh_int_p num_quadrilaterals_,
				     wmesh_int_p num_tetrahedra_,
				     wmesh_int_p num_pyramids_,
				     wmesh_int_p num_wedges_,
				     wmesh_int_p num_hexahedra_);

  wmesh_status_t bms_read_medit_open(wmesh_str_t 	filename_,
				     int64_t*		inm_,
				     int32_t * 		version_,
				     int32_t * 		dim_);
  
  wmesh_status_t bms_write_medit_open(wmesh_str_t 	filename_,
				      int64_t*		inm_);  

  wmesh_status_t bms_medit_close(int64_t		inm_);
  
  wmesh_status_t
  bms_write_medit_topology
  (int64_t				inm_,
   
   const wmesh_int_t 			c2n_size,
   const_wmesh_int_p 			c2n_ptr,
   const_wmesh_int_p 			c2n_m,
   const_wmesh_int_p 			c2n_n,
   const_wmesh_int_p 			c2n_v,
   const_wmesh_int_p 			c2n_ld,
   
   const wmesh_int_t 			c_c_size,
   const_wmesh_int_p 			c_c_ptr,
   const_wmesh_int_p 			c_c_m,
   const_wmesh_int_p 			c_c_n,
   const_wmesh_int_p 			c_c_v,
   const_wmesh_int_p 			c_c_ld);
  
  
  wmesh_status_t
  bms_write_medit_geometry
  (int64_t 				inm_,
   wmesh_int_t				coo_m_,
   wmesh_int_t				coo_n_,
   const double *__restrict__ 		coo_,
   wmesh_int_t 				coo_ld_,
   const_wmesh_int_p  			nflags_,
   wmesh_int_t 				nflags_ld_);


  wmesh_status_t
  bms_read_medit_topology
  (int64_t				inm_,
   
   const wmesh_int_t 			c2n_size,
   const_wmesh_int_p 			c2n_ptr,
   const_wmesh_int_p 			c2n_m,
   const_wmesh_int_p 			c2n_n,
   wmesh_int_p 				c2n_v,
   const_wmesh_int_p 			c2n_ld,
			  
   const wmesh_int_t 			c_c_size,
   const_wmesh_int_p 			c_c_ptr,
   const_wmesh_int_p 			c_c_m,
   const_wmesh_int_p 			c_c_n,
   wmesh_int_p 				c_c_v,
   const_wmesh_int_p 			c_c_ld);

  wmesh_status_t
  bms_read_medit_geometry
  (int64_t 				inm_,
   wmesh_int_t				coo_m_,
   wmesh_int_t				coo_n_,
   double *__restrict__ 		coo_,
   wmesh_int_t 				coo_ld_,
   wmesh_int_p  			nflags_,
   wmesh_int_t 				nflags_ld_);


  
#ifdef __cplusplus
}
#endif


#endif
