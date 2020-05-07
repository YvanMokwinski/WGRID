
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"
#include "bms.h"
#include "wmesh.hpp"
#include <chrono>
#include <iostream>

using namespace std::chrono;
extern "C"
{

  wmesh_status_t wmesh_shape_eval_basis(wmesh_int_t     cell_type_,
					wmesh_int_t 	coo_m_,
					wmesh_int_t 	coo_n_,
					double * 	coo_x_,
					wmesh_int_t 	coo_ld_,
					double * 	x_,
					wmesh_int_t 	x_ld_)
  {
    
    static const double s_zero = ((double)0.0);
    static const double s_one = ((double)1.0);
    static constexpr double one = 1.0;
    if (cell_type_==0)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];
	    
	    x_[k*x_ld_+0] = s_one - (r+s+t);
	    x_[k*x_ld_+1] = r;
	    x_[k*x_ld_+2] = s;
	    x_[k*x_ld_+3] = t;
	  }		  
      }
    else if (cell_type_==1)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];
	    
	    const auto
	      tscale = (t < s_one) ? s_one / (s_one - t) : s_zero,
	      lts = s_one - t - s,
	      ltr = s_one - t - r;
	    
	    x_[k*x_ld_+0] = ltr  * lts * tscale;
	    x_[k*x_ld_+1] = r * lts * tscale;
	    x_[k*x_ld_+2] = r * s * tscale;
	    x_[k*x_ld_+3] = ltr * s  * tscale;
	    x_[k*x_ld_+4] = t;
	  }		  
      }
    else if (cell_type_==2)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];
	   
	    const auto
	      l = one - (r+s),
	      lt = one - t;
	    
	    x_[k*x_ld_+0] = l * lt;
	    x_[k*x_ld_+1] = r * lt;
	    x_[k*x_ld_+2] = s * lt;
	    x_[k*x_ld_+3] = l * t;
	    x_[k*x_ld_+4] = r * t;
	    x_[k*x_ld_+5] = s * t;
	  }		  
      }
    else if (cell_type_==3)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];

	    const auto
	      lr = one - r,
	      ls = one - s,
	      lt = one - t;
	    
	    x_[k*x_ld_+0] = lr* ls * lt;
	    x_[k*x_ld_+1] = r * ls * lt;
	    x_[k*x_ld_+2] = r *  s * lt;
	    x_[k*x_ld_+3] = lr*  s * lt;
	    x_[k*x_ld_+4] = lr* ls * t;
	    x_[k*x_ld_+5] =  r* ls * t;
	    x_[k*x_ld_+6] =  r*  s * t;
	    x_[k*x_ld_+7] = lr*  s * t;
	  }
      }
    else
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
      }
    return WMESH_STATUS_SUCCESS;
  }

  




  
  
  
  wmesh_status_t wmesh_analysis_edges(wmesh_t* 		self_,
				      wmesh_int_p 	num_edges_)
  {
    wmesh_int_t edge_idx = 0;
    
    wmesh_int_t work_n;
    wmesh_status_t status =  wmesh_indexing_edges_buffer_size(4,
							      self_->m_c2n.m_n,
							      &work_n);      
    WMESH_STATUS_CHECK(status);
    
    wmesh_int_t * work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
    if (!work)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }      
    status =  wmesh_indexing_edges(WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2n),
				   WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2e),
				   WMESH_INT_SPARSEMAT_FORWARD(self_->m_s_e2n),				     
				   &edge_idx,
				   work_n,
				   work);
    free(work);
    WMESH_STATUS_CHECK(status);
    num_edges_[0] = edge_idx;
    return WMESH_STATUS_SUCCESS;
  }
  

  wmesh_status_t wmesh_analysis_faces(wmesh_t* 		self_)
  {
    wmesh_status_t status;
    
    wmesh_int_t work_n;
    status =  wmesh_indexing_triangles_buffer_size(self_->m_c2n.m_size,
						   self_->m_c2n.m_n,
						   &work_n);
    WMESH_STATUS_CHECK(status);
    wmesh_int_p work = (wmesh_int_p)malloc(work_n*sizeof(wmesh_int_t));
    if (!work)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }      
    
    //
    // Enumerates triangles faces.
    //
    wmesh_int_t face_idx = 0;
    status =  wmesh_indexing_triangles(WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2n),
				       WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2f_t),
				       WMESH_INT_SPARSEMAT_FORWARD(self_->m_s_t2n),				     
				       &face_idx,
				       work_n,
				       work);      
    
    WMESH_STATUS_CHECK(status);      
    face_idx = 0;
    status = wmesh_indexing_quadrilaterals(WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2n),
					   WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2f_q),
					   WMESH_INT_SPARSEMAT_FORWARD(self_->m_s_q2n),				     
					   &face_idx,
					   work_n,
					   work);      
    WMESH_STATUS_CHECK(status);
    return WMESH_STATUS_SUCCESS;
  }

  
  wmesh_status_t wmesh_analysis_neighbors(wmesh_t* 		self_,
					  wmesh_int_p 		num_faces_,
					  wmesh_int_p 		num_boundary_faces_)
  {    
        
    wmesh_int_t work_n;
    wmesh_int_t * work  = nullptr;
    wmesh_status_t status;
    
    status =  wbms_c2c_calculate_buffer_size(self_->m_c2n.m_size,
					     self_->m_c2n.m_n,
					     &work_n);    
    WMESH_STATUS_CHECK(status);    
    work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
    if (!work)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }      

    num_boundary_faces_[0] 	= 0;
    num_boundary_faces_[1] 	= 0;
    num_faces_[0] 		= 0;
    num_faces_[1] 		= 0;
    
    status = wbms_c2c_calculate(WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2n),
				WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2c_t),
				WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2c_q),
				WMESH_INT_SPARSEMAT_FORWARD(self_->m_s_t2n),
				WMESH_INT_SPARSEMAT_FORWARD(self_->m_s_q2n),
				work_n,
				work,
				num_faces_,
				num_boundary_faces_);

    free(work);      
    WMESH_STATUS_CHECK(status);    
    return WMESH_STATUS_SUCCESS;
    
#undef WINT_SPARSE_MAT_PARAMS2    
  }



  
  using timing_t = high_resolution_clock::time_point;
  inline timing_t timing_stop(){ return high_resolution_clock::now(); }
  inline timing_t timing_start(){ return high_resolution_clock::now(); }
  inline double timing_seconds(timing_t&t,timing_t&t1){ return duration_cast<duration<double>>(t1-t).count(); }




  
  wmesh_status_t wmesh_analysis(wmesh_t*    self_,
				wmesh_int_t degree)
  {
    const wmesh_int_t topology_dimension = self_->m_topology_dimension;
    
    wmesh_status_t 	status;
    wmesh_int_t 	num_faces[2];
    wmesh_int_t 	num_bfaces[2];
    
#ifndef NDEBUG
    std::cerr << "wmesh_analysis: topology_dimension " << topology_dimension << std::endl;
#endif

    
#if 0
    status    = wmesh_analysis_neighbors(self_,
					 num_faces,
					 num_bfaces);
    
    self_->m_num_triangles      = num_faces[0];
    self_->m_num_quadrilaterals = num_faces[1];
    self_->m_num_bfaces[0]      = num_bfaces[0];
    self_->m_num_bfaces[1]      = num_bfaces[1];
    
    //
    //
    //
    status = wmesh_adjacencies(self_);

#endif

    //
    // 1/ Create the cells to cells.
    //
    //
    // 1.a/ Use c2f since this is the same layout.
    //
    status = wmesh_init_c2f	(&self_->m_c2n,
				 &self_->m_c2c);
    WMESH_STATUS_CHECK(status);
	
    if (topology_dimension == 3)
      {
	status = wmesh_init_c2f_t(&self_->m_c2c,
				  &self_->m_c2c_t);
	WMESH_STATUS_CHECK(status);
	
	//
	// 1.a/ Define the part through quadrilaterals.
	//
	status = wmesh_init_c2f_q	(&self_->m_c2c,
					 &self_->m_c2c_q);
	WMESH_STATUS_CHECK(status);      
      }
  


    //
    // 3/ Create the cells to edges.
    //    
    {
      wmesh_status_t status;
      status = wmesh_init_c2e(&self_->m_c2n,
			      &self_->m_c2e,
			      self_->m_topology_dimension);
      WMESH_STATUS_CHECK(status);
    }
    
    //
    // 4/ Create the cells to faces.
    //    
    {
      wmesh_status_t status;
      status = wmesh_init_c2f(&self_->m_c2n,
			      &self_->m_c2f);
      WMESH_STATUS_CHECK(status);
      
      //
      // 4.a/ Define the part of triangles.
      //      
      status = wmesh_init_c2f_t	(&self_->m_c2f,
				 &self_->m_c2f_t);
      WMESH_STATUS_CHECK(status);
      
      //
      // 4.b/ Define the part of quadrilaterals.
      //      
      status = wmesh_init_c2f_q	(&self_->m_c2f,
				 &self_->m_c2f_q);
      WMESH_STATUS_CHECK(status);      
    }
    
    //
    // Initialize the reference shape edges to nodes.
    //
    {
      wmesh_status_t status;
#if 1
      {
	auto start = timing_start();      
	status    = wmesh_analysis_neighbors(self_,
					     num_faces,
					     num_bfaces);
	
	self_->m_num_triangles      = num_faces[0];
	self_->m_num_quadrilaterals = num_faces[1];
	self_->m_num_bfaces[0]      = num_bfaces[0];
	self_->m_num_bfaces[1]      = num_bfaces[1];
	auto stop = timing_stop();
	std::cout << "neighbors detection, elapsed " << timing_seconds(start,stop) << std::endl;
      }
#endif      
      //      WMESH_STATUS_CHECK(status);

      {
	auto start = timing_start();      
	status = wmesh_analysis_edges(self_, &self_->m_num_edges);
	auto stop = timing_stop();
	std::cout << "edge analysis, elapsed " << timing_seconds(start,stop) << std::endl;
	WMESH_STATUS_CHECK(status);
      }
      
      {
	std::cout << "// wmesh_analysis ... faces indexing ..."  << std::endl;
	auto start = timing_start();      
	status    = wmesh_analysis_faces(self_);
	auto stop = timing_stop();
	std::cout << "faces analysis, elapsed " << timing_seconds(start,stop) << std::endl;
	WMESH_STATUS_CHECK(status);
      }

      WMESH_STATUS_CHECK(status);
    }

    std::cout << "boundary faces (triangles, quadrilaterals):  "
	      << num_bfaces[0]
	      << ", "
	      << num_bfaces[1]
	      << std::endl;
    
    std::cout << "total num faces (triangles, quadrilaterals):  "
	      << num_faces[0]
	      << ", "
	      << num_faces[1] << std::endl;
    
    std::cout << "total num edges:  " << self_->m_num_edges  << std::endl;




#if 0
    
    wmesh_int_t ndofs[4]= { ((degree+1)*(degree+2)*(degree+3)) / 6,
			    ((degree+2)*(degree+1)*(2*degree+3)) / 6,
			    ((degree+1)*(degree+1)*(degree+2)) / 2,
			    (degree+1)*(degree+1)*(degree+1)};
      
      
    wmesh_int_t ndofs_n[4] = {1,1,1,1};
    wmesh_int_t ndofs_e[4] = {degree-1,degree-1,degree-1,degree-1};

    wmesh_int_t ndofs_t[4] = {(degree>2) ? ((degree-1)*(degree-2))/2 : 0,
			      (degree>2) ? ((degree-1)*(degree-2))/2 : 0,
			      (degree>2) ? ((degree-1)*(degree-2))/2 : 0,
			      (degree>2) ? ((degree-1)*(degree-2))/2 : 0};

    wmesh_int_t ndofs_q[4] = {(degree>1) ? (degree-1)*(degree-1) : 0,
			      (degree>1) ? (degree-1)*(degree-1) : 0,
			      (degree>1) ? (degree-1)*(degree-1) : 0,
			      (degree>1) ? (degree-1)*(degree-1) : 0};
            
    wmesh_int_t ndofs_i[4] = { (degree>0) ? ((degree-1)*(degree-2)*(degree-3)) / 6 : 1,
			       (degree>0) ? ( (degree-2)*(degree-1)*(2*degree-3) ) / 6 : 1,
			       (degree>0) ? ((degree-1)*(degree-1)*(degree-2)) / 2 : 1, //				 (degree>2) ? (degree-2)*(degree-1) : 0,
			       (degree>0) ? (degree-1)*(degree-1)*(degree-1) : 1 };
#if 1
    //
    // cells to dofs.
    //
    {
      wmesh_status_t status;
      
      //
      // Use c2f since this is the same layout.
      //
      wmesh_int_t shifts[4] = {0,0,0,0};
      
      status = wmesh_init_c2d	(&self_->m_c2n,
				 &self_->m_c2d,
				 ndofs,
				 ndofs);
      WMESH_STATUS_CHECK(status);

      //
      // Node-based dofs.
      //
      status = wmesh_init_c2d_n	(&self_->m_c2d,
				 &self_->m_c2d_n,
				 ndofs_n,
				 shifts);
      WMESH_STATUS_CHECK(status);

      //
      // edge-based dofs.
      //
      status = wmesh_init_c2d_e	(&self_->m_c2d,
				 &self_->m_c2d_e,
				 ndofs_e,
				 shifts);
      WMESH_STATUS_CHECK(status);

      //
      // triangle-based dofs.
      //
      status = wmesh_init_c2d_t	(&self_->m_c2d,
				 &self_->m_c2d_t,
				 ndofs_t,
				 shifts);
      WMESH_STATUS_CHECK(status);
      
      //
      // quadrilateral-based dofs.
      //
      status = wmesh_init_c2d_q	(&self_->m_c2d,
				 &self_->m_c2d_q,
				 ndofs_q,
				 shifts);
      WMESH_STATUS_CHECK(status);      
      
      //
      // cell-based dofs.
      //
      status = wmesh_init_c2d_i(&self_->m_c2d,
				&self_->m_c2d_i,
				ndofs_i,
				shifts);
      WMESH_STATUS_CHECK(status);
    }
#endif
#if 0
    {
      auto start = timing_start();      
      wmesh_status_t status = wmeshspace_analysis(self_,
						  degree);
      auto stop = timing_stop();
      std::cout << "analysis space, elapsed " << timing_seconds(start,stop) << std::endl;
      WMESH_STATUS_CHECK(status);
    }
#endif
    {
#if 0
      
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_n,stdout);
      fprintf(stdout,"------\n");
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_e,stdout);
      fprintf(stdout,"------\n");
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_t,stdout);
      fprintf(stdout,"------\n");
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_q,stdout);
      fprintf(stdout,"------\n");
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_i,stdout);
#endif
    }


#if 0
    {
      for (wmesh_int_t degree = 2;degree<=11;++degree)
	{
	  
	  wmesh_int_t ndofs[4]= { ((degree+1)*(degree+2)*(degree+3)) / 6,
				  ((degree+2)*(degree+1)*(2*degree+3)) / 6,
				  ((degree+1)*(degree+1)*(degree+2)) / 2,
				  (degree+1)*(degree+1)*(degree+1)};
	  wmesh_int_t ndofs_per_node = 1;
	  wmesh_int_t ndofs_per_edge = degree-1;
	  wmesh_int_t ndofs_per_faces[2] = {(degree>2) ? ((degree-1)*(degree-2))/2 : 0,
					    (degree>1) ? (degree-1)*(degree-1) : 0};
	  
	  wmesh_int_t ndofs_per_vol[4] = { (degree>2) ? ((degree+1)*(degree+2)*(degree+3)) / 6 - ndofs_per_faces[0]*4-4-6*ndofs_per_edge : 0,
					   (degree>0) ? ( (degree-2)*(degree-1)*(2*degree-3) ) / 6 : 0,
					   (degree>2) ? (degree-2)*(degree-1) : 0,
					   (degree>1) ? (degree-1)*(degree-1)*(degree-1) : 0 };
#if 0
	  std::cout << "edge part "<<ndofs_per_edge * num_edges[0] << std::endl;
	  std::cout << "face part "<<num_faces[0] * ndofs_per_faces[0] +
	    num_faces[1] * ndofs_per_faces[1]  << std::endl;
	  std::cout << "volume part "<<self_->m_c2n.m_n[0] * ndofs_per_vol[0] +
	    self_->m_c2n.m_n[1] * ndofs_per_vol[1] +
	    self_->m_c2n.m_n[2] * ndofs_per_vol[2] +
	    self_->m_c2n.m_n[3] * ndofs_per_vol[3] << std::endl;
#endif
	  std::cout << self_->m_num_nodes << std::endl;
	  std::cout << "P1 " <<  ndofs_per_node * self_->m_num_nodes << std::endl;
	  std::cout << "P2 " <<  ndofs_per_node * self_->m_num_nodes + ndofs_per_edge * self_->m_num_edges << std::endl;
	  std::cout << "P3 " <<
	    num_faces[0] * ndofs_per_faces[0] +
	    num_faces[1] * ndofs_per_faces[1] +
	    ndofs_per_node * self_->m_num_nodes +
	    ndofs_per_edge * self_->m_num_edges << std::endl;

	  for (int i=0;i<4;++i) std::cout << "$#$# " << ndofs_per_vol[i]<< std::endl;
	  unsigned long long int total_ndofs =
	    ndofs_per_node * self_->m_num_nodes +
	    ndofs_per_edge * self_->m_num_edges +
	    num_faces[0] * ndofs_per_faces[0] +
	    num_faces[1] * ndofs_per_faces[1] +
	    self_->m_c2n.m_n[0] * ndofs_per_vol[0] +
	    self_->m_c2n.m_n[1] * ndofs_per_vol[1] +
	    self_->m_c2n.m_n[2] * ndofs_per_vol[2] +
	    self_->m_c2n.m_n[3] * ndofs_per_vol[3];

	  unsigned long long int total_dndofs =
	    self_->m_c2n.m_n[0] * ndofs[0] +
	    self_->m_c2n.m_n[1] * ndofs[1] +
	    self_->m_c2n.m_n[2] * ndofs[2] +
	    self_->m_c2n.m_n[3] * ndofs[3];
	  
	  std::cout << "total ndofs[" << degree << "]:  " << total_ndofs  << " " << total_dndofs  << std::endl;
	}
      
    }
#endif
#endif
    return WMESH_STATUS_SUCCESS;
  }
  


};
