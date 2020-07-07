#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "bms.h"
#include "wmesh_t.hpp"
#include <chrono>
#include <iostream>

using namespace std::chrono;
extern "C"
{

  wmesh_status_t wmesh_c2e_calculate(const wmesh_int_sparsemat_t* 	c2n_,
				     wmesh_int_sparsemat_t* 		c2e_,
				     const wmesh_int_sparsemat_t* 	s_e2n_,
				     wmesh_int_p 			num_edges_)
  {
    WMESH_CHECK_POINTER(c2n_);
    WMESH_CHECK_POINTER(c2e_);
    WMESH_CHECK_POINTER(s_e2n_);
    WMESH_CHECK_POINTER(num_edges_);    
    wmesh_int_t edge_idx = 0;    
    wmesh_int_t work_n;
    wmesh_status_t status =  bms_c2e_buffer_size(c2n_->m_size,
						 c2n_->m_n,
						 &work_n);      
    WMESH_STATUS_CHECK(status);
    wmesh_int_t * work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
    if (!work)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      } 
    status =  bms_c2e(WMESH_INT_SPARSEMAT_FORWARD(*c2n_),
		      WMESH_INT_SPARSEMAT_FORWARD(*c2e_),
		      WMESH_INT_SPARSEMAT_FORWARD(*s_e2n_),				     
		      &edge_idx,
		      work_n,
		      work);
    free(work);
    WMESH_STATUS_CHECK(status);
    num_edges_[0] = edge_idx;
    return WMESH_STATUS_SUCCESS;
  }

  

  wmesh_status_t wmesh_c2f_calculate(const wmesh_int_sparsemat_t* 	c2n_,
				     const wmesh_int_sparsemat_t* 	s_t2n_,
				     wmesh_int_sparsemat_t* 		c2f_t_,
				     const wmesh_int_sparsemat_t* 	s_q2n_,
				     wmesh_int_sparsemat_t* 		c2f_q_)
  {
    WMESH_CHECK_POINTER(c2n_);
    WMESH_CHECK_POINTER(s_t2n_);
    WMESH_CHECK_POINTER(c2f_t_);
    WMESH_CHECK_POINTER(s_q2n_);
    WMESH_CHECK_POINTER(c2f_q_);

    wmesh_status_t status;    
    wmesh_int_t work_n0;
    wmesh_int_t work_n;

    status =  bms_c2f_t_buffer_size(c2n_->m_size,
				    c2n_->m_n,
				    &work_n0);
    WMESH_STATUS_CHECK(status);    
    status =  bms_c2f_q_buffer_size(c2n_->m_size,
				    c2n_->m_n,
				    &work_n);
    WMESH_STATUS_CHECK(status);
    work_n = (work_n0 > work_n) ? work_n0 : work_n;
    wmesh_int_p work 	= (wmesh_int_p)malloc(work_n*sizeof(wmesh_int_t));
    if (!work)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }      
    
    //
    // Enumerates triangles faces.
    //
    wmesh_int_t face_idx = 0;
    status =  bms_c2f_t(WMESH_INT_SPARSEMAT_FORWARD(*c2n_),
			WMESH_INT_SPARSEMAT_FORWARD(*c2f_t_),
			WMESH_INT_SPARSEMAT_FORWARD(*s_t2n_),				     
			&face_idx,
			work_n,
			work);      
    
    WMESH_STATUS_CHECK(status);      
    face_idx = 0;
    status = bms_c2f_q(WMESH_INT_SPARSEMAT_FORWARD(*c2n_),
		       WMESH_INT_SPARSEMAT_FORWARD(*c2f_q_),
		       WMESH_INT_SPARSEMAT_FORWARD(*s_q2n_),				     
		       &face_idx,
		       work_n,
		       work);      
    WMESH_STATUS_CHECK(status);
    return WMESH_STATUS_SUCCESS;
  }



  wmesh_status_t wmesh_c2c_calculate(const wmesh_int_t 			topodim_,
				     const wmesh_int_sparsemat_t* 	c2n_,
				     const wmesh_int_sparsemat_t* 	s_e2n_,
				     wmesh_int_sparsemat_t* 		c2c_e_,
				     const wmesh_int_sparsemat_t* 	s_t2n_,
				     wmesh_int_sparsemat_t* 		c2c_t_,
				     const wmesh_int_sparsemat_t* 	s_q2n_,
				     wmesh_int_sparsemat_t* 		c2c_q_,
				     wmesh_int_p 			num_hyperfaces_,
				     wmesh_int_p 			num_boundary_hyperfaces_)
  {
    wmesh_int_t status;
    wmesh_int_t work_n;
    wmesh_int_t * work  = nullptr;
    if (topodim_==3)
      {	
	status =  bms_c2c_buffer_size(c2n_->m_size,
				      c2n_->m_n,
				      &work_n);    
	WMESH_STATUS_CHECK(status);
	
	work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
	if (!work)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	  }      
	
	num_boundary_hyperfaces_[0] 	= 0;
	num_boundary_hyperfaces_[1] 	= 0;
	num_hyperfaces_[0] 		= 0;
	num_hyperfaces_[1] 		= 0;
	
	status = bms_c2c(WMESH_INT_SPARSEMAT_FORWARD(*c2n_),
			 WMESH_INT_SPARSEMAT_FORWARD(*c2c_t_),
			 WMESH_INT_SPARSEMAT_FORWARD(*c2c_q_),
			 WMESH_INT_SPARSEMAT_FORWARD(*s_t2n_),
			 WMESH_INT_SPARSEMAT_FORWARD(*s_q2n_),
			 work_n,
			 work,
			 num_hyperfaces_,
			 num_boundary_hyperfaces_);
	
	free(work);      
	WMESH_STATUS_CHECK(status);    
	return WMESH_STATUS_SUCCESS;
      }
    else if (topodim_==2)      
      {
	status =  bms_c2c_e_buffer_size(c2n_->m_size,
					c2n_->m_n,
					&work_n);    
	WMESH_STATUS_CHECK(status);
	
	work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
	if (!work)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	  }      
	
	num_boundary_hyperfaces_[0] = 0;
		
	status = bms_c2c_e(WMESH_INT_SPARSEMAT_FORWARD(*c2n_),
			   WMESH_INT_SPARSEMAT_FORWARD(*c2c_e_),
			   WMESH_INT_SPARSEMAT_FORWARD(*s_e2n_),
			   num_boundary_hyperfaces_,
			   work_n,
			   work);
	
	num_hyperfaces_[0] = (num_boundary_hyperfaces_[0] + c2n_->m_n[0]*3+c2n_->m_n[1]*4);
	WMESH_CHECK(num_hyperfaces_[0] % 2 == 0);
	num_hyperfaces_[0] /= 2;
	free(work);      
	WMESH_STATUS_CHECK(status);    
	return WMESH_STATUS_SUCCESS;

      }
    else
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_NOT_IMPLEMENTED);    
      }
  }


  wmesh_status_t wmesh_c2c(const wmesh_t* 		self_,
			   wmesh_int_sparsemat_t* 	c2c_,
			   wmesh_int_sparsemat_t* 	c2c_t_,
			   wmesh_int_sparsemat_t* 	c2c_q_,
			   wmesh_int_p 			num_facets_,
			   wmesh_int_p 			num_bfacets_)
  {
    const wmesh_int_t topodim = self_->m_topology_dimension;
    wmesh_int_sparsemat_t* c2e = nullptr;

    //
    // Initialize c2c.
    //
    wmesh_status_t status = wmesh_init_c2c	(topodim,
						 &self_->m_c2n,
						 c2c_);
    WMESH_STATUS_CHECK(status);	
    switch(topodim)
      {
      case 3:
	{
	  //
	  // Initialize sub c2c for triangles.
	  //
	  status = wmesh_init_c2f_t(c2c_,
				    c2c_t_);
	  WMESH_STATUS_CHECK(status);
	  
	  //
	  // Initialize sub c2c for quadrilaterals.
	  //
	  status = wmesh_init_c2f_q(c2c_,
				    c2c_q_);
	  WMESH_STATUS_CHECK(status);
	  break;	  
	}

      case 2:
	{
	  c2e = c2c_;
	  break;
	}
	
      }
    
    return  wmesh_c2c_calculate(topodim,
				&self_->m_c2n,
				&self_->m_s_e2n,
				c2e,
				&self_->m_s_t2n,
				c2c_t_,
				&self_->m_s_q2n,
				c2c_q_,
				num_facets_,
				num_bfacets_);
  }

  
  wmesh_status_t wmesh_c2e(const wmesh_t* 		self_,
			   wmesh_int_sparsemat_t* 	c2e_,				     
			   wmesh_int_p 			num_edges_)
  {
    return wmesh_c2e_calculate(&self_->m_c2n,
			       c2e_,
			       &self_->m_s_e2n,
			       num_edges_);
  }

  wmesh_status_t wmesh_c2f(const wmesh_t* 		self_,
			   wmesh_int_sparsemat_t* 	c2f_t_,
			   wmesh_int_sparsemat_t* 	c2f_q_)
  {
    return  wmesh_c2f_calculate(&self_->m_c2n,
				&self_->m_s_t2n,
				c2f_t_,
				&self_->m_s_q2n,
				c2f_q_);
  }

  using timing_t = high_resolution_clock::time_point;
  inline timing_t timing_stop(){ return high_resolution_clock::now(); }
  inline timing_t timing_start(){ return high_resolution_clock::now(); }
  inline double timing_seconds(timing_t&t,timing_t&t1){ return duration_cast<duration<double>>(t1-t).count(); }
  
  wmesh_status_t wmesh_analysis(wmesh_t*    self_)
  {
    const wmesh_int_t topology_dimension = self_->m_topology_dimension;
    wmesh_status_t 	status;
    wmesh_int_t 	num_facets[2];
    wmesh_int_t 	num_bfacets[2];
    
#ifndef NDEBUG
    std::cerr << "wmesh_analysis: topology_dimension " << topology_dimension << std::endl;
#endif

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
    if (topology_dimension == 3)
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


    

      auto start = timing_start();      
      status    = wmesh_c2c(self_,
			    &self_->m_c2c,
			    &self_->m_c2c_t,
			    &self_->m_c2c_q,
			    num_facets,
			    num_bfacets);

	if (self_->m_topology_dimension==3)
	  {
	    self_->m_num_triangles      = num_facets[0];
	    self_->m_num_quadrilaterals = num_facets[1];
	    self_->m_num_bfaces[0]      = num_bfacets[0];
	    self_->m_num_bfaces[1]      = num_bfacets[1];
	  }
	else if (self_->m_topology_dimension==2)
	  {
	    self_->m_num_edges      = num_facets[0];
	    self_->m_num_bfaces[0]  = num_bfacets[0];
	  }
	else if (self_->m_topology_dimension==1)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_NOT_IMPLEMENTED);
	  }
	
	auto stop = timing_stop();
	std::cout << "neighbors detection, elapsed " << timing_seconds(start,stop) << std::endl;

	//      WMESH_STATUS_CHECK(status);
	
	if (self_->m_topology_dimension==3)
	  {
	    start = timing_start();      
	    status = wmesh_c2e(self_, &self_->m_c2e, &self_->m_num_edges);
	    stop = timing_stop();
	    std::cout << "edge analysis, elapsed " << timing_seconds(start,stop) << std::endl;
	    WMESH_STATUS_CHECK(status);

	    start = timing_start();      
	    status    = wmesh_c2f(self_,
				  &self_->m_c2f_t,
				  &self_->m_c2f_q);
	    stop = timing_stop();
	    std::cout << "faces analysis, elapsed " << timing_seconds(start,stop) << std::endl;
	    WMESH_STATUS_CHECK(status);
	  }	
	else if (self_->m_topology_dimension==2)
	  {
	    start = timing_start();      
	    status = wmesh_c2e(self_, &self_->m_c2e, &self_->m_num_edges);
	    stop = timing_stop();
	    std::cout << "edge analysis, elapsed " << timing_seconds(start,stop) << std::endl;
	    WMESH_STATUS_CHECK(status);
	  }
	

	WMESH_STATUS_CHECK(status);
    
    if (self_->m_topology_dimension==3)
      {

      std::cout << "boundary faces (triangles, quadrilaterals):  "
		<< num_bfacets[0]
		<< ", "
		<< num_bfacets[1]
		<< std::endl;
    
      std::cout << "total num faces (triangles, quadrilaterals):  "
		<< num_facets[0]
		<< ", "
		<< num_facets[1] << std::endl;
    
      std::cout << "total num edges:  " << self_->m_num_edges  << std::endl;
    }
    else if (self_->m_topology_dimension==2)
      {

      std::cout << "boundary edges:  "
		<< num_bfacets[0]
		<< std::endl;
    
      std::cout << "total num edges:  "
		<< num_facets[0]
		<< std::endl;
    
      std::cout << "total num edges:  " << self_->m_num_edges  << std::endl;
    }
    return WMESH_STATUS_SUCCESS;
  }
  


};
