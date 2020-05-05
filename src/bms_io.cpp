#include <array>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <valarray>
#include <iostream>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"

#include <chrono>
//#include "wfe_element.h"
#include "bms.h"
#include "libmeshb7.h"
using namespace std::chrono;

extern "C"
{
  wmesh_status_t wmesh_write_medit(const wmesh_t* 		self_,
				   const char * 		filename_,
				   ...)
  {
    WMESH_POINTER_CHECK(self_);
    WMESH_POINTER_CHECK(filename_);


    wmesh_status_t status;    
    wmesh_str_t filename;
    
    { va_list args;
      va_start(args,filename_);
      vsprintf(filename,filename_,args);
      va_end(args); }

    
    int64_t inm;

    status = bms_write_medit_open(filename,&inm);
    WMESH_STATUS_CHECK(status);
    
    //
    // Write the geometry.
    //
    wmesh_int_t 		coo_m  	= 3;
    wmesh_int_t 		coo_n  	= self_->m_num_nodes;
    wmesh_int_t 		coo_ld 	= 3;
    const double*__restrict__ 	coo_v 	= self_->m_coo;
    
    status = bms_write_medit_geometry(inm,
				      coo_m,
				      coo_n,
				      coo_v,
				      coo_ld,
				      self_->m_n_c.v,
				      self_->m_n_c.ld);
    WMESH_STATUS_CHECK(status);
    GmfSetKwd(inm,
	      GmfTetrahedra,
	      8);


    status = bms_write_medit_topology(inm,

				      self_->m_c2n.m_size,
				      self_->m_c2n.m_ptr,
				      self_->m_c2n.m_m,
				      self_->m_c2n.m_n,
				      self_->m_c2n.m_data,
				      self_->m_c2n.m_ld,
				      
				      self_->m_c_c.m_size,
				      self_->m_c_c.m_ptr,
				      self_->m_c_c.m_m,
				      self_->m_c_c.m_n,
				      self_->m_c_c.m_data,
				      self_->m_c_c.m_ld);    
    WMESH_STATUS_CHECK(status);

    
    if (self_->m_bf2n.m_size > 0)
      {
	status = bms_write_medit_topology(inm,					
					  self_->m_bf2n.m_size,
					  self_->m_bf2n.m_ptr,
					  self_->m_bf2n.m_m,
					  self_->m_bf2n.m_n,
					  self_->m_bf2n.m_data,
					  self_->m_bf2n.m_ld,
					  self_->m_bf_c.m_size,
					  self_->m_bf_c.m_ptr,
					  self_->m_bf_c.m_m,
					  self_->m_bf_c.m_n,
					  self_->m_bf_c.m_data,
					  self_->m_bf_c.m_ld);
	
	WMESH_STATUS_CHECK(status);
      }

    //
    // Close.
    //
    status = bms_medit_close(inm);
    WMESH_STATUS_CHECK(status);
    
    return WMESH_STATUS_SUCCESS;  
  };

  wmesh_status_t wmesh_read_medit(wmesh_t** 		self__,
				  const char * 		filename_,
				  ...)
  {

    wmesh_status_t status;
    char filename[256];    
    { va_list args;
      va_start(args,filename_);
      vsprintf(filename,filename_,args);
      va_end(args); }
    
    int32_t dim;
    int32_t version;
    int64_t inm;
    
    status =  bms_read_medit_open(filename,
				  &inm,
				  &version,
				  &dim);
    WMESH_STATUS_CHECK(status);

    
    wmesh_int_t
      num_nodes = 0,
      num_edges = 0,
      num_triangles = 0,
      num_quadrilaterals = 0,
      num_tetrahedra = 0,
      num_pyramids = 0,
      num_wedges = 0,
      num_hexahedra = 0;
    
    status =  bms_read_medit_stat(inm,
				  &num_nodes,
				  &num_edges,
				  &num_triangles,
				  &num_quadrilaterals,
				  &num_tetrahedra,
				  &num_pyramids,
				  &num_wedges,
				  &num_hexahedra);

    
    WMESH_STATUS_CHECK(status);

    //    std::cout << "numnodes " << num_nodes << std::endl;

    wmesh_int_t 	coo_m 		= dim;
    wmesh_int_t 	coo_n 		= num_nodes;
    wmesh_int_t 	coo_ld 		= coo_m;
    double * 		coo_v   	= nullptr;
  
    wmesh_int_t 	nflags_ld 	= 1;  
    wmesh_int_p   	nflags_v;

    wmesh_int_t c2n_size = 0;
    wmesh_int_t c2n_ptr[4+1];
    wmesh_int_t c2n_m[4]{4,5,6,8};
    wmesh_int_t c2n_n[4];
    wmesh_int_t c2n_ld[4];
    wmesh_int_p c2n_v = nullptr;

    wmesh_int_t c_c_size = 0;
    wmesh_int_t c_c_ptr[4+1];
    wmesh_int_t c_c_m[4]{1,1,1,1};
    wmesh_int_t c_c_n[4];
    wmesh_int_t c_c_ld[4];
    wmesh_int_p c_c_v = nullptr;
  
    c2n_n[0] = num_tetrahedra;
    c2n_n[1] = num_pyramids;
    c2n_n[2] = num_wedges;
    c2n_n[3] = num_hexahedra;

    wmesh_int_t bf2n_size = 0;
    wmesh_int_t bf2n_ptr[2+1];
    wmesh_int_t bf2n_m[2]{3,4};
    wmesh_int_t bf2n_n[2];
    wmesh_int_t bf2n_ld[2];
    wmesh_int_p bf2n_v = nullptr;


    wmesh_int_t bf_c_size = 0;
    wmesh_int_t bf_c_ptr[4+1];
    wmesh_int_t bf_c_m[4]{1,1,1,1};
    wmesh_int_t bf_c_n[4];
    wmesh_int_t bf_c_ld[4];
    wmesh_int_p bf_c_v = nullptr;

    wmesh_int_t total_num_bfaces = 0;
    wmesh_int_t total_num_cells = c2n_n[0] + c2n_n[1] + c2n_n[2] + c2n_n[3];
    if (total_num_cells > 0)
      {
	c2n_size = 4;

	bf2n_n[0] = num_triangles;
	bf2n_n[1] = num_quadrilaterals;
	total_num_bfaces = bf2n_n[0] + bf2n_n[1];
	if (total_num_bfaces > 0)
	  {
	    bf2n_size = 2;
	  }      
      }
    else
      {
	c2n_m[0] = 3;
	c2n_m[1] = 4;
      
	c2n_n[0] = num_triangles;
	c2n_n[1] = num_quadrilaterals;
	
	total_num_cells = c2n_n[0] + c2n_n[1];
	if (total_num_cells > 0)
	  {
	    c2n_size = 2;
	    bf2n_n[0] = num_edges;
	    total_num_bfaces = bf2n_n[0];
	    if (total_num_bfaces > 0)
	      {
		bf2n_m[0] = 2;
		bf2n_size = 1;
	      }
	  }
	else
	  {
	    c2n_n[0] = num_edges;
	    total_num_cells = c2n_n[0];
	    if (total_num_cells > 0)
	      {
		c2n_size = 1;	  
	      }
	    else
	      {
		return WMESH_STATUS_INVALID_CONFIG;
	      }
	  }
      }

    c_c_size = c2n_size;
    for (wmesh_int_t i=0;i<c2n_size;++i) c_c_n[i] = c2n_n[i];
    for (wmesh_int_t i=0;i<c2n_size;++i) c_c_ld[i] = c_c_m[i];
    c_c_ptr[0] = 0;
    for (wmesh_int_t i=0;i<c2n_size;++i) c_c_ptr[i+1] = c_c_ptr[i] + c_c_n[i] * c_c_ld[i];

    bf_c_size = bf2n_size;
    for (wmesh_int_t i=0;i<bf2n_size;++i) bf_c_n[i] = bf2n_n[i];
    for (wmesh_int_t i=0;i<bf2n_size;++i) bf_c_ld[i] = bf_c_m[i];
    bf_c_ptr[0] = 0;
    for (wmesh_int_t i=0;i<bf2n_size;++i) bf_c_ptr[i+1] = bf_c_ptr[i] + bf_c_n[i] * bf_c_ld[i];

  
    for (wmesh_int_t i=0;i<c2n_size;++i)
      {
	c2n_ld[i] = c2n_m[i];
      }

    c2n_ptr[0] = 0;
    for (wmesh_int_t i=0;i<c2n_size;++i)
      {
	c2n_ptr[i+1] = c2n_ptr[i] + c2n_ld[i] * c2n_n[i];
      }
  
    for (wmesh_int_t i=0;i<c2n_size;++i)
      {
	bf2n_ld[i] = bf2n_m[i];
      }
  
    bf2n_ptr[0] = 0;
    for (wmesh_int_t i=0;i<bf2n_size;++i)
      {
	bf2n_ptr[i+1] = bf2n_ptr[i] + bf2n_ld[i] * bf2n_n[i];
      }

    if (bf2n_size > 0)
      {
	bf2n_v 	= (wmesh_int_p)malloc(bf2n_ptr[bf2n_size] * sizeof(wmesh_int_t));
	bf_c_v 	= (wmesh_int_p)malloc(bf_c_ptr[bf_c_size] * sizeof(wmesh_int_t));
      }
  
    c2n_v 	= (wmesh_int_p)malloc(c2n_ptr[c2n_size] * sizeof(wmesh_int_t));
    c_c_v 	= (wmesh_int_p)malloc(c_c_ptr[c_c_size]*sizeof(wmesh_int_t));
  
    nflags_v 	= (wmesh_int_p)malloc(num_nodes*sizeof(wmesh_int_t));
    coo_v   	= (double*)malloc((3*num_nodes)*sizeof(double));
    if (!coo_v)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }
    
    status = bms_read_medit_geometry(inm,
				     coo_m,
				     coo_n,
				     coo_v,
				     coo_ld,
				     nflags_v,
				     nflags_ld);
    WMESH_STATUS_CHECK(status);

    status = bms_read_medit_topology(inm,
				   
				     c2n_size,
				     c2n_ptr,
				     c2n_m,
				     c2n_n,
				     c2n_v,
				     c2n_ld,
				   
				     c_c_size,
				     c_c_ptr,
				     c_c_m,
				     c_c_n,
				     c_c_v,
				     c_c_ld);
    WMESH_STATUS_CHECK(status);

    if (bf2n_size > 0)
      {      
	status = bms_read_medit_topology(inm,
				       
					 bf2n_size,
					 bf2n_ptr,
					 bf2n_m,
					 bf2n_n,
					 bf2n_v,
					 bf2n_ld,
					 
					 bf_c_size,
					 bf_c_ptr,
					 bf_c_m,
					 bf_c_n,
					 bf_c_v,
					 bf_c_ld);
	WMESH_STATUS_CHECK(status);      
      }
  
    status = bms_medit_close(inm);
    WMESH_STATUS_CHECK(status);      
    
    status =  wmesh_factory(self__,
			    num_nodes,

			    c2n_size,
			    c2n_ptr,
			    c2n_m,
			    c2n_n,
			    c2n_v,
			    c2n_ld,
			    
			    c_c_size,
			    c_c_ptr,
			    c_c_m,
			    c_c_n,
			    c_c_v,
			    c_c_ld,
				 
			    bf2n_size,
			    bf2n_ptr,
			    bf2n_m,
			    bf2n_n,
			    bf2n_v,
			    bf2n_ld,

			    bf_c_size,
			    bf_c_ptr,
			    bf_c_m,
			    bf_c_n,
			    bf_c_v,
			    bf_c_ld,				 
			    
			    dim,
			    num_nodes,
			    coo_v,
			    coo_ld,

			    nflags_v,
			    1);
    
#if 0
    status = wmesh_factory	(self__,
				 num_nodes,
				 c2n_size,
				 c2n_ptr,
				 c2n_m,
				 c2n_n,
				 c2n_v,
				 c2n_ld,
				 c_c_v,
				 1,//c_c_ld,
				 coo_v,
				 coo_ld,
				 nflags_v,
				 1, // nflags_ld,
				 bf2n_size,
				 bf2n_ptr,
				 bf2n_m,
				 bf2n_n,
				 bf2n_v,
				 bf2n_ld,
				 bf_c_v,
				 1);//bf_c_ld);

#if 0    
    wmesh_int_t c_c_size = self_->m_c2n.m_size;
    wmesh_int_t c_c_ld[4] = {1,1,1,1};
    wmesh_int_t c_c_ptr[5];
    wmesh_int_t c_c_m[4] = {1,1,1,1};
    wmesh_int_t c_c_n[4] = {0,0,0,0};
    for (wmesh_int_t i=0;i<self_->m_c2n.m_size;++i)
      {
	c_c_n[i] = self_->m_c2n.m_n[i];
      }
    c_c_ptr[0] = 0;
    for (wmesh_int_t i=0;i<self_->m_c2n.m_size;++i)
      {
	c_c_ptr[i+1] = c_c_ptr[i+1] + c_c_n[i] * c_c_ld[i];
      }
    
    wmesh_int_t bf_c_size = self_->m_bf2n.m_size;
    wmesh_int_t bf_c_ld[4] = {1,1,1,1};
    wmesh_int_t bf_c_ptr[5];
    wmesh_int_t bf_c_m[4] = {1,1,1,1};
    wmesh_int_t bf_c_n[4] = {0,0,0,0};
    wmesh_int_p bf_c_v = self_->m_flag_bfaces;	
    for (wmesh_int_t i=0;i<self_->m_bf2n.m_size;++i)
      {
	bf_c_n[i] = self_->m_bf2n.m_n[i];
      }
	
    bf_c_ptr[0] = 0;
    for (wmesh_int_t i=0;i<self_->m_bf2n.m_size;++i)
      {
	bf_c_ptr[i+1] = bf_c_ptr[i+1] + bf_c_n[i] * bf_c_ld[i];
      }
#endif
#endif
    return WMESH_STATUS_SUCCESS;  
  };


};
