
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"
#include "bms.h"

#include <chrono>
#include <iostream>

#define WINT_SPARSE_MAT_PARAMS(_f)				\
    _f.m_ptr,							\
      _f.m_m,							\
      _f.m_n,							\
      _f.m_data,						\
      _f.m_ld

using namespace std::chrono;
extern "C"
{
  


  

  wmesh_status_t wmeshspace_sublinearmesh(wmeshspace_t * 	self_,
					  wmesh_t ** 		mesh__)
  {

    wmesh_status_t status;
    double * refevals[4];
    for (wmesh_int_t l=0;l<4;++l)
      {

	refevals[l]	= nullptr;
	if (self_->m_mesh->m_c2n.m_n[l]>0)
	  {
#ifndef NDEBUG
	    //	    std::cerr << "// wmesh.ndebug.verbose: wmesh_analysis, refmesh celltype_=" << l << ", degree_ " << degree_ << std::endl;
#endif
	    //wmesh_write_medit(self_->m_patterns[l],"ahaha.mesh");
	    
	    if (l==0)
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->m_patterns[l]->m_num_nodes*4);
		for (wmesh_int_t k=0;k<self_->m_patterns[l]->m_num_nodes;++k)
		  {
		    double r = self_->m_patterns[l]->m_coo[3*k+0];
		    double s = self_->m_patterns[l]->m_coo[3*k+1];
		    double t = self_->m_patterns[l]->m_coo[3*k+2];
		    b[k*4+0] = ((double)1.0)-(r+s+t);
		    b[k*4+1] = r;
		    b[k*4+2] = s;
		    b[k*4+3] = t;
		  }		  
		refevals[l] = b;		  
	      }

	    if (l==1)
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->m_patterns[l]->m_num_nodes*5);
		for (wmesh_int_t k=0;k<self_->m_patterns[l]->m_num_nodes;++k)
		  {
		    double r = self_->m_patterns[l]->m_coo[3*k+0];
		    double s = self_->m_patterns[l]->m_coo[3*k+1];
		    double t = self_->m_patterns[l]->m_coo[3*k+2];
		    // std::cout << "rst = " << r << " " << s << " " << t << " " << std::endl;

#if 0
		    double one = 1.0;
		    b[k*5+0] = ( one - r )* ( one - s ) * ( one - t );
		    b[k*5+1] = ( ( r )* ( one - s ) * ( one - t ) );
		    b[k*5+2] = ( ( r )* ( s ) * ( one - t ) );
		    b[k*5+3] = ( one - r )* ( s ) * ( one - t );
		    b[k*5+4] = ( one - r )* ( one - s ) * ( t ) +
		      ( r )* ( one - s ) * ( t ) +
		      ( r )* ( s ) * ( t ) +
		      ( one - r )* ( s ) * ( t );
#endif

#if 0
		    r = 2.0*r-1.0;
		    s = 2.0*s-1.0;
		    t = 2.0*t-1.0;
		    double hr = (1.0+r)*( (1.0-t) / 2.0 ) -1.0;
		    double hs = (1.0+s)*( (1.0-t) / 2.0 ) -1.0;
		    double ht = t;
#endif
		    if (t<1.0)
		      {
			  
			double h = 1.0-t;
			  
			double phir0 = (h-r)/h;
			double phir1 = r/h;
		      
			double phis0 = (h-s)/h;
			double phis1 = s/h;
			  
			b[k*5+0] = phir0 * phis0 * (1.0-t);
			b[k*5+1] = phir1 * phis0 * (1.0-t);
			b[k*5+2] = phir1 * phis1 * (1.0-t);
			b[k*5+3] = phir0 * phis1 * (1.0-t);
			b[k*5+4] = t;

			b[k*5+0] = (1.0 - t  - r) * (1.0-t-s) / (1.0-t);
			b[k*5+1] = (r * (1.0-t-s)) / (1.0-t);
			b[k*5+2] = (r * s) / (1.0-t);
			b[k*5+3] = ((1.0 - t - r) * s) / (1.0-t);
			b[k*5+4] = t;
			  
#if 0
			r = 2.0*r-1.0;
			s = 2.0*s-1.0;
			t = 2.0*t-1.0;
			double c = (1.0-t)/2.0;
			b[k*5+0] = ( (c-t) * (c-s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			b[k*5+1] = ( (c+t) * (c-s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			b[k*5+2] = ( (c+t) * (c+s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			b[k*5+3] = ( (c-t) * (c+s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			b[k*5+4] = (t+1.0)/2.0;
			  

#else
			b[k*5+0] = (1.0 - t  - r) * (1.0-t-s) / (1.0-t);
			b[k*5+1] = (r * (1.0-t-s)) / (1.0-t);
			b[k*5+2] = (r * s) / (1.0-t);
			b[k*5+3] = ((1.0 - t - r) * s) / (1.0-t);
			b[k*5+4] = t;
#endif
		      }
		    else
		      {
			b[k*5+0] = 0.0;
			b[k*5+1] = 0.0;
			b[k*5+2] = 0.0;
			b[k*5+3] = 0.0;
			b[k*5+4] = 1.0;			  
		      }
		      
		    //
		    // Split into tetrahedra.
		    //
		    // 4
		    //
		    // 3 2
		    // 0 1
		    //
#if 0
		    if (r > s)
		      {
			b[k*5+0] = (1.0-r-s-t);
			b[k*5+1] = r;
			b[k*5+2] = s;
			b[k*5+3] = 0.0;
			b[k*5+4] = t;
		      }
		    else
		      {
			b[k*5+0] = (1.0-r-s-t);
			b[k*5+1] = 0.0;
			b[k*5+2] = r;
			b[k*5+3] = s;
			b[k*5+4] = t;
		      }
#endif
		      
		      
#if 0
		    b[k*5+0] = ((r+t)*(s+t))/(2.0*(1.0-t));
		    b[k*5+1] = -((r+t)*(s+1.0))/(2.0*(1.0-t));
		    b[k*5+2] = -((r+t)*(s+t))/(2.0*(1.0-t));
		    b[k*5+3] = ((r+1.0)*(s+1.0))/(2.0*(1.0-t));
		    b[k*5+4] = (1.0+t)/2.0;
		    r=hr;
		    s=hs;
		    t=ht;
		    b[k*5+0] = ((r+t)*(s+t))/(2.0*(1.0-t));
		    b[k*5+1] = -((r+t)*(s+1.0))/(2.0*(1.0-t));
		    b[k*5+2] = -((r+t)*(s+t))/(2.0*(1.0-t));
		    b[k*5+3] = ((r+1.0)*(s+1.0))/(2.0*(1.0-t));
		    b[k*5+4] = (1.0+t)/2.0;
		    b[k*5+0] = (1.0-r) * (1.0-s) * (1.0-t);
		    b[k*5+1] = (r) * (1.0-s) * (1.0-t);
		    b[k*5+2] = (r) * s * (1.0-t);
		    b[k*5+3] = (1.0-r) * s * (1.0-t);
		    b[k*5+4] = t;
#endif

#if 0
		    b[k*5+5] = ( one - r -t  )* ( one - s -t ) / (4.0 * (1-t));
		    b[k*5+1] = 
		      b[k*5+2] = ( ( r )* ( s ) * ( one - t ) );
		    b[k*5+3] = ( one - r )* ( s ) * ( one - t );
		    b[k*5+4] = ( one - r )* ( one - s ) * ( t );
#endif
		  }

		refevals[l] = b;		  
#if 0
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->m_patterns[l]->m_num_nodes*6);
		for (wmesh_int_t k=0;k<self_->m_patterns[l]->m_num_nodes;++k)
		  {
		    double r = self_->m_patterns[l]->m_coo[3*k+0];
		    double s = self_->m_patterns[l]->m_coo[3*k+1];
		    double t = self_->m_patterns[l]->m_coo[3*k+2];
		    double one=1.0;
		    b[k*8+0] = ( one - r - s) * ( one - t );
		    b[k*8+1] = r * ( one - t );
		    b[k*8+2] = s * ( one - t );
		    b[k*8+4] = ( one - r - s) * ( t );
		    b[k*8+5] = r * ( t );
		  }		  
		refevals[l] = b;
#endif
	      }
	      
	    if (l==2)
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->m_patterns[l]->m_num_nodes*6);
		for (wmesh_int_t k=0;k<self_->m_patterns[l]->m_num_nodes;++k)
		  {
		    double r = self_->m_patterns[l]->m_coo[3*k+0];
		    double s = self_->m_patterns[l]->m_coo[3*k+1];
		    double t = self_->m_patterns[l]->m_coo[3*k+2];
		    double one=1.0;
		    b[k*6+0] = ( one - r - s) * ( one - t );
		    b[k*6+1] = r * ( one - t );
		    b[k*6+2] = s * ( one - t );
		    b[k*6+3] = ( one - r - s) * ( t );
		    b[k*6+4] = r * ( t );
		    b[k*6+5] = s * ( t );
		  }		  
		refevals[l] = b;		  
	      }

	    if (l==3)
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->m_patterns[l]->m_num_nodes*8);
		//		  std::cout << " " << std::endl;
		for (wmesh_int_t k=0;k<self_->m_patterns[l]->m_num_nodes;++k)
		  {
		    double r = self_->m_patterns[l]->m_coo[3*k+0];
		    double s = self_->m_patterns[l]->m_coo[3*k+1];
		    double t = self_->m_patterns[l]->m_coo[3*k+2];
		    //		      std::cout << "rst = " << r << " " << s << " " << t << " " << std::endl;
		    double one=1.0;
		    b[k*8+0] = ( one - r )* ( one - s ) * ( one - t );
		    b[k*8+1] = ( ( r )* ( one - s ) * ( one - t ) );
		    b[k*8+2] = ( ( r )* ( s ) * ( one - t ) );
		    b[k*8+3] = ( one - r )* ( s ) * ( one - t );
		    b[k*8+4] = ( one - r )* ( one - s ) * ( t );
		    b[k*8+5] = ( ( r )* ( one - s ) * ( t ) );
		    b[k*8+6] = ( ( r )* ( s ) * ( t ) );
		    b[k*8+7] = ( one - r )* ( s ) * ( t );
		  }

		refevals[l] = b;		  
	      }
	      
	  }

      }

#if 1
    
    double * coo_dofs = (double*)malloc(3*sizeof(double)*self_->m_ndofs);
    
    {
      double cell_xyz[32];
      wmesh_int_t cell_ld=3;
      // double * 	c_xyz 	= (double*)malloc(sizeof(double)*2048);
      // wmesh_int_p 	c_dofs 	= (wmesh_int_p)malloc(sizeof(wmesh_int_t)*2048);
      for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++l)
	{
	  double * 	refeval = refevals[l];

	    //
	    // Local c2d.
	    //	    
	    // wmesh_int_t c2d_n	= self_->m_c2d.m_n[l];
	    wmesh_int_t c2d_m 		= self_->m_c2d.m_m[l];
	    wmesh_int_t c2d_ld 		= self_->m_c2d.m_ld[l];
	    wmesh_int_p c2d_v 		= self_->m_c2d.m_data + self_->m_c2d.m_ptr[l];

	    //
	    // Local c2n.
	    //
	    wmesh_int_t c2n_n 		= self_->m_mesh->m_c2n.m_n[l];
	    wmesh_int_t c2n_m 		= self_->m_mesh->m_c2n.m_m[l];
	    wmesh_int_t c2n_ld		= self_->m_mesh->m_c2n.m_ld[l];
	    wmesh_int_p c2n_v 		= self_->m_mesh->m_c2n.m_data + self_->m_mesh->m_c2n.m_ptr[l];
	    
	    for (wmesh_int_t j=0;j<c2n_n;++j)
	      {
		
		//
		// Get the coordinates of the cell.
		//		
		for (wmesh_int_t i=0;i<c2n_m;++i)
		  {
		    wmesh_int_t idx = c2n_v[c2n_ld * j + i] - 1;
		    for (wmesh_int_t k=0;k<3;++k)
		      {
			cell_xyz[cell_ld * i + k] = self_->m_mesh->m_coo[3 * idx + k];
		      }
		  }
		
		
		//
		// Get the physical coordinates of the dofs.
		//
		//		std::cout << "c2d_m " << c2d_m << std::endl;
		for (wmesh_int_t i=0;i<c2d_m;++i)
		  {
		    double x = 0.0, y =0.0, z = 0.0;
		    for (wmesh_int_t k=0;k<c2n_m;++k)
		      x += refeval[c2n_m * i + k] * cell_xyz[cell_ld * k + 0];		    
		    for (wmesh_int_t k=0;k<c2n_m;++k)
		      y += refeval[c2n_m * i + k] * cell_xyz[cell_ld * k + 1];		    
		    for (wmesh_int_t k=0;k<c2n_m;++k)
		      z += refeval[c2n_m * i + k] * cell_xyz[cell_ld * k + 2];		    

		    wmesh_int_t idx = c2d_v[c2d_ld * j + i] - 1;
		    // std::cout << "idx " << idx << std::endl;
		    coo_dofs[3 * idx + 0] = x;
		    coo_dofs[3 * idx + 1] = y;
		    coo_dofs[3 * idx + 2] = z;
		  }
	      }
	  }
      }

      printf("generate codes.\n");

      wmesh_int_p dofs_cod = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*self_->m_ndofs);      
      for (wmesh_int_t i=0;i<self_->m_mesh->m_num_nodes;++i)
	{
	  dofs_cod[i] = self_->m_mesh->m_n_c.v[i];
	}
      for (wmesh_int_t i=self_->m_mesh->m_num_nodes;i<self_->m_ndofs;++i)
	{
	  dofs_cod[i] = 1001;
	}



      printf("generate sublinear connectivity.\n");
      {
	
	wmesh_int_t c2n_size = 4;
	wmesh_int_t c2n_m[4]{4,5,6,8};
	wmesh_int_t c2n_ld[4]{4,5,6,8};
	wmesh_int_t c2n_n[4]{0,0,0,0};	

	for (int i=0;i<4;++i)
	  {
	    if (self_->m_patterns[i]!=nullptr)
	      {
		if (i==1)
		  {
		    c2n_n[0] += self_->m_mesh->m_c2n.m_n[i] * self_->m_patterns[i]->m_c2n.m_n[0];
		    c2n_n[1] += self_->m_mesh->m_c2n.m_n[i] * self_->m_patterns[i]->m_c2n.m_n[1];
		  }
		else
		  {
		    c2n_n[i] = self_->m_mesh->m_c2n.m_n[i] * self_->m_patterns[i]->m_num_cells;
		  }
	      }
	  }
	

	wmesh_int_t c2n_ptr[5];
	c2n_ptr[0]=0;
	for (int i=0;i<4;++i)
	  {
	    //	    std::cout << c2n_ptr[i] << std::endl;
	    c2n_ptr[i+1] = c2n_ptr[i] + c2n_ld[i] * c2n_n[i];
	  }
	//	    std::cout << c2n_ptr[4] << std::endl;
	
	wmesh_int_p c2n_v = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*c2n_ptr[4]);
	
	wmesh_int_t mxdofs = 0;
	for (int i=0;i<4;++i)
	  //	  if (c2n_n[i] > 0) mxdofs = (self_->m_patterns[i]->m_num_nodes > mxdofs) ? self_->m_patterns[i]->m_num_nodes : mxdofs;
	  if (self_->m_patterns[i])
	  mxdofs = (self_->m_patterns[i]->m_num_nodes > mxdofs) ? self_->m_patterns[i]->m_num_nodes : mxdofs;
	wmesh_int_p dofs = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mxdofs);
	wmesh_int_p lidx = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mxdofs);

	printf("generate connectivity.\n");
	wmesh_int_t idx[4]{0,0,0,0};
	for (wmesh_int_t l=0;l<4;++l)
	  {
	    auto ref_c2n = &self_->m_patterns[l]->m_c2n;

	    wmesh_int_t ncells = self_->m_c2d.m_n[l];
	    for (wmesh_int_t j=0;j<ncells;++j)
	      {
		
		//
		// extract dofs.
		//
		for (wmesh_int_t i=0;i<self_->m_c2d.m_m[l];++i)
		  {
		    dofs[i] = self_->m_c2d.m_data[self_->m_c2d.m_ptr[l] + j * self_->m_c2d.m_ld[l] + i];
		  }

		
		//
		// For each.
		//
		for (wmesh_int_t ref_cell_type=0;ref_cell_type<ref_c2n->m_size;++ref_cell_type)
		  {
		    wmesh_int_t ref_ncells = ref_c2n->m_n[ref_cell_type];
		    wmesh_int_t ref_nnodes_per_cell = ref_c2n->m_m[ref_cell_type];
		    for (wmesh_int_t sj=0;sj<ref_ncells;++sj)
		      {
			
			for (wmesh_int_t si=0;si<ref_nnodes_per_cell;++si)
			  {
			    lidx[si] = ref_c2n->m_data[ref_c2n->m_ptr[ref_cell_type] + sj * ref_c2n->m_ld[ref_cell_type]  + si ] - 1;
			  }
			
			for (wmesh_int_t si=0;si<ref_c2n->m_m[ref_cell_type];++si)
			  {
			    c2n_v[ c2n_ptr[ref_cell_type] + c2n_ld[ref_cell_type] * idx[ref_cell_type] + si] = dofs[lidx[si]];
			  }
			++idx[ref_cell_type];		       
		      }
		  }
	      }
	  }
      printf("generate connectivity done.\n");

      status =  wmesh_def(mesh__,
			  self_->m_mesh->m_topology_dimension,				 
			  self_->m_ndofs,
			  c2n_size,
			  c2n_ptr,
			  c2n_m,
			  c2n_n,
			  c2n_v,
			  c2n_ld,
			  coo_dofs,
			  3);
      
      for (wmesh_int_t i=0;i<self_->m_ndofs;++i)
	{
	  mesh__[0]->m_n_c.v[i] = dofs_cod[i];
	}
      
      WMESH_STATUS_CHECK(status);
      //	printf("write mesh.\n");
      
      //	wmesh_write(space_sublinear,"roger.mesh");
      //	printf("write mesh done.\n");
      }
#endif
      return WMESH_STATUS_SUCCESS;
  }

};
