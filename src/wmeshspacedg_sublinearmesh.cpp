//
// treat 2d
//
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-math.hpp"
#include "wmesh-status.h"
#include "wmesh.hpp"
#include "wmesh_utils.hpp"
#include "wmesh-blas.h"
#include "bms.h"

#include <chrono>
#include <iostream>
#include <math.h>
#include "bms_templates.hpp"

using namespace std::chrono;
extern "C"
{

  wmesh_status_t wmeshspacedg_generate_coodofs(wmeshspacedg_t * 	self_,
					       wmesh_int_t 	coo_storage_,
					       wmesh_int_t 	coo_m_,
					       wmesh_int_t 	coo_n_,
					       double * 	coo_,
					       wmesh_int_t 	coo_ld_)
  {    
    static constexpr double r0 = 0.0;
    static constexpr double r1 = 1.0;
    wmesh_status_t 	status;
    wmesh_t * 		mesh 	= self_->m_mesh;
    wmesh_int_t         topodim = mesh->m_topology_dimension;
    wmesh_int_t 	coo_m  	= mesh->m_coo_m;
    
    wmesh_int_t 	num_types;
    wmesh_int_t 	elements[4];
    double * 		refevals[4] {};
    double 		cell_xyz[32];
    wmesh_int_t 	cell_xyz_ld = coo_m;

    status = bms_topodim2elements(topodim,
				  &num_types,
				  elements);
    WMESH_STATUS_CHECK(status);        
    for (wmesh_int_t l=0;l<num_types;++l)
      {
	if (mesh->m_c2n.m_n[l]>0)
	  {
	    wmesh_t* 		rmacro			= self_->m_patterns[l];
	    const double * 	rmacro_coo 		= wmesh_get_coo(rmacro);
	    wmesh_int_t 	rmacro_coo_ld 		= rmacro->m_coo_ld;
	    wmesh_int_t 	rmacro_num_nodes 	= rmacro->m_num_nodes;
	    if ((topodim==0)&&(l==0))
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*rmacro_num_nodes*4);
		for (wmesh_int_t k=0;k<rmacro_num_nodes;++k)
		  {
		    double r = rmacro_coo[rmacro_coo_ld*k+0];
		    double s = rmacro_coo[rmacro_coo_ld*k+1];
		    double lr = ((double)1.0)-r;
		    double ls = ((double)1.0)-s;
		    b[k*4+0] = lr* ls;
		    b[k*4+1] = r * ls;
		    b[k*4+2] = r *  s;
		    b[k*4+3] = lr * s;
		  }		  
		refevals[l] = b;		  
	      }
	    if ((topodim==2)&&(l==0))
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*rmacro_num_nodes*3);
		for (wmesh_int_t k=0;k<rmacro_num_nodes;++k)
		  {
		    double r = rmacro_coo[rmacro_coo_ld*k+0];
		    double s = rmacro_coo[rmacro_coo_ld*k+1];
		    
		    b[k*3+0] = ((double)1.0)-(r+s);
		    b[k*3+1] = r;
		    b[k*3+2] = s;
		  }		  
		refevals[l] = b;		  
	      }
	    
	    if ((topodim==2)&&(l==1))
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*rmacro_num_nodes*4);
		for (wmesh_int_t k=0;k<rmacro_num_nodes;++k)
		  {
		    double r = rmacro_coo[rmacro_coo_ld*k+0];
		    double s = rmacro_coo[rmacro_coo_ld*k+1];
		    double lr = ((double)1.0)-r;
		    double ls = ((double)1.0)-s;
		    b[k*4+0] = lr* ls;
		    b[k*4+1] = r * ls;
		    b[k*4+2] = r *  s;
		    b[k*4+3] = lr * s;
		  }		  
		refevals[l] = b;		  
	      }

	    if ((topodim==3)&&(l==0))
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*rmacro_num_nodes*4);
		for (wmesh_int_t k=0;k<rmacro_num_nodes;++k)
		  {
		    double r = rmacro_coo[rmacro_coo_ld*k+0];
		    double s = rmacro_coo[rmacro_coo_ld*k+1];
		    double t = rmacro_coo[rmacro_coo_ld*k+2];
		    
		    b[k*4+0] = ((double)1.0)-(r+s+t);
		    b[k*4+1] = r;
		    b[k*4+2] = s;
		    b[k*4+3] = t;
		  }		  
		refevals[l] = b;		  
	      }

	    if ((topodim==3)&&(l==1))
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->m_patterns[l]->m_num_nodes*5);
		for (wmesh_int_t k=0;k<self_->m_patterns[l]->m_num_nodes;++k)
		  {
		    double r = rmacro_coo[rmacro_coo_ld*k+0];
		    double s = rmacro_coo[rmacro_coo_ld*k+1];
		    double t = rmacro_coo[rmacro_coo_ld*k+2];
		    // std::cout << "rst = " << r << " " << s << " " << t << " " << std::endl;
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
	      
	    if ((topodim==3)&&(l==2))
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
	    if ((topodim==3)&&(l==3))
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

    
    wmesh_int_t rw_n = 0;
    wmesh_int_t iw_n = 0;
    for (wmesh_int_t l=0;l<self_->m_mesh->m_c2n.m_size;++l)
      {
	wmesh_int_t k = self_->m_dofs_m[l]*topodim;
	rw_n = (rw_n < k) ? k : rw_n;
	iw_n = (iw_n < self_->m_dofs_m[l]) ? self_->m_dofs_m[l] : iw_n;
      }
    
    double * rw = (double*)malloc(sizeof(double)*rw_n);
    
    wmesh_int_p iw = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*iw_n);


    wmesh_int_t lc2d_inc 		= 1;
    wmesh_int_p lc2d 		= iw;
    for (wmesh_int_t l=0;l<self_->m_mesh->m_c2n.m_size;++l)
      {
	double * 	refeval = refevals[l];
	
	//
	// Local c2d.
	//	    
	// wmesh_int_t c2d_n	= self_->m_c2d.m_n[l];
	wmesh_int_t lc2d_n 		= self_->m_dofs_m[l];
	
	//
	// Local c2n.
	//
	wmesh_int_t c2n_n 		= mesh->m_c2n.m_n[l];
	wmesh_int_t c2n_m 		= mesh->m_c2n.m_m[l];
	//	wmesh_int_t c2n_ld		= mesh->m_c2n.m_ld[l];
	//	wmesh_int_p c2n_v 		= mesh->m_c2n.m_data + mesh->m_c2n.m_ptr[l];
	
	for (wmesh_int_t j=0;j<c2n_n;++j)
	  {
	    //
	    // Get the global indices of the dofs.
	    //
	    status = wmeshspacedg_get_dofs(self_, l, j, lc2d, lc2d_inc);
	    WMESH_STATUS_CHECK(status);    
#if 0
	    //
	    // Get the coordinates of the cell.
	    //		
	    for (wmesh_int_t i=0;i<c2n_m;++i)
	      {
		wmesh_int_t idx = c2n_v[c2n_ld * j + i] - 1;
		for (wmesh_int_t k=0;k<self_->m_mesh->m_coo_m;++k)
		  {
		    cell_xyz[cell_xyz_ld * i + k] = self_->m_mesh->m_coo[self_->m_mesh->m_coo_ld * idx + k];
		  }
	      }
#endif


	    status =  wmesh_get_cooelm(mesh,
				       l,
				       j,
				       WMESH_STORAGE_INTERLEAVE,
				       self_->m_mesh->m_coo_m,
				       c2n_m,
				       cell_xyz,				       
				       cell_xyz_ld);
	    WMESH_STATUS_CHECK(status);    

	    
	    dgemm("N",
		  "N",
		  &coo_m ,
		  &lc2d_n,
		  &c2n_m ,
		  &r1,
		  cell_xyz,
		  &cell_xyz_ld,
		  refeval,
		  &c2n_m,
		  &r0,
		  rw,
		  &coo_m);

	    //
	    // Copy.
	    //
	    if (WMESH_STORAGE_INTERLEAVE == coo_storage_)
	      {
		for (wmesh_int_t i=0;i<lc2d_n;++i)
		  {
		    wmesh_int_t idx = lc2d[lc2d_inc*i];
		    for (wmesh_int_t k = 0;k<coo_m;++k)
		      {
			coo_[coo_ld_ * idx + k] = rw[coo_m * i + k];
		      }
		  }
	      }
	    else
	      {
		for (wmesh_int_t i=0;i<lc2d_n;++i)
		  {
		    wmesh_int_t idx = lc2d[lc2d_inc*i];
		    for (wmesh_int_t k = 0;k<coo_m;++k)
		      {
			coo_[coo_ld_ * k + idx] = rw[coo_m * i + k];
		      }
		  }
	      }	    
	  }
      }
    for (wmesh_int_t l=0;l<num_types;++l)
      {
	if (refevals[l])
	  {
	    free(refevals[l]);
	  }
      }
    return WMESH_STATUS_SUCCESS;
  };



  wmesh_status_t wmeshspacedg_sublinearmesh(wmeshspacedg_t * 	self_,
					    wmesh_t ** 		mesh__)
  {    
    wmesh_status_t 	status;
    wmesh_t * 		mesh 	= self_->m_mesh;
    wmesh_int_t         topodim = mesh->m_topology_dimension;
    wmesh_int_t 	coo_m  	= mesh->m_coo_m;
    
    wmesh_int_t 	num_types;
    wmesh_int_t 	elements[4];
    double * 		refevals[4] {};
    double 		cell_xyz[32];
    wmesh_int_t 	cell_xyz_ld = coo_m;

    status = bms_topodim2elements(topodim,
				  &num_types,
				  elements);
    WMESH_STATUS_CHECK(status);    
    
    wmesh_int_t coo_dofs_m  	= coo_m;
    wmesh_int_t coo_dofs_n  	= self_->m_num_dofs;
    wmesh_int_t coo_dofs_ld 	= coo_dofs_m;
    double * 	coo_dofs 	= (double*)malloc(sizeof(double) * coo_dofs_n * coo_dofs_ld);
    
    status =  wmeshspacedg_generate_coodofs(self_,
					  WMESH_STORAGE_INTERLEAVE,
					  coo_dofs_m,
					  coo_dofs_n,
					  coo_dofs,
					  coo_dofs_ld);
    WMESH_STATUS_CHECK(status);

    //
    // GENERATE SUBLINEAR CONNECTIVITY
    //
    {
      
      wmesh_int_t c2n_size = num_types;
      wmesh_int_t c2n_ptr[5];
      wmesh_int_t c2n_m[4]{};
      wmesh_int_t c2n_n[4]{};	
      wmesh_int_t c2n_ld[4]{};
      
      status = bms_elements_num_nodes(num_types,
				      elements,
				      c2n_m);
      WMESH_STATUS_CHECK(status);
      
      for (int i=0;i<num_types;++i)
	{
	  if (self_->m_patterns[i]!=nullptr)
	    {
	      if (i==1 && topodim == 3)
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

      status = wmesh_int_sparsemat_init(c2n_size,
					c2n_ptr,
					c2n_m,
					c2n_n,
					c2n_ld);
      WMESH_STATUS_CHECK(status);      
      wmesh_int_p c2n_v = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*c2n_ptr[c2n_size]);
      if (!c2n_v)
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);      
	}

      wmesh_int_t mxdofs = 0;
      for (int i=0;i<num_types;++i)
	//	  if (c2n_n[i] > 0) mxdofs = (self_->m_patterns[i]->m_num_nodes > mxdofs) ? self_->m_patterns[i]->m_num_nodes : mxdofs;
	if (self_->m_patterns[i])
	  mxdofs = (self_->m_patterns[i]->m_num_nodes > mxdofs) ? self_->m_patterns[i]->m_num_nodes : mxdofs;
      wmesh_int_p dofs = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mxdofs);
      wmesh_int_p lidx = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mxdofs);

      //
      // Easy part.
      //
      wmesh_int_t idx[4]{0,0,0,0};
      for (wmesh_int_t l=0;l<num_types;++l)
	{
	  //	  wmesh_int_t lc2d_n = self_->m_dofs_m[l];
	  auto ref_c2n = &self_->m_patterns[l]->m_c2n;	  
	  wmesh_int_t ncells = mesh->m_c2n.m_n[l];

	  for (wmesh_int_t j=0;j<ncells;++j)
	    {
	      
	      status = wmeshspacedg_get_dofs(self_, l, j, dofs, 1);
	      WMESH_STATUS_CHECK(status);

#if 0
	      //
	      // extract dofs.
	      //
	      for (wmesh_int_t i=0;i<self_->m_c2d.m_m[l];++i)
		{
		  dofs[i] = self_->m_c2d.m_data[self_->m_c2d.m_ptr[l] + j * self_->m_c2d.m_ld[l] + i];
		}
#endif
	      
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
			  c2n_v[ c2n_ptr[ref_cell_type] + c2n_ld[ref_cell_type] * idx[ref_cell_type] + si] = dofs[lidx[si]] + 1;
			}
		      ++idx[ref_cell_type];
		    }
		}
	    }
	}
      
      //
      // Define the mesh.
      //
      status =  wmesh_def(mesh__,
			  self_->m_mesh->m_topology_dimension,				 
			  c2n_size,
			  c2n_ptr,
			  c2n_m,
			  c2n_n,
			  c2n_v,
			  c2n_ld,
			  coo_dofs_m,
			  coo_dofs_n,
			  coo_dofs,
			  coo_dofs_ld);

      //
      // Copy the dof codes.
      // It's a bit tricky, need to be changed.
      //
#if 0
      for (wmesh_int_t i=0;i<self_->m_ndofs;++i)
	{
	  mesh__[0]->m_n_c.v[i] = self_->m_dof_codes[i];
	}
#endif
      
      WMESH_STATUS_CHECK(status);
    }
    return WMESH_STATUS_SUCCESS;
  }

};




