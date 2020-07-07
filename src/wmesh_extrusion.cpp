#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include <chrono>
#include <iostream>
#include "wmesh_t.hpp"
#include "GenericEncoding.hpp"
#include "wmesh_utils.hpp"
#include <math.h>
#include "wmesh_bspline.h"
extern "C"
{

  wmesh_status_t  wmesh_spline_extrusion(wmesh_t *x,
					 const wmesh_t *	surface_,
					 const wmesh_bspline_t * spline)
  {
    wmesh_int_t n = surface_->m_num_nodes;
    return wmesh_bspline_transform_coordinates(spline,
					       x->m_num_nodes,
					       x->m_coo,
					       3,
					       x->m_n_c.v,
					       &n);
  };

  wmesh_status_t  wmesh_curve_extrusion(wmesh_t *	x,
					const wmesh_t *	surface_,
					const wmesh_t * curve_)
  {
    
    wmesh_int_t n = surface_->m_num_nodes;

    
    wmesh_bspline_t* spline;
    wmesh_bspline_ddef(&spline,
		       3,
		       curve_->m_num_nodes,
		       curve_->m_coo,
		       3);

    // wmesh_bspline_write_medit(spline,"testcurve.mesh");

    return wmesh_bspline_transform_coordinates(spline,
					       x->m_num_nodes,
					       x->m_coo,
					       3,
					       x->m_n_c.v,
					       &n);
    
  };

  wmesh_status_t wmesh_def_extrusion	(wmesh_t ** 			self__,
					 const wmesh_t * 		surface_,

					 wmesh_int_t 			nz_,
					 wmesh_int_t 			ndz_,
					 const double*__restrict__	dz_,
					 const_wmesh_int_p		bfaces_ids_)
  {
    wmesh_status_t status;
    const wmesh_int_t surface_numVertices 	= surface_->m_num_nodes;
    const wmesh_int_t nbVerticesOnOneFrame_	= surface_->m_num_nodes;    
    const wmesh_int_t nframes_ 			= nz_;
    const wmesh_int_t mesh3d_nbVertices  	= surface_numVertices * (nz_+1);
    
    wmesh_int_t c2n_size 	= 4;
    wmesh_int_t c2n_m	[4] 	= {4,5,6,8};
    wmesh_int_t c2n_ld	[4] 	= {4,5,6,8};
    wmesh_int_t c2n_n	[4];
    wmesh_int_t c2n_ptr	[5];
    c2n_n[0] = 0;
    c2n_n[1] = 0;
    c2n_n[2] = surface_->m_c2n.m_n[0] * nframes_;
    c2n_n[3] = surface_->m_c2n.m_n[1] * nframes_;
    
    c2n_ptr[0] = 0;
    for (wmesh_int_t  i = 0;i<4;++i)
      {
	c2n_ptr[i+1] = c2n_ptr[i] + c2n_n[i] * c2n_ld[i];
      }    

    wmesh_int_p c2n_v 			= (wmesh_int_p)malloc(sizeof(wmesh_int_t)*c2n_ptr[4]);
    double * 	node_coo 		= (double*)malloc(mesh3d_nbVertices * 3 * sizeof(double));
    const wmesh_int_t node_coo_ld 	= 3;    
    status =  wmesh_def	(self__,
			 3,				 
			 c2n_size,
			 c2n_ptr,
			 c2n_m,
			 c2n_n,
			 c2n_v,
			 c2n_ld,
			 3,
			 mesh3d_nbVertices,			 
			 node_coo,
			 node_coo_ld);
    WMESH_STATUS_CHECK(status);
    wmesh_t* self_ = self__[0];

    wmesh_int_p	node_ids 		= self_->m_n_c.v;
    const wmesh_int_t node_ids_ld 	= 1;
    
    double 	cooVertex[3];
    wmesh_int_t idTopology;
    
    double izValue = 0.0;
    for (wmesh_int_t iz=0;iz<=nz_;++iz)
      {
	for (wmesh_int_t ivertex=0;ivertex<surface_numVertices;++ivertex)
	  {
	    cooVertex[0] = surface_->m_coo[surface_->m_coo_ld*ivertex + 0];
	    cooVertex[1] = surface_->m_coo[surface_->m_coo_ld*ivertex + 1];
	    cooVertex[2] = izValue;
		
#if 1
#if 0
	    cooVertex[2] *= 2.0;
	    cooVertex[2] -= 1.0;
	    /* POUR LA SPLINE */
	    cooVertex[0] += 1.0;
	    cooVertex[1] += 1.0;
	    cooVertex[0] *= 0.5;
	    cooVertex[1] *= 0.5;
#endif
#endif
	    idTopology = surface_->m_n_c.v[ivertex * surface_->m_n_c.ld + 0];
	    node_coo[node_coo_ld* ( iz * surface_numVertices + ivertex ) + 0] = cooVertex[0];
	    node_coo[node_coo_ld* ( iz * surface_numVertices + ivertex ) + 1] = cooVertex[1];
	    node_coo[node_coo_ld* ( iz * surface_numVertices + ivertex ) + 2] = cooVertex[2];
		
	    if (idTopology==1)
	      {
		if (iz == 0)
		  {
		    node_ids[node_ids_ld * (iz * surface_numVertices + ivertex) + 0] = idTopology;
		  }
		else if (iz == nz_)
		  {
		    node_ids[node_ids_ld * (iz * surface_numVertices + ivertex) + 0] = idTopology;
		  }
		else
		  {
		    node_ids[node_ids_ld * (iz * surface_numVertices + ivertex) + 0] = bfaces_ids_[0];
		  }
	      }
	    else if (idTopology == 100)
	      {
		if (iz == 0)
		  {
		    node_ids[node_ids_ld * (iz * surface_numVertices + ivertex) + 0] = idTopology;
		  }
		else if (iz == nz_)
		  {
		    node_ids[node_ids_ld * (iz * surface_numVertices + ivertex) + 0] = bfaces_ids_[1];
		  }
		else
		  {
		    node_ids[node_ids_ld * (iz * surface_numVertices + ivertex) + 0] = 1000;
		  }
	      }		
	  } 	     	      
	izValue += (ndz_>1) ? dz_[iz] : dz_[0];
      }     
	
    
    wmesh_int_t f2n[4];
    wmesh_int_t v2n[8];
    wmesh_int_t face2volume[2] {2, 3};
    

    for (wmesh_int_t face_type=0;face_type<2;++face_type)
      {
	const wmesh_int_t nfaces 	= surface_->m_c2n.m_n[face_type];
	const wmesh_int_t nbNodesInFace = surface_->m_c2n.m_m[face_type];
	for (wmesh_int_t iframe=0;iframe<nframes_;++iframe)
	  {
	    for (wmesh_int_t iface=0;iface<nfaces;++iface)
	      {

		for (wmesh_int_t k=0;k<nbNodesInFace;++k)
		  {
		    f2n[k] = surface_->m_c2n.m_data[surface_->m_c2n.m_ptr[face_type] + surface_->m_c2n.m_ld[face_type] * iface + k] - 1;
		  }
		
		wmesh_int_t id = surface_->m_c_c.m_data[surface_->m_c_c.m_ptr[face_type] + iface * surface_->m_c_c.m_ld[face_type] + 0];
		for (wmesh_int_t inode=0;inode<nbNodesInFace;++inode)
		  {
		    v2n[inode] = f2n[inode] + iframe * nbVerticesOnOneFrame_;
		  } 
		
		for (wmesh_int_t inode=0;inode<nbNodesInFace;++inode)
		  {
		    v2n[nbNodesInFace + inode] = f2n[inode] + (iframe+1) * nbVerticesOnOneFrame_;
		  } 

		//
		// Copy volume information.
		//
		{
		  wmesh_int_t idx = nfaces * iframe + iface;
		  for (wmesh_int_t k=0;k<c2n_m[face2volume[face_type]];++k)
		    {		      
		      c2n_v[c2n_ptr[face2volume[face_type]] + c2n_ld[face2volume[face_type]] * idx + k] = v2n[k] + 1;
		    }
		  self_->m_c_c.m_data[self_->m_c_c.m_ptr[face2volume[face_type]] + iface * self_->m_c_c.m_ld[face2volume[face_type]] + 0] = id;
		} 
	      }
	  }
      }

    return WMESH_STATUS_SUCCESS;
  }

  wmesh_status_t wmesh_def_polar_extrusion(wmesh_t ** 		self_,
					   wmesh_int_t 		dim_,
					   const_wmesh_int_p 	n_,
					   const double * 	x_,
					   const_wmesh_int_p  	nbRotations_)
  {
    const wmesh_int_t topodim = 2;
    
    double cooVertex[3];
    wmesh_int_t cell2Nodes[4];
    double MNS_TWICEPI = acos(-1.0) * 2.0;
    const wmesh_int_t 	nbRotations		= nbRotations_[0];
    const wmesh_int_t 	n			= n_[0];
    const wmesh_int_t 	n_m1			= n-1;
    const wmesh_int_t 	numNodes		= (n-1) * nbRotations + 1;
    const wmesh_int_t 	numTriangles 		= nbRotations;
    const wmesh_int_t 	numQuadrilaterals 	= nbRotations * (n-2);
    
    wmesh_int_t c2n_size = 2;
    wmesh_int_t c2n_ptr[2+1] {0,numTriangles*3, numTriangles*3+numQuadrilaterals*4} ;
    wmesh_int_t c2n_m[2] {3,4} ;
    wmesh_int_t c2n_n[2] {numTriangles, numQuadrilaterals} ;
    wmesh_int_p c2n_v = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*c2n_ptr[2]);
    wmesh_int_t c2n_ld[2]{3,4};
    
    const double 	dtheta	= MNS_TWICEPI / ((double)nbRotations);
    
    double* coo = (double*)malloc(dim_*sizeof(double)*numNodes);
  
    //!
    //! @brief Get the number of entities with a specific dimension.
    //!
    wmesh_status_t status = wmesh_def	(self_,
					 topodim,
					 c2n_size,
					 c2n_ptr,
					 c2n_m,
					 c2n_n,
					 c2n_v,
					 c2n_ld,
					 dim_,
					 numNodes,
					 coo,
					 dim_);
    WMESH_STATUS_CHECK(status);
    
    wmesh_int_p		coo_c 		= self_[0]->m_n_c.v;
    wmesh_int_p		cell_c 		= self_[0]->m_c_c.m_data;
    
    /* COMPUTE THE GEOMETRY */
    cooVertex[2] = ((double)0.0);  
    
    wmesh_int_t vertexIndex = 0;
    cooVertex[0] = ((double)0.0);
    cooVertex[1] = ((double)0.0);
    {
      wmesh_int_t idx = vertexIndex++;  
      for(wmesh_int_t k=0;k<dim_;++k)
	{
	  coo[dim_*idx+k] = cooVertex[k];
	}
      coo_c[idx] = 100;
    }
    
  for (wmesh_int_t irot=0;irot<nbRotations;++irot)
    {
      const double theta = ((double)irot)*dtheta;
      for (wmesh_int_t i=1;i<n;++i)
	{
	  cooVertex[0] = x_[i] * cos(theta);
	  cooVertex[1] = x_[i] * sin(theta);
	    
	  {
	    wmesh_int_t idx = vertexIndex++;  
	    for(wmesh_int_t k=0;k<dim_;++k)
	      {
		coo[dim_*idx+k] = cooVertex[k];
	      }
	    coo_c[idx] = (i<n-1) ? 100 : 1;
	  }
	}	      
    }
  
  /* COMPUTE THE TRIANGLES */
  wmesh_int_t triangleIndex = 0;
  for (wmesh_int_t irot=0;irot<nbRotations;++irot)
    {
      const wmesh_int_t jrot = (irot+1)%nbRotations;
      
      cell2Nodes[0] = 0;
      cell2Nodes[1] = 1 + irot * n_m1;
      cell2Nodes[2] = 1 + jrot * n_m1;
      
      {
	wmesh_int_t idx = triangleIndex++;  
	for(wmesh_int_t k=0;k<3;++k)
	  {
	    c2n_v[c2n_ptr[0] + c2n_ld[0]*idx+k] = cell2Nodes[k]+1;
	  }
	cell_c[idx]=0;
      }
      
    }
  
  /* COMPUTE THE QUAD */
  wmesh_int_t quadrilateralIndex = 0;
  for (wmesh_int_t  irot=0;irot<nbRotations;++irot)
    {
      const wmesh_int_t jrot = (irot+1)%nbRotations;
      for (wmesh_int_t i=1;i<n_m1;++i)
	{
	  cell2Nodes[0] = i + irot * n_m1;
	  cell2Nodes[2] = (i+1) + jrot * n_m1;
	  
	  cell2Nodes[1] = cell2Nodes[0]+1;
	  cell2Nodes[3] = cell2Nodes[2]-1;
	  
	  {
	    wmesh_int_t idx = quadrilateralIndex++;  
	    for(wmesh_int_t k=0;k<4;++k)
	      {
		c2n_v[c2n_ptr[1] + c2n_ld[1]*idx+k] = cell2Nodes[k]+1;
	      }
	    cell_c[triangleIndex + idx]=1;
	  }
	  
	}
    } 
  
  return WMESH_STATUS_SUCCESS;
}
  
#ifdef __cplusplus
};
#endif
