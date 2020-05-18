#include <stdlib.h>
#include <iostream>
#include "wmesh.hpp"
#include "bms.h"
#include "wmesh_utils.hpp"

extern "C"
{  
  wmesh_status_t wmesh_rmacro_buffer_size(wmesh_int_t 	element_,
					  wmesh_int_t 	degree_,
					  wmesh_int_p	work_n_,
					  wmesh_int_p	num_entities_)
  {
    wmesh_status_t status;
    status = bms_rmacro_buffer_size(element_,
				    degree_,
				    work_n_,
				    num_entities_);
    WMESH_STATUS_CHECK(status);
    return WMESH_STATUS_SUCCESS;    
  }
  
  wmesh_status_t wmesh_rmacro_def(wmesh_t ** 	mesh__,
				  wmesh_int_t 	element_,
				  wmesh_int_t 	degree_)
  {
    
    WMESH_CHECK_POINTER(mesh__);
    wmesh_status_t status;
    
    wmesh_int_t
      topodim,
      numtypes,
      work_n,
      num_elements_of_topodim,
      elements_of_topodim[WMESH_ELEMENT_ALL],
      num_entities[WMESH_ELEMENT_ALL];
    

    //
    //
    //
    status = wmesh_element2topodim(element_,
				   &topodim);
    WMESH_STATUS_CHECK(status);

    //
    //
    //
    status = wmesh_topodim2numtypes(topodim,
				    &numtypes);
    WMESH_STATUS_CHECK(status);

    //
    //
    //
    status = wmesh_topodim2elements(topodim,
				    &num_elements_of_topodim,
				    elements_of_topodim);
    WMESH_STATUS_CHECK(status);


    //
    //
    //
    for (wmesh_int_t i=0;i<WMESH_ELEMENT_ALL;++i)
      {
	num_entities[i] = 0;
      }
    status = bms_rmacro_buffer_size(element_,
				    degree_,
				    &work_n,
				    num_entities);
    WMESH_STATUS_CHECK(status);
    
    wmesh_int_t num_nodes = num_entities[WMESH_ELEMENT_NODE];    
    wmesh_int_t c2n_size = num_elements_of_topodim;
    wmesh_int_t c2n_ptr[5];
    wmesh_int_t c2n_m[4];
    wmesh_int_t c2n_n[4];
    wmesh_int_t c2n_ld[4];
    wmesh_int_p c2n_v;    

        
    double * coo = (double*)malloc(sizeof(double) * num_nodes * topodim);
    if (!coo)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);	
      }

    
    //
    // Get the number of nodes per element.
    //
    status = wmesh_elements_num_nodes(c2n_size,
				      elements_of_topodim,
				      c2n_m);

    //
    // Get the elements.
    //
    for (wmesh_int_t i=0;i<c2n_size;++i)
      {
	c2n_n[i] = num_entities[elements_of_topodim[i]];	
      }

    status = wmesh_int_sparsemat_init(c2n_size,
				      c2n_ptr,
				      c2n_m,
				      c2n_n,
				      c2n_ld);
    WMESH_STATUS_CHECK(status);
    
    c2n_v = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*c2n_ptr[c2n_size]);
    if (!c2n_v)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }
    
    
    wmesh_int_p work = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*work_n);
    if (!work)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }	    

    wmesh_int_t icoo_n 		= num_nodes * topodim;
    wmesh_int_t work_n_rem 	= work_n - icoo_n;
    
    wmesh_int_t icoo_ld		= topodim;
    wmesh_int_p icoo 		= work;
    wmesh_int_p work_rem 	= work + icoo_n;
    
    status = bms_rmacro(element_,
			degree_,
			num_nodes,
			
			c2n_size,
			c2n_ptr,
			c2n_m,
			c2n_n,
			c2n_v,
			c2n_ld,
			
			topodim,
			num_nodes,
			icoo,
			icoo_ld,
			
			work_n_rem,
			work_rem);

    //
    //
    //
    wmesh_int_t family = 0;
    switch(family)
      {
      case 0:
	{
	  double idegree = ((double)1.0)/((double)degree_);
	  for (wmesh_int_t i=0;i<num_nodes*topodim;++i)
	    {
	      coo[i] = icoo[i]*idegree;
	    }
	  break;
	}
      default:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
	}
      }
    
#if 0
    //
    //
    //
    wmesh_int_t family = WFE_NODES_FAMILY_LAGRANGE;
    switch(family)
      {
      case WFE_NODES_FAMILY_LAGRANGE:
	{
	  double idegree = ((double)1.0)/((double)degree_);
	  for (wmesh_int_t i=0;i<num_nodes*topodim;++i)
	    {
	      coo[i] = icoo[i]*idegree;
	    }
	  break;
	}
	
      case WFE_NODES_FAMILY_LAGRANGEBUBBLE:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_NOT_IMPLEMENTED);
	  break;
	}
	
      case WFE_NODES_FAMILY_GAUSSLOBATTO:
	{
	  //
	  //
	  //
	  wmesh_nodes_icoo_gauss_lobatto(dim,
					 num_nodes,
					 icoo,
					 dim);
	  
	  //
	  //
	  //
	  const wmesh_int_t k_ld = (degree_+1)*(degree_+1);
	  const wmesh_int_t i_ld = (degree_+1);
	  wmesh_nodes_coo_gauss_lobatto(dim,
					num_nodes,
					coo,
					dim,
					dim,
					num_nodes,
					icoo,
					dim);
	  
	  for (wmesh_int_t idx=0;idx<num_nodes;++idx)
	    {
	      for (wmesh_int_t idim=0;idim<topodim;++idim)
		{
		  ijk[idim] = icoo[topodim*idx + idim];
		}
	      if (topodim==3)
		perm[NxN * ijk[2] + N * ijk[0] + ijk[1]] = 1+idx;
	      else if (topodim==2)
		perm[N * ijk[0] + ijk[1]] = 1+idx;
	      else if (topodim==1)
		perm[ijk[0]] = 1+idx;
	    }

	  //
	  // 
	  //	  
	  double idegree = ((double)1.0)/((double)degree_);
	  for (wmesh_int_t i=0;i<num_nodes*topodim;++i)
	    {	      
	      const wmesh_int_t ip	= icoo[3*l+0];
	      const wmesh_int_t jp	= icoo[3*l+1];
	      const wmesh_int_t kp	= icoo[3*l+2];
	      wmesh_int_t id 		= perm[k_ld*kp+i_ld*ip+jp];

	      //
	      //
	      //
	      coo[i] 			= coogl[i]*idegree;
	    }
	  
	  break;
	}
	
      }
#endif
    
    
    status = wmesh_def(mesh__,
		       topodim,
		       c2n_size,
		       c2n_ptr,
		       c2n_m,
		       c2n_n,
		       c2n_v,
		       c2n_ld,
		       topodim,
		       num_nodes,
		       coo,
		       topodim);
    
    return WMESH_STATUS_SUCCESS;
  };

}
