#include "wmesh.hpp"
#include "bms.hpp"
extern "C"
{

  wmesh_status_t wmesh_def_rmacro(wmesh_t ** 	mesh__,
				  wmesh_int_t 	element_,
				  wmesh_int_t 	nodes_family_,
				  wmesh_int_t 	degree_)
  {
    
    WMESH_CHECK_POINTER(mesh__);
    wmesh_status_t status;
    wmesh_int_t 
      topodim,
      numtypes,
      num_elements_of_topodim,
      elements_of_topodim[WMESH_ELEMENT_ALL],
      num_entities[WMESH_ELEMENT_ALL];

    //
    // Get topology dimension.
    //
    status = bms_element2topodim(element_,
				 &topodim);
    WMESH_STATUS_CHECK(status);

    //
    // Get number of types.
    //
    status = bms_topodim2numtypes(topodim,
				    &numtypes);
    WMESH_STATUS_CHECK(status);

    //
    // Get elements.
    //
    status = bms_topodim2elements(topodim,
				  &num_elements_of_topodim,
				  elements_of_topodim);
    WMESH_STATUS_CHECK(status);


    
    for (wmesh_int_t i=0;i<WMESH_ELEMENT_ALL;++i)
      {
	num_entities[i] = 0;
      }

    wmesh_int_t iwork_n;
    status = bms_rmacro_buffer_size(element_,
				    degree_,
				    &iwork_n,
				    num_entities);
    WMESH_STATUS_CHECK(status);

    wmesh_int_t num_nodes = num_entities[WMESH_ELEMENT_NODE];    
    wmesh_int_t c2n_size = num_elements_of_topodim;
    wmesh_int_t c2n_ptr[5];
    wmesh_int_t c2n_m[4];
    wmesh_int_t c2n_n[4];
    wmesh_int_t c2n_ld[4];
    wmesh_int_p c2n_v;    

    wmesh_int_t c_storage = WMESH_STORAGE_INTERLEAVE;
    wmesh_int_t c_n = num_nodes;
    wmesh_int_t c_m = topodim;
    wmesh_int_t c_ld = c_m;
    double * c_v = (double*)malloc(sizeof(double) * c_m * c_n);
    if (!c_v)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);	
      }

    //
    // Get the number of nodes per element.
    //
    status = bms_elements_num_nodes(c2n_size,
				    elements_of_topodim,
				    c2n_m);

    //
    // Set c2n_n.
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
    
   
    wmesh_int_p iwork = (iwork_n > 0) ? (wmesh_int_p)malloc(sizeof(wmesh_int_t)*iwork_n) : nullptr;
    if (iwork_n > 0 && !iwork)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }	    

    const wmesh_int_t 	b_storage 	= WMESH_STORAGE_INTERLEAVE;
    const wmesh_int_t 	b_m 		= topodim;
    const wmesh_int_t 	b_n 		= num_nodes;
    const wmesh_int_t 	b_ld 		= b_m; 
    wmesh_int_p 	b_v		= (wmesh_int_p)malloc((b_m * b_n)*sizeof(wmesh_int_t));
    
    status = bms_rmacro(element_,
			degree_,
			num_nodes,
			
			c2n_size,
			c2n_ptr,
			c2n_m,
			c2n_n,
			c2n_v,
			c2n_ld,

			b_m,
			b_n,
			b_v,
			b_ld,
			
			iwork_n,
			iwork);
    
    WMESH_STATUS_CHECK(status);
        
    wmesh_int_t rwork_n;
    wmesh_int_t next_iwork_n;

    status = bms_nodes_buffer_sizes(element_,
				    nodes_family_,
				    degree_,
				    &next_iwork_n,
				    &rwork_n);
   
    WMESH_STATUS_CHECK(status);
    if (next_iwork_n > iwork_n)
      {
	iwork_n = next_iwork_n;
	if (iwork)
	  {
	    wmesh_int_p next_iwork = (wmesh_int_p)realloc(iwork,iwork_n*sizeof(wmesh_int_t));
	    if (!next_iwork)
	      {
		WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	      }
	    iwork = next_iwork;
	  }
	else
	  {
	    iwork = (wmesh_int_p)malloc(iwork_n*sizeof(wmesh_int_t));	    
	  }
      }
    
    double* __restrict__ rwork = (rwork_n>0)?(double* __restrict__ )malloc(sizeof(double)*rwork_n):nullptr;
    if (rwork_n > 0 && !rwork)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }	    

    
    status = bms_dnodes(element_,
			nodes_family_,
			degree_,
			
			b_storage,
			b_m,
			b_n,
			b_v,
			b_ld,
			
			c_storage,
			c_m,
			c_n,
			c_v,
			c_ld,

			iwork_n,
			iwork,
			rwork_n,
			rwork);

    WMESH_STATUS_CHECK(status);

    if (rwork) { free(rwork); } rwork_n = 0; rwork = nullptr;
    if (iwork) { free(iwork); } iwork_n = 0; iwork = nullptr;
    
    status = wmesh_def(mesh__,
		       topodim,
		       c2n_size,
		       c2n_ptr,
		       c2n_m,
		       c2n_n,
		       c2n_v,
		       c2n_ld,
		       c_m,
		       c_n,
		       c_v,
		       c_ld);



    
    //
    // Set some topological code for fun.
    //
    if (topodim==3)
      {
	for (wmesh_int_t k=0;k<mesh__[0]->m_c_c.m_size;++k)
	  {
	    for (wmesh_int_t j=0;j<mesh__[0]->m_c_c.m_n[k];++j)
	      {
		mesh__[0]->m_c_c.m_data[mesh__[0]->m_c_c.m_ptr[k] + mesh__[0]->m_c_c.m_ld[k]*j+0] = 4+k;
	      }
	  }
      }
    else if (topodim==2)
      {
	for (wmesh_int_t k=0;k<mesh__[0]->m_c_c.m_size;++k)
	  {
	    for (wmesh_int_t j=0;j<mesh__[0]->m_c_c.m_n[k];++j)
	      {
 		mesh__[0]->m_c_c.m_data[mesh__[0]->m_c_c.m_ptr[k] + mesh__[0]->m_c_c.m_ld[k]*j+0] = 2+k;
	      }
	  }
      }
    else if (topodim==1)
      {
	for (wmesh_int_t k=0;k<mesh__[0]->m_c_c.m_size;++k)
	  {
	    for (wmesh_int_t j=0;j<mesh__[0]->m_c_c.m_n[k];++j)
	      {
			mesh__[0]->m_c_c.m_data[mesh__[0]->m_c_c.m_ptr[k] + mesh__[0]->m_c_c.m_ld[k]*j+0] = 1+k;
	      }
	  }
      }

    //
    // Set the topological dimension of the freedom's support.
    //
    status = bms_ordering_topoid(element_,
				 degree_,
				 mesh__[0]->m_n_c.n,
				 mesh__[0]->m_n_c.v,
				 mesh__[0]->m_n_c.ld);
    WMESH_STATUS_CHECK( status );
    
    return WMESH_STATUS_SUCCESS;
  };

}
