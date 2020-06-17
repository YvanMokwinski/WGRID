#pragma once

#include "bms_template_shape_lagrange_splz.hpp"

template<wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_template_shape_lagrange
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	degree_,
				     const_wmesh_int_p	diff_,				     
				     wmesh_int_t 	c_storage_,					  
				     wmesh_int_t 	c_m_,
				     wmesh_int_t 	c_n_,
				     const T * 	c_,
				     wmesh_int_t 	c_ld_,
		      
				     wmesh_int_t 	b_storage_,
				     wmesh_int_t 	b_m_,
				     wmesh_int_t 	b_n_,
				     T* 		b_,
				     wmesh_int_t 	b_ld_,
		      
				     wmesh_int_t   iw_n_,
				     wmesh_int_p   iw_,
				     wmesh_int_t   rw_n_,
				     T* rw_)
  {
#ifdef FORWARD
#error FORWARD already defined
#endif
#define FORWARD					\
    diff_,					\
      c_storage_,				\
      c_m_,					\
      c_n_,					\
      c_,					\
      c_ld_,					\
						\
      b_storage_,				\
      b_m_,					\
      b_n_,					\
      b_,					\
      b_ld_,					\
    						\
      iw_n_,					\
      iw_,					\
      rw_n_,					\
      rw_

    //    std::cout << "sssssssssssssssssssssss  " <<  degree_ << std::endl;
    switch(degree_)
      {
      case 0:
	{
	  return bms_template_shape_eval(FORWARD,
					 bms_template_shape_lagrange_splz<0,ELEMENT_,T>::basis);
	}
      case 1:
	{
	  return bms_template_shape_eval(FORWARD,
					 bms_template_shape_lagrange_splz<1,ELEMENT_,T>::basis);
	}
      default:
	{

	  //    std::cout << "sssssssssssssssssssssss GRAAL  " <<  degree_ << std::endl;
	  
	  wmesh_int_t dim;
	  wmesh_status_t status = bms_element2topodim(ELEMENT_,
				       &dim);
	  WMESH_STATUS_CHECK(status);
	  wmesh_int_t ev_n = (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_n_ : c_m_;
	  wmesh_int_t ndofs;
	  status = bms_ndofs(ELEMENT_,
			     degree_,
			     &ndofs);
	  WMESH_STATUS_CHECK(status);
	  
	  wmesh_int_t dof_storage = WMESH_STORAGE_INTERLEAVE;
	  wmesh_int_t dof_m 	= dim;
	  wmesh_int_t dof_n 	= ndofs;
	  wmesh_int_t dof_ld 	= dim;

	  wmesh_int_t v_storage = WMESH_STORAGE_INTERLEAVE;
	  wmesh_int_t v_m 	= ndofs;
	  wmesh_int_t v_n 	= ndofs;
	  wmesh_int_t v_ld 	= ndofs;
	  
	  wmesh_int_t o_storage = WMESH_STORAGE_INTERLEAVE;
	  wmesh_int_t o_m 	= dim;
	  wmesh_int_t o_n 	= ndofs;
	  wmesh_int_t o_ld 	= dim;
	  wmesh_int_p o_v    	= (wmesh_int_p)malloc(sizeof(wmesh_int_t) * o_n * o_ld);  
	  status = bms_ordering(ELEMENT_,
				degree_,
				o_storage,
				o_m,
				o_n,
				o_v,
				o_ld);
	  WMESH_STATUS_CHECK( status );
	  
	  T * dof    = (T*)malloc(sizeof(T) * dof_n * dof_ld);
	  T * v    = (T*)malloc(sizeof(T) * dof_n * dof_n);

	  wmesh_int_t rwork_n;
	  wmesh_int_t iwork_n;

	  
	  status = bms_nodes_buffer_sizes(ELEMENT_,
					  WMESH_NODES_FAMILY_LAGRANGE,
					  degree_,
					  &iwork_n,
					  &rwork_n);	  
	  WMESH_STATUS_CHECK(status);
	  
	  wmesh_int_p iwork = (iwork_n > 0) ? (wmesh_int_p)malloc(sizeof(wmesh_int_t)*iwork_n) : nullptr;
	  T* __restrict__ rwork = (rwork_n>0)?(T* __restrict__ )malloc(sizeof(T)*rwork_n):nullptr;
	  
	  status = bms_nodes(ELEMENT_,
			     WMESH_NODES_FAMILY_LAGRANGE,
			     degree_,
			     o_storage,
			     o_m,
			     o_n,
			     o_v,
			     o_ld,
			     dof_storage,
			     dof_m,
			     dof_n,
			     dof,
			     dof_ld,
			     iwork_n,
			     iwork,
			     rwork_n,
			     rwork);
	  
	  if (iwork) free(iwork);
	  if (rwork) free(rwork);
	  
	  WMESH_STATUS_CHECK( status );


	  if (ELEMENT_==WMESH_ELEMENT_QUADRILATERAL || ELEMENT_==WMESH_ELEMENT_HEXAHEDRON)
	    {
	      for (wmesh_int_t i=0;i<dof_n*dof_m;++i)
		{
		  dof[i] = 2.0*dof[i]-1.0;
		}
	    }
	  else
	    {

	    }
	  //
	  // Compute basis over dofs coordinates.
	  //
	  wmesh_int_t diff_vandermonde[3] = {0,0,0};
	  //	   std::cout << "yo0 " << std::endl;	  
	   status =  bms_template_shape_jacobi<ELEMENT_,T>::eval(degree_,
								 diff_vandermonde,
								 dof_storage,
								 dof_m,
								 dof_n,
								 dof,
								 dof_ld,
								 v_storage,
								 v_m,
								 v_n,
								 v,
								 v_ld,							     
								 iw_n_,
								 iw_,
								 rw_n_,
								 rw_);
	   WMESH_STATUS_CHECK( status );
	   //	   std::cout << "yo1 " << std::endl;	   
	   wmesh_int_t info_lapack;
#if 0
	   
	   for (wmesh_int_t i=0;i<c_m_;++i)
	     {
	       for (wmesh_int_t j=0;j<c_n_;++j)
		 {
		   fprintf(stdout," %e",c_[c_ld_*j+i]);
		 }
	       fprintf(stdout,"\n");	       
	     }
	   
	   fprintf(stdout,"------------------\n");	       

	   for (wmesh_int_t i=0;i<dof_m;++i)
	     {
	       for (wmesh_int_t j=0;j<dof_n;++j)
		 {
		   fprintf(stdout," %e",dof[dof_ld*j+i]);
		 }
	       fprintf(stdout,"\n");	       
	     }
	   
	   fprintf(stdout,"------------------\n");	       
	   

	   for (wmesh_int_t i=0;i<v_m;++i)
	     {
	       for (wmesh_int_t j=0;j<v_n;++j)
		 {
		   fprintf(stdout," %e",v[v_ld*j+i]);
		 }
	       fprintf(stdout,"\n");	       
	     }

	   fprintf(stdout,"------------------\n");	       
	   for (wmesh_int_t i=0;i<dof_m;++i)
	     {
	       for (wmesh_int_t j=0;j<dof_n;++j)
		 {
		   fprintf(stdout," %e",dof[dof_ld*j+i]);
		 }
	       fprintf(stdout,"\n");	       
	     }
	   
	   fprintf(stdout,"------------------\n");	       
#endif
	   
	   
	   status =  bms_template_shape_jacobi<ELEMENT_,T>::eval(degree_,
								 diff_,
								 c_storage_,
								 c_m_,
								 c_n_,
								 c_,
								 c_ld_,
								 b_storage_,
								 b_m_,
								 b_n_,
								 b_,
								 b_ld_,							     
								 iw_n_,
								 iw_,
								 rw_n_,
								 rw_);



	   

	   
	   WMESH_STATUS_CHECK( status );
#if 0
	   std::cout << "v_m " << v_m << std::endl;
	   std::cout << "v_ld " << v_ld << std::endl;
	   std::cout << "ev_n " << ev_n << std::endl;
	   std::cout << "b_ld_ " << b_ld_ << std::endl;
	   wmesh_int_t perm[32];
#endif
	   xgesv((wmesh_int_p)&v_m,
		 (wmesh_int_p)&ev_n,
		 v,
		 (wmesh_int_p)&v_ld,
		 iw_,
		 b_,
		 (wmesh_int_p)&b_ld_,
		 (wmesh_int_p)&info_lapack);
	   
	   //	   std::cout << " " << info_lapack << std::endl;
	   //	   	   for (wmesh_int_t l=0;l<v_m;++l)	   std::cout << " " << perm[l] << std::endl;
	   WMESH_CHECK( 0 == info_lapack );
#if 0
	   
	   for (wmesh_int_t i=0;i<c_m_;++i)
	     {
	       for (wmesh_int_t j=0;j<c_n_;++j)
		 {
		   fprintf(stdout," %e",c_[c_ld_*j+i]);
		 }
	       fprintf(stdout,"\n");	       
	     }
	   
	   fprintf(stdout,"------------------\n");	       

	   for (wmesh_int_t i=0;i<dof_m;++i)
	     {
	       for (wmesh_int_t j=0;j<dof_n;++j)
		 {
		   fprintf(stdout," %e",dof[dof_ld*j+i]);
		 }
	       fprintf(stdout,"\n");	       
	     }
	   
	   fprintf(stdout,"------------------\n");	       

	   for (wmesh_int_t i=0;i<v_m;++i)
	     {
	       for (wmesh_int_t j=0;j<v_n;++j)
		 {
		   fprintf(stdout," %e",v[v_ld*j+i]);
		 }
	       fprintf(stdout,"\n");	       
	     }
	   fprintf(stdout,"------------------\n");	       
	   

	   for (wmesh_int_t i=0;i<b_m_;++i)
	     {
	       for (wmesh_int_t j=0;j<b_n_;++j)
		 {
		   fprintf(stdout," %e",b_[b_ld_*j+i]);
		 }
	       fprintf(stdout,"\n");	       
	     }

	   fprintf(stdout,"------------------\n");	       
#endif

	   //
	   // Compute basis over c_
	   //
	   return WMESH_STATUS_SUCCESS;
	}
      }
    
    return WMESH_STATUS_INVALID_ARGUMENT;
#undef FORWARD
    
  }

};

