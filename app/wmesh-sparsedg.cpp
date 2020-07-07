
#include "app.hpp"
#include <iostream>
#include <math.h>
#include <algorithm>
#include "bms.h"
#include "wmeshspacedg.h"

#if 0
template<typename T>
static  std::ostream& operator<<(std::ostream&out_,
				 const wmesh_mat_t<T>&that_)
{
  for (wmesh_int_t i=0;i<that_.m;++i)
    {
      for (wmesh_int_t j=0;j<that_.n;++j)
	{
	  out_ << " " << that_.v[that_.ld * j + i];
	}
      out_ << std::endl;
    }
  return out_;
};
#endif


int main(int argc, char ** argv)
{
  wmesh_t* 		mesh 		= nullptr;
  wmesh_status_t 	status;
  
  //
  // Parameters.
  //
  WCOMMON::cmdline::str_t 	ofilename;
  const char * 			ifilename 	= nullptr;
  bool 				verbose 	= false;
  wmesh_int_t 			degree		= 0;
  
  {
    WCOMMON::cmdline cmd(argc,
			 argv);
    //
    // Get verbose.
    //
    verbose = cmd.option("-v");
    if (verbose)
      {
	cmd.disp_header(stdout);
      }     
   
    //
    // Get the degree of the dg space.
    //
    if (false == cmd.option("-d", &degree))
      {
	fprintf(stderr,"missing output file, '-d' option.\n");
	return WMESH_STATUS_INVALID_ARGUMENT;
      }
    
    //
    // Get output filename.
    //
    if (false == cmd.option("-o", ofilename))
      {
	fprintf(stderr,"missing output file, '-o' option.\n");
	return WMESH_STATUS_INVALID_ARGUMENT;
      }
    
    if (cmd.get_nargs() == 1)
      {
	fprintf(stderr,"no file found.\n");
	return WMESH_STATUS_INVALID_ARGUMENT;
      }
    
    ifilename = cmd.get_arg(1);
  }
  
  //
  // Read the mesh.
  //
  status = wmesh_read(&mesh,
		      ifilename);
  WMESH_STATUS_CHECK(status);


  //
  // Analyze the mesh.
  //
  status = wmesh_analysis(mesh);
  
  WMESH_STATUS_CHECK(status);

  wmesh_shape_info_t* shape_info_element;
  wmesh_shape_info_t* shape_info_trial;
  wmesh_shape_info_t* shape_info_u;
  wmesh_shape_info_t* shape_info_test;
  

  const wmesh_int_t shape_element_degree 		= 1;
  const wmesh_int_t shape_element_family 		= WMESH_SHAPE_FAMILY_LAGRANGE;
  //  const wmesh_int_t shape_element_nodes_family		= WMESH_NODES_FAMILY_LAGRANGE;

  const wmesh_int_t shape_u_degree 			= 2;
  const wmesh_int_t shape_u_family 			= WMESH_SHAPE_FAMILY_LAGRANGE;
  const wmesh_int_t shape_u_nodes_family 		= WMESH_NODES_FAMILY_LAGRANGE;

  const wmesh_int_t shape_trial_degree 			= degree;
  const wmesh_int_t shape_trial_family 		= WMESH_SHAPE_FAMILY_LAGRANGE;
  const wmesh_int_t shape_trial_nodes_family 		= WMESH_NODES_FAMILY_LAGRANGE;

  const wmesh_int_t shape_test_degree 		= degree;
  const wmesh_int_t shape_test_family 		= WMESH_SHAPE_FAMILY_LAGRANGE;
  //  const wmesh_int_t shape_test_nodes_family 		= WMESH_NODES_FAMILY_LAGRANGE;

  const wmesh_int_t cubature_family 		= WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE;
  const wmesh_int_t cubature_degree 		= shape_u_degree + shape_trial_degree + shape_test_degree;

  status = wmesh_shape_info_def(&shape_info_element,shape_element_family, shape_element_degree);
  WMESH_STATUS_CHECK(status);
  
  status = wmesh_shape_info_def(&shape_info_u,shape_u_family, shape_u_degree);
  WMESH_STATUS_CHECK(status);
  
  status = wmesh_shape_info_def(&shape_info_trial,shape_trial_family, shape_trial_degree);
  WMESH_STATUS_CHECK(status);
  
  status = wmesh_shape_info_def(&shape_info_test,shape_test_family, shape_test_degree);
  WMESH_STATUS_CHECK(status);

  wmesh_cubature_info_t * cubature_info;
  status = wmesh_cubature_info_def(&cubature_info,cubature_family, cubature_degree);
  WMESH_STATUS_CHECK(status);
  
  //
  // Generate the finite element space for the velocity.
  // 
  wmeshspace_t * velocity_space;
  status = wmeshspace_def(&velocity_space,
			  shape_u_nodes_family,
			  shape_u_degree,
			  mesh);
  WMESH_STATUS_CHECK(status);
  
  wmesh_int_t csr_size 	= 0;
  wmesh_int_p csr_ptr 	= nullptr;
  wmesh_int_p csr_ind 	= nullptr;
  wmeshspacedg_sparse(mesh,
		      degree,
		      &csr_size,
		      &csr_ptr,
		      &csr_ind);

  double *		rhs 	= (double*) malloc(sizeof(double) * csr_size);
  double *		csr_val = (double*) malloc(sizeof(double) * csr_ptr[csr_size]);
  
#if 0
  
  //
  // Spy the symbolic matrix.
  //
  std::cout << "size " << csr_size << std::endl;
  std::cout << "nnz  " << csr_ptr[csr_size] << std::endl;

  FILE * f = fopen("out.txt","w");
  for (wmesh_int_t i=0;i<csr_size;++i)
    {
      for (wmesh_int_t s=csr_ptr[i];s<csr_ptr[i+1];++s)
	{
	  wmesh_int_t j = csr_ind[s];
	  fprintf(f," " WMESH_INT_FORMAT " " WMESH_INT_FORMAT "\n",csr_size - i, j);
	  //	  std::cout << csr_size - i << " " << j << std::endl;
	}
    }
  fclose(f);
#endif

  wmesh_int_t 	topodim;
  status = wmesh_get_topodim(mesh,
			     &topodim);
  WMESH_STATUS_CHECK(status);
  
  const wmesh_int_t 	velocity_storage 	= WMESH_STORAGE_INTERLEAVE;
  const wmesh_int_t	velocity_m 		= topodim;

  wmesh_int_t	velocity_n;  
  status = wmeshspace_get_ndofs(velocity_space,
				&velocity_n);
  WMESH_STATUS_CHECK(status);

  const wmesh_int_t	velocity_ld 		= velocity_m;  
  double *		velocity 		= (double*)malloc(sizeof(double) * velocity_n * velocity_m);

  
  //
  // Interpolate the velocity.
  //
  {
    wmesh_int_t coo_dofs_m  	= topodim;
    wmesh_int_t coo_dofs_n  	= velocity_n;
    wmesh_int_t coo_dofs_ld 	= coo_dofs_m;
    double * 	coo_dofs 	= (double*)malloc(sizeof(double) * coo_dofs_n * coo_dofs_m);    
    status 			=  wmeshspace_generate_dcoodofs(velocity_space,
							       WMESH_STORAGE_INTERLEAVE,
							       coo_dofs_m,
							       coo_dofs_n,
							       coo_dofs,
							       coo_dofs_ld);
    WMESH_STATUS_CHECK(status);
    double x[3];
    double u[3];
    for (wmesh_int_t j=0;j<coo_dofs_n;++j)
      {
	
	for (wmesh_int_t i=0;i<coo_dofs_m;++i)
	  {
	    x[i] = coo_dofs[coo_dofs_ld*j+i];
	  }

	for (wmesh_int_t i=0;i<coo_dofs_m;++i)
	  {	    
	    u[i] = 0.0 * x[i];
	  }	
	u[0] = 1.0;
	
	if (velocity_storage == WMESH_STORAGE_INTERLEAVE)
	  {
	    for (wmesh_int_t i=0;i<coo_dofs_m;++i)
	      {
		velocity[velocity_ld * j + i] = u[i];
	      }
	  }
	else 
	  {
	    for (wmesh_int_t i=0;i<coo_dofs_m;++i)
	      {
		velocity[velocity_ld * i + j] = u[i];
	      }
	  }	
      }
    
    free(coo_dofs);    
  }

  
  wmeshspacedg_t * spacedg;            
  status = wmeshspacedg_def(&spacedg,
			    shape_trial_nodes_family,
			    shape_trial_degree,
			    mesh);  
  WMESH_STATUS_CHECK(status);
  
  status = wmeshspacedg_advection(spacedg,
				  cubature_info,
				  shape_info_element,
				  shape_info_trial,
				  shape_info_test,
				  shape_info_u,
				  
				  velocity_space,
				  velocity_storage,
				  
				  velocity_m,
				  velocity_n,
				  velocity,
				  velocity_ld,
				  
				  csr_size,
				  csr_ptr,
				  csr_ind,
				  csr_val,				  
				  rhs);
  WMESH_STATUS_CHECK(status);

  //
  // Output.
  //  
  status =  bms_matrix_market_dense_dwrite(csr_size,
					   1,
					   rhs,
					   csr_size,
					   "%s_rhs.mtx",
					   ofilename);
  WMESH_STATUS_CHECK(status);
  
  status =  bms_matrix_market_csr_dwrite(csr_size,
					 csr_size,
					 csr_ptr[csr_size],
					 csr_ptr,
					 csr_ind,
					 csr_val,
					 "%s.mtx",
					 ofilename);
  WMESH_STATUS_CHECK(status);
  
  return WMESH_STATUS_SUCCESS;
}
