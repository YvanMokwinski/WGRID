#include <stdlib.h>
#include "wmesh.h"
#include "cmdline.hpp"
#include <iostream>
#include "app.hpp"
#include "bms.h"

int main(int argc, char ** argv)
{
  
  WCOMMON::cmdline cmd(argc,argv);

  wmesh_t* mesh;
  wmesh_status_t status;
  wmesh_int_t degree;
  wmesh_int_t nodes_family;
  bool old = cmd.option("--old");
  //
  // Get the output file name.
  //
    if (false == cmd.option("-d", &degree))
    {
      fprintf(stderr,"missing output file, '-d' option.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }

  WCOMMON::cmdline::str_t nodes_family_name;
  if (false == cmd.option("-n", nodes_family_name))
    {
      fprintf(stderr,"// bms_nodes.tests::error: missing nodes family name, '-n' option.\n");
	return WMESH_STATUS_INVALID_ARGUMENT;
    }
  
  status = app_str2nodesfamily(nodes_family_name,
			       &nodes_family);
  if (status != WMESH_STATUS_SUCCESS)
    {
      fprintf(stderr,"wrong nodes family name '%s'\n",nodes_family_name);
      fprintf(stderr," available nodes are:\n");
      for  (wmesh_int_t i=0;i<WMESH_NODES_FAMILY_ALL;++i)
	{	  
	  status = app_nodesfamily2str(i,
				       nodes_family_name);
	  WMESH_STATUS_CHECK(status);
	  fprintf(stderr," - '%s'\n",nodes_family_name);
	}      
      return WMESH_STATUS_INVALID_ARGUMENT;
    }



  WCOMMON::cmdline::str_t ofilename;
  if (false == cmd.option("-o", ofilename))
    {
      fprintf(stderr,"// %s::error: missing output file name, '-o' option.\n",argv[0]);
      return WMESH_STATUS_INVALID_ARGUMENT;
    }

  
  if (cmd.get_nargs() == 1)
    {
      fprintf(stderr,"no file found.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  
  const char * filename = cmd.get_arg(1);

  //
  // Read the mesh.
  //
  status = wmesh_read(&mesh,
		      filename);
  WMESH_STATUS_CHECK(status);

  status = wmesh_analysis(mesh);
  WMESH_STATUS_CHECK(status);

  //
  // Define the finite element space.
  //
  wmeshspace_t * meshspace;
  status = wmeshspace_def(&meshspace,
			  nodes_family,
			  degree,
			  mesh);
  WMESH_STATUS_CHECK(status);
  //
  // Build the sublinear mesh.
  //
#if 0
  wmesh_t * refined_mesh;
  status = wmeshspace_sublinearmesh	(meshspace,
					 &refined_mesh);
  WMESH_STATUS_CHECK(status);
  
  //
  // Write the refined mesh.
  //
  status = wmesh_write			(refined_mesh,
					 "toto.mesh");
  WMESH_STATUS_CHECK(status);
#endif
  
  //
  // Get the sparsity pattern of a scalar equation.
  //
  std::cout << "create sparse" << std::endl;
  wmesh_int_t csr_size 	= 0;
  wmesh_int_p csr_ptr 	= nullptr;
  wmesh_int_p csr_ind 	= nullptr;
  status =  wmeshspace_sparse	(meshspace,
				 &csr_size,
				 &csr_ptr,
				 &csr_ind);
  WMESH_STATUS_CHECK(status);

  
  //
  // Spy the symbolic matrix.
  //
  std::cout << "size " << csr_size << std::endl;
  std::cout << "nnz  " << csr_ptr[csr_size] << std::endl;
#if 0
  FILE * f = fopen("out.txt","wb");
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

  //
  // Now compute Laplace equation.
  //

  
  //
  // Quadrature for each element.
  // Finite element evaluated on quadrature for each element
  // Loop and assembly.
  //

  double * csr_val 	= (double*)calloc(csr_ptr[csr_size],sizeof(double));
  double * rhs 		= (double*)calloc(csr_size,sizeof(double));
  std::cout << "compute laplace ... " << std::endl;


  if (old)
    {
  status =  wmeshspace_laplace_old(meshspace,
				   csr_size,
				   csr_ptr,
				   csr_ind,
				   csr_val,
				   rhs);
    }
  
  else
    {
  wmesh_cubature_info_t 	cubature_info;
  wmesh_shape_info_t		shape_info_element;	
  wmesh_shape_info_t		shape_info_trial;	
  wmesh_shape_info_t		shape_info_test;	
  wmesh_shape_info_t		shape_info_a;	

  wmesh_shape_info_def(&shape_info_a,WMESH_SHAPE_FAMILY_LAGRANGE,1);
  wmesh_shape_info_def(&shape_info_element,WMESH_SHAPE_FAMILY_LAGRANGE,1);
  wmesh_shape_info_def(&shape_info_test,WMESH_SHAPE_FAMILY_LAGRANGE, degree);
  wmesh_shape_info_def(&shape_info_trial,WMESH_SHAPE_FAMILY_LAGRANGE, degree);
  wmesh_cubature_info_def(&cubature_info,WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE, ( (degree-1) + (degree-1) + 1 + 1));
  
   status = wmeshspace_laplace(meshspace,
			       &cubature_info,
			       &shape_info_element,
			       &shape_info_trial,
			       &shape_info_test,
			       &shape_info_a,
			       csr_size,
			       csr_ptr,
			       csr_ind,
			       csr_val,
			       rhs);
    }


   WMESH_STATUS_CHECK(status);
  std::cout << "compute laplace done. " << std::endl;

  const char * extension = file_extension(ofilename);
      if ((nullptr != extension) && !strcmp(extension,".mtx"))
	{
	  wmesh_str_t obasename;
	  status = wmesh_basename(ofilename,obasename);
	  WMESH_STATUS_CHECK(status);
	  std::cout << "output ... " << obasename << std::endl;
	  status =  bms_matrix_market_dense_dwrite(csr_size,
						   1,
						   rhs,
						   csr_size,
						   "%s.rhs.mtx",
						   obasename);
	  WMESH_STATUS_CHECK(status);
	  
	  status =  bms_matrix_market_csr_dwrite(csr_size,
						 csr_size,
						 csr_ptr[csr_size],
						 csr_ptr,
						 csr_ind,
						 csr_val,
						 "%s",
						 ofilename);
	}
      else
	{
	  std::cout << "output ... " << std::endl;
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
      }
	WMESH_STATUS_CHECK(status);
  
  
  //
  // Save the linear system.
  //
  
  if (csr_ptr)
    {
      free(csr_ptr);
      csr_ptr = nullptr;
    }

  if (csr_ind)
    {
      free(csr_ind);
      csr_ind = nullptr;
    }

  return WMESH_STATUS_SUCCESS;
}
