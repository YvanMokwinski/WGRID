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

  status = wmeshspace_laplace(meshspace,
			      csr_size,
			      csr_ptr,
			      csr_ind,
			      csr_val,
			      rhs);
  WMESH_STATUS_CHECK(status);

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
