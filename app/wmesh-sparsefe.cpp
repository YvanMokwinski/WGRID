#include <stdlib.h>
#include "wmesh.h"
#include "WCOMMON/cmdline.hpp"
#include <iostream>

int main(int argc, char ** argv)
{
  
  WCOMMON::cmdline cmd(argc,argv);

  wmesh_t* mesh;
  wmesh_status_t status;
  wmesh_int_t degree;
  //
  // Get the output file name.
  //
  char ofilename[256];
  if (false == cmd.option("-o", ofilename))
    {
      fprintf(stderr,"missing output file, '-o' option.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  
  if (false == cmd.option("-d", &degree))
    {
      fprintf(stderr,"missing output file, '-d' option.\n");
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

#if 0
  //
  // Mesh analysis.
  //
  status = wmesh_reorder(mesh);
  WMESH_STATUS_CHECK(status);
#endif

  wmesh_int_t csr_size 	= 0;
  wmesh_int_p csr_ptr 	= nullptr;
  wmesh_int_p csr_ind 	= nullptr;
  status = wmesh_fespace_endomorphism	(mesh,
					 degree,
					 &csr_size,
					 &csr_ptr,
					 &csr_ind);
  WMESH_STATUS_CHECK(status);
  for (wmesh_int_t i=0;i<csr_size;++i)
    {
      for (wmesh_int_t s=csr_ptr[i];s<csr_ptr[i+1];++s)
	{
	  wmesh_int_t j = csr_ind[s];
	  std::cout << csr_size - i << " " << j + 1 << std::endl;
	}
    }
#if 0  
  //
  // Create the mesh.
  //
  status = wmesh_write(mesh,
		       ofilename);
  WMESH_STATUS_CHECK(status);
#endif
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
