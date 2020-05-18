#include "wmesh.h"
#include "cmdline.hpp"

int main(int argc, char ** argv)
{
  
  WCOMMON::cmdline cmd(argc,argv);

  wmesh_t* mesh;
  wmesh_status_t status;
  wmesh_int_t nparts = 0;
  //
  // Get the output file name.
  //
  char ofilename[256];

  //
  // Get verbose.
  //
  bool verbose = cmd.option("-v");
  if (verbose)
    {
      cmd.disp_header(stdout);
    }     

  //
  // Get the number of partitions.
  //
  if (false == cmd.option("-n", &nparts))
    {
      fprintf(stderr,"missing output file, '-n' option.\n");
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
  
  const char * filename = cmd.get_arg(1);

  //
  // Read the mesh.
  //
  status = wmesh_read(&mesh,
		      filename);
  WMESH_STATUS_CHECK(status);

  //
  // Mesh analysis.
  //
  status = wmesh_partitioning(mesh, nparts);
  WMESH_STATUS_CHECK(status);

  
  //
  // Create the mesh.
  //
  status = wmesh_write(mesh,
		       ofilename);
  WMESH_STATUS_CHECK(status);
  
  return WMESH_STATUS_SUCCESS;
}
