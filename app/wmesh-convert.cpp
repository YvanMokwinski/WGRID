#include "wmesh.h"
#include "cmdline.hpp"
#include <math.h>
#include <stdlib.h>

void usage(const char * appname_)
{
  fprintf(stderr,"//\n");
  fprintf(stderr,"// %s <filename> -o <filename>\n",appname_);
  fprintf(stderr,"//\n");
  fprintf(stderr,"// Convert a mesh file.\n");
  fprintf(stderr,"// Example: %s example.mesh -o example.vtk\n",appname_);
  fprintf(stderr,"//\n");
  fprintf(stderr,"// Note: this is mainly used as a parrot for debugging.\n");
  fprintf(stderr,"//\n");
}

int main(int 		argc,
	 char** 	argv)
{
  //
  // Command line.
  //
  WCOMMON::cmdline cmd(argc,
		       argv);

  if (cmd.get_nargs() == 1)
    {
      usage(argv[0]);
      return 1;
    }  
  else if ( cmd.option("-h") || cmd.option("--help") )
    {
      usage(argv[0]);
      return 0;
    }

  //
  // Get verbose.
  //
  bool verbose = cmd.option("-v");
  if (verbose)
    {
      cmd.disp_header(stdout);
    }    

  //
  // Get output filename.
  //
  WCOMMON::cmdline::str_t ofilename;
  if (false == cmd.option("-o", ofilename))
    {
      fprintf(stderr,"// wmesh-convert::error: missing output file, '-o' option.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }

  //
  // Get input filename.
  //
  if (cmd.get_nargs() == 1)
    {
      fprintf(stderr,"// wmesh-convert::error: no file found.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }  
  const char * ifilename = cmd.get_arg(1);

  //
  // Read the mesh.
  //
  wmesh_t* mesh;
  wmesh_status_t status = wmesh_read(&mesh,
				     ifilename);
  WMESH_STATUS_CHECK(status);
  
  //
  // Write the mesh.
  //
  status = wmesh_write(mesh,
		       ofilename);
  WMESH_STATUS_CHECK(status);

  //
  // Kill the mesh.
  //
  
  status 	= wmesh_kill(mesh);
  WMESH_STATUS_CHECK(status);
  mesh 		= nullptr;
  return 0;
}
