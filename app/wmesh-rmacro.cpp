#include "app.hpp"
#include <math.h>
#include <stdlib.h>

void usage(const char * appname_)
{
  fprintf(stderr,"//\n");
  fprintf(stderr,"// %s -s <element> -d <integer> -o <filename>\n",appname_);
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
  wmesh_status_t status;
  wmesh_int_t degree;
  wmesh_int_t element;

  
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

  if (false == cmd.option("-d", &degree))
    {
      fprintf(stderr,"missing option, '-d <integer>'.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  else if (degree < 1)
    {
      fprintf(stderr,"wrong value from option, '-d <integer>', must be >=1.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;      
    }
  
  //
  // Get output filename.
  //
  WCOMMON::cmdline::str_t element_name;
  if (false == cmd.option("-e", element_name))
    {
      fprintf(stderr,"// wmesh-convert::error: missing element name, '-e' option.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  


  status = app_str2element(element_name,
			   &element);
  if (status != WMESH_STATUS_SUCCESS)
    {
      fprintf(stderr,"wrong element name '%s'\n",element_name);
      fprintf(stderr," available elements are:\n");
      for  (wmesh_int_t i=WMESH_ELEMENT_EDGE;i<WMESH_ELEMENT_ALL;++i)
	{
	  
	  status = app_element2str(i,
				     element_name);
	  WMESH_STATUS_CHECK(status);
	  fprintf(stderr," - '%s'\n",element_name);
	}      
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  
  WCOMMON::cmdline::str_t ofilename;
  if (false == cmd.option("-o", ofilename))
    {
      fprintf(stderr,"// wmesh-convert::error: missing output file, '-o' option.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }


  //
  // Read the mesh.
  //
  wmesh_t* mesh;
  status = wmesh_rmacro_def(&mesh,
			    element,
			    degree);
  WMESH_STATUS_CHECK(status);
  
  //
  // Write the mesh.
  //
  status = wmesh_write(mesh,
		       ofilename);
  WMESH_STATUS_CHECK(status);

  return 0;
}
