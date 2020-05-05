#include "wmesh.h"

#include <math.h>
#include <stdlib.h>
#include "cmdline.hpp"

void usage(const char * appname_)
{
  fprintf(stderr,"//\n");
  fprintf(stderr,"// %s  --nteta <value> --nradius <value> -o <filename>\n",appname_);
  fprintf(stderr,"//\n");
  fprintf(stderr,"// Generate a circular face mesh\n");
  fprintf(stderr,"// Example: %s  --nteta 5 --nradius 7 -o example.mesh\n",appname_);
  fprintf(stderr,"//\n");
}

int main(int 		argc,
	  char** 	argv)
{
  wmesh_t* 			mesh = nullptr;
  wmesh_status_t 		status;
  WCOMMON::cmdline::str_t 	ofilename;
  const char * 			ifilename 	= nullptr;
  bool 				verbose 	= false;

  wmesh_int_t
    opt_nteta = 0,
    opt_nradius = 0;
  
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
  verbose = cmd.option("-v");
  
  //
  // Get the number of partitions.
  //
  if (false == cmd.option("--nteta", &opt_nteta))
    {
      fprintf(stderr,"missing option '--nteta <integer-value>' option.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  else if (opt_nteta < 3)
    {
      fprintf(stderr,"invalid value from '--nteta <integer-value>' option, must be > 2.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  
  if (false == cmd.option("--nradius", &opt_nradius))
    {
      fprintf(stderr,"missing option '--nradius <integer-value>' option.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  else if (opt_nradius < 2)
    {
      fprintf(stderr,"invalid value from '--nradius <integer-value>' option, must be > 1.\n");
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

  
  double*x1d = new double[opt_nradius];
  for (wmesh_int_t i=0;i<opt_nradius;++i)
    {
      x1d[i] = ((double)i) / ((double)(opt_nradius-1));
    }
  status = wmesh_def_polar_extrusion(&mesh,
				     3,
				     &opt_nradius,
				     x1d,
				     &opt_nteta);
  delete [] x1d;
  WMESH_STATUS_CHECK(status);

  status = wmesh_write(mesh,
		       ofilename);
  WMESH_STATUS_CHECK(status);
  
  return WMESH_STATUS_SUCCESS;
}
