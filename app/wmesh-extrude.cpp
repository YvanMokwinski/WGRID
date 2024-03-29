#include "wmesh.h"


#include <math.h>
#include <stdlib.h>

#include "cmdline.hpp"

int main(int 		argc,
	  char** 	argv)
{
  wmesh_t* 		surface 	= nullptr;
  wmesh_t* 		mesh 		= nullptr;
  wmesh_status_t 	status;

  //
  // Parameters.
  //
  WCOMMON::cmdline::str_t 	ofilename;
  WCOMMON::cmdline::str_t 	ifilename_curve;
  const char * 			ifilename 	 = nullptr;
  bool 				verbose 	 = false;
  wmesh_int_t
    nbRotations = 0;
  
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
  // Get the number of partitions.
  //
  if (false == cmd.option("-n", &nbRotations))
    {
      fprintf(stderr,"missing face option, '-n <integer>' option.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }

  bool has_curve = cmd.option("-s", ifilename_curve);

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
  
  //
  // Read the mesh.
  //
  status = wmesh_read	(&surface,
			 ifilename);
  WMESH_STATUS_CHECK(status);
  
#if 0

  const double ctrlpts[] = {-2.0,-0.5,-1.0,
			    -1.0,0.5,0.0,
			    0.0625,0.6,0.0,
			    0.0625,4.3,0.0,
			    0.75,2.8,2.0,
			    2.0,0.5,2.0,
			    2.0,-0.5,2.0,
			    3.0,-1.5,4.0,
			    4.0,2.5,4.0,
			    5.0,1.0,4.0,
			    7.0,0.0,7.0};

  wmesh_bspline_t* spline;
  wmesh_bspline_ddef(&spline,
		     3,
		     11,
		     ctrlpts,
		     3);

  wmesh_bspline_write_medit(spline,"testcurve.mesh");
  
  /*    wmesh_int_t nxy[2]   = {431,11};  */
  
  static const wmesh_int_t 	n1d 		= 6;
  static const double 		x1d[] 		= {((double)0.0),((double)0.2),((double)0.4),((double)0.6),((double)0.8),((double)1.0)};
  static const wmesh_int_t 	nbRotations1d 	= 83;

  wmesh_t * face;
  status = wmesh_def_polar_extrusion(&face,
				     3,
				     &n1d,
				     x1d,
				     &nbRotations1d);
  WMESH_STATUS_CHECK(status);
  status = wmesh_write_medit(face,ofilename);
  WMESH_STATUS_CHECK(status);
#endif

  wmesh_int_t bface_ids[3] 	= {101,102,103};
  
  double dz[1]= {1.0/((double)nbRotations)};

  status =  wmesh_def_extrusion(&mesh,
				surface,
				nbRotations,
				1,
				dz,
				bface_ids);
  
  WMESH_STATUS_CHECK(status);

  wmesh_t* curve = nullptr;
  if (has_curve)
    {
      status = wmesh_read	(&curve,
				 ifilename_curve);
      WMESH_STATUS_CHECK(status);

      status = wmesh_curve_extrusion(mesh, surface, curve);
      WMESH_STATUS_CHECK(status);
    }
  
  status = wmesh_extract_boundary(mesh);
  WMESH_STATUS_CHECK(status);
  status = wmesh_write(mesh,ofilename);
  WMESH_STATUS_CHECK(status);

  return 0;
}
