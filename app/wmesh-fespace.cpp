#include "wmesh.h"
#include "WCOMMON/cmdline.hpp"

int main(int argc, char ** argv)
{  
  wmesh_t* 		mesh = nullptr;
  wmesh_t* 		refined_mesh = nullptr;
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
    
    //
    // Get the number of partitions.
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
  // Refine the mesh.
  //
  status = wmesh_refine(mesh,
			degree,
			&refined_mesh);
  
  WMESH_STATUS_CHECK(status);
  
  //
  // Write the refined mesh.
  //
  status = wmesh_write(refined_mesh,
		       ofilename);
  WMESH_STATUS_CHECK(status);

  
  return WMESH_STATUS_SUCCESS;
}
