#include "wmesh.h"
#include "cmdline.hpp"
#include <math.h>
#include <stdlib.h>
#include "bms.h"
#if 0
static const char * file_extension(const char * filename_)
{
  int i = -1;
  int len = 0;
  for (;;)
    {
      if (filename_[len] == '.')
	{
	  i = len;
	}
      if (filename_[len++] == '\0')
	{
	  break;
	}
    }
  return (i>=0) ? &filename_[i] : nullptr;
}
#endif

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

  const char * ifilename_extension = file_extension(ifilename);
  if (!ifilename_extension)
    {      
      fprintf(stderr,"// wmesh-convert::error: no extension found from file name '%s'.\n", ifilename);
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  
  if (!strcmp(ifilename_extension,".mesh") || !strcmp(ifilename_extension,".meshb") )
    {
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
      return WMESH_STATUS_SUCCESS;
    }
  
  if (!strcmp(ifilename_extension,".mtx") )
    {
      //
      // Read the mesh.
      //

      wmesh_int_t dense_m;
      wmesh_int_t dense_n;
      double * dense_v;
      wmesh_int_t dense_ld;
      
      wmesh_status_t status =  bms_matrix_market_dense_dread(&dense_m,
							     &dense_n,
							     &dense_v,
							     &dense_ld,
							     ifilename);
      WMESH_STATUS_CHECK(status);      

      FILE * f = fopen(ofilename,"w");
      fprintf(f,"MeshVersionFormatted\n1\nDimension\n2\n\n");
      fprintf(f,"HOSolutionAtQuadrilateralsQ1\n4\n1 1\n2 9" );
      for (wmesh_int_t i=0;i<dense_m;++i)
	{
	  for (wmesh_int_t j=0;j<dense_n;++j)
	    {
	      fprintf(f, " %8.15e", dense_v[dense_ld*j+i]);
	    }
	  fprintf(f,"\n");
	}
      fclose(f);

#if 0
      FILE * f = fopen(ofilename,"w");
      fprintf(f,"2 " WMESH_INT_FORMAT " " WMESH_INT_FORMAT " 2\n" , dense_n,dense_m);
      for (wmesh_int_t i=0;i<dense_m;++i)
	{
	  for (wmesh_int_t j=0;j<dense_n;++j)
	    {
	      fprintf(f, " %8.15e", dense_v[dense_ld*j+i]);
	    }
	  fprintf(f,"\n");
	}
      fclose(f);
#endif
      free(dense_v);


      
      return WMESH_STATUS_SUCCESS;
    }
  
}
