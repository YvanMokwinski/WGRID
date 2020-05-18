#include "app.hpp"
#include <math.h>
#include <stdlib.h>
#include "bms.h"
#include "wmesh-blas.h"
#include "wmesh-math.hpp"
#include "wmesh_nodes.h"
#include "wmesh_nodes.h"
#include "wmesh_utils.hpp"

#include <iostream>
void usage(const char * appname_)
{
  fprintf(stderr,"//\n");
  fprintf(stderr,"// %s -e <element> -d <integer> -n <integer>\n",appname_);
  fprintf(stderr,"//\n");
}

int main(int 		argc,
	 char** 	argv)
{
  wmesh_status_t status;

  wmesh_int_t element;
  wmesh_int_t nodes_family;
  wmesh_int_t degree;
  
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

  if (false == cmd.option("-d", &degree))
    {
      fprintf(stderr,"missing option, '-d <integer>'.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  else if (degree < 0)
    {
      fprintf(stderr,"wrong value from option, '-d <integer>', must be >= 0.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;      
    }
  
  //
  // Get output filename.
  //
  WCOMMON::cmdline::str_t element_name;
  if (false == cmd.option("-e", element_name))
    {
      fprintf(stderr,"// bms_nodes.tests::error: missing element name, '-e' option.\n");
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
      fprintf(stderr," available elements are:\n");
      for  (wmesh_int_t i=0;i<WMESH_NODES_FAMILY_ALL;++i)
	{	  
	  status = app_nodesfamily2str(i,
				       nodes_family_name);
	  WMESH_STATUS_CHECK(status);
	  fprintf(stderr," - '%s'\n",nodes_family_name);
	}      
      return WMESH_STATUS_INVALID_ARGUMENT;
    }

  wmesh_int_t dim;
  status = wmesh_element2topodim(element,
				 &dim);
  WMESH_STATUS_CHECK(status);

  wmesh_int_t ndofs;
  status = wfe_ndofs(element,
		     degree,
		     &ndofs);
  WMESH_STATUS_CHECK(status);
  
  wmesh_int_t c_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_int_t c_m 	= dim;
  wmesh_int_t c_n 	= ndofs;
  wmesh_int_t c_ld 	= dim;

  wmesh_int_t o_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_int_t o_m 	= dim;
  wmesh_int_t o_n 	= ndofs;
  wmesh_int_t o_ld 	= dim;
  wmesh_int_p o_v    	= (wmesh_int_p)malloc(sizeof(wmesh_int_t) * o_n * o_ld);  
  status = bms_ordering(element,
			degree,
			o_storage,
			o_m,
			o_n,
			o_v,
			o_ld);
  WMESH_STATUS_CHECK( status );
  
  double * c    = (double*)malloc(sizeof(double) * c_n * c_ld);  
  wmesh_int_t work_n = 1024;
  double * work = (double*)malloc(sizeof(double) * work_n);  
  status = bms_dnodes(element,
		      nodes_family,
		      degree,
		      o_storage,
		      o_m,
		      o_n,
		      o_v,
		      o_ld,
		      c_storage,
		      c_m,
		      c_n,
		      c,
		      c_ld,
		      work_n,
		      work);
  free(work);  
  WMESH_STATUS_CHECK( status );

  //
  // Print.
  //
  std::cout << "MeshVersionFormatted" << std::endl;
  std::cout << "1" << std::endl;
  std::cout << "Dimension" << std::endl;
  std::cout << ((dim < 2) ? 2 : dim) << std::endl;
  std::cout << "Vertices" << std::endl;
  std::cout << ndofs << std::endl;
  for (wmesh_int_t j=0;j<c_n;++j)
    {
      for (wmesh_int_t i=0;i<c_m;++i)
	{
	  std::cout << " " << c[c_ld*j+i];
	}
      if (c_m == 1 )
	{
	  std::cout << " 0";
	}
      std::cout << " 0" << std::endl;
    }
  std::cout << " End" << std::endl;  
  return 0;
}
