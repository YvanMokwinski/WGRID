#include "app.hpp"
#include <math.h>
#include <stdlib.h>
#include "bms.h"
#include "wmesh-blas.h"
#include "wmesh-math.hpp"

#include <iostream>

void usage(const char * appname_)
{
  fprintf(stderr,"//\n");
  fprintf(stderr,"// %s -e <element> -d <integer>\n",appname_);
  fprintf(stderr,"//\n");
  fprintf(stderr,"// Testing quadrature.\n");
  fprintf(stderr,"//\n");
}

int main(int 		argc,
	 char** 	argv)
{
  wmesh_status_t status;
  wmesh_int_t degree;
  wmesh_int_t element;
  wmesh_int_t family;  

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

  WCOMMON::cmdline::str_t family_name;
  if (false == cmd.option("-f", family_name))
    {
      fprintf(stderr,"// error: missing nodes family name, '-f' option.\n");
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  
  status = app_str2cubaturefamily(family_name,
				  &family);
  if (status != WMESH_STATUS_SUCCESS)
    {
      fprintf(stderr,"wrong nodes family name '%s'\n",family_name);
      fprintf(stderr," available elements are:\n");
      for  (wmesh_int_t i=0;i<WMESH_CUBATURE_FAMILY_ALL;++i)
	{	  
	  status = app_cubaturefamily2str(i,
					  family_name);
	  WMESH_STATUS_CHECK(status);
	  fprintf(stderr," - '%s'\n",family_name);
	}      
      return WMESH_STATUS_INVALID_ARGUMENT;
    }
  
  
  wmesh_int_t topodim;
  status = bms_element2topodim(element,&topodim);
  WMESH_STATUS_CHECK(status);
  

  wmesh_int_t 	p_storage;
  wmesh_int_t 	p_m;
  wmesh_int_t 	p_n;
  double *  	p;
  wmesh_int_t 	p_ld;
  
  wmesh_int_t 	w_n;
  double *  	w;
  wmesh_int_t 	w_inc = 1;
  
  wmesh_int_t p1d_n;
  status = bms_cubature_num_nodes(WMESH_ELEMENT_EDGE,
				  family,
				  degree,
				  &p1d_n);
  WMESH_STATUS_CHECK(status);
  status = bms_cubature_num_nodes(element,
				  family,
				  degree,
				  &p_n);
    
  WMESH_STATUS_CHECK(status);

  w_n		= p_n;
  w_inc     	= 1;    

  p_storage 	= WMESH_STORAGE_INTERLEAVE;
  p_m         	= topodim;
  p_ld        	= p_m;

  p   	= (double*)malloc(sizeof(double)* p_n * p_m);
  w   	= (double*)malloc(sizeof(double)* p_n);
    
  wmesh_int_t rwork_n;
  status =  bms_cubature_buffer_size(element,
				     family,
				     p1d_n,		
				     &rwork_n);
  
  WMESH_STATUS_CHECK(status);
  double * rwork = (rwork_n > 0) ? (double*)malloc(sizeof(double)*rwork_n) : nullptr;
  if (rwork_n > 0 && !rwork)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
    }
  
  status = bms_cubature(element,
			family,
			p1d_n,
			  
			p_storage,
			p_m,
			p_n,
			p,
			p_ld,
			  
			w_n,
			w,
			w_inc,
			  
			rwork_n,
			rwork);

  double sum = 0.0;
  for (wmesh_int_t i=0;i<w_n;++i)
    sum+= w[w_inc*i];
  std::cerr << "sum " << sum << std::endl;
  std::cout << "MeshVersionFormatted" << std::endl;
  std::cout << "1" << std::endl;
  std::cout << "Dimension" << std::endl;
  std::cout << topodim << std::endl;
  std::cout << "Vertices" << std::endl;
  std::cout << p_n << std::endl;
  for (wmesh_int_t j=0;j<p_n;++j)
    {
      for (wmesh_int_t i=0;i<p_m;++i)
	{	    
	  std::cout << " " << p[p_ld*j+i];
	}
      std::cout << " 0" << std::endl;
    }
  std::cout << "End" << std::endl;
    
  if (rwork)
    {
      free(rwork);
      rwork = nullptr;
    }
  if (w)
    {
      free(w);
    }
  if (p)
    {
      free(p);
    }
  WMESH_STATUS_CHECK(status);
  return WMESH_STATUS_SUCCESS;
}
