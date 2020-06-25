
#include "app.hpp"

#include <iostream>
#include "wmesh-blas.hpp"
#include <algorithm>
#include "bms.h"
#include "wmesh.hpp"
#include "wmesh-enums.hpp"


template<typename T>
static  std::ostream& operator<<(std::ostream&out_,
				 const wmesh_mat_t<T>&that_)
{
  for (wmesh_int_t i=0;i<that_.m;++i)
    {
      for (wmesh_int_t j=0;j<that_.n;++j)
	{
	  out_ << " " << that_.v[that_.ld * j + i];
	}
      out_ << std::endl;
    }
  return out_;
};


int main(int argc, char ** argv)
{
#if 0  
  for(auto e : wmesh_enum_element_t::all)
    {
      std::cout << e << std::endl;
    }
  
  wmesh_enum_element_t a = wmesh_enum_element_t::TRIANGLE;
  std::cout << "ddddddd " << std::endl;
  std::cout << a << std::endl;
  std::cout << "ddddddd " << std::endl;
  a = wmesh_enum_element_t::WEDGE;
  std::cout << "ddddddd " << std::endl;
  std::cout << a << std::endl;
  std::cout << "ddddddd " << std::endl;
  std::cout << sizeof(wmesh_enum_element_t) << std::endl;
  std::cout << "ddddddd " << std::endl;
  wmesh_int_t h = 14;
std::cout << h << std::endl;
 h=a;
  std::cout << "ddddddd " << std::endl;
#endif  
  WCOMMON::cmdline::str_t 	ofilename;
  bool 				verbose 	= false;
  wmesh_int_t 			degree		= 0,element = 0, shape_family;
  wmesh_status_t status;
  {
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
    WCOMMON::cmdline::str_t element_name;
    if (false == cmd.option("-e", element_name))
      {
	fprintf(stderr,"// wmesh_restrict.tests::error: missing element name, '-e' option.\n");
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
    
    //
    // Get the degree of the dg space.
    //
    if (false == cmd.option("-d", &degree))
      {
	fprintf(stderr,"missing degree, '-d' option.\n");
	return WMESH_STATUS_INVALID_ARGUMENT;
      }


    
    WCOMMON::cmdline::str_t shape_family_name;
    if (false == cmd.option("-n", shape_family_name))
      {
	fprintf(stderr,"// bms_shape.tests::error: missing shape family name, '-n' option.\n");
	return WMESH_STATUS_INVALID_ARGUMENT;
      }
  
    status = app_str2shapefamily(shape_family_name,
				 &shape_family);
    if (status != WMESH_STATUS_SUCCESS)
      {
	fprintf(stderr,"wrong shape family name '%s'\n",shape_family_name);
	fprintf(stderr," available shape are:\n");
	for  (wmesh_int_t i=0;i<WMESH_SHAPE_FAMILY_ALL;++i)
	  {	  
	    status = app_shapefamily2str(i,
					 shape_family_name);
	    WMESH_STATUS_CHECK(status);
	    fprintf(stderr," - '%s'\n",shape_family_name);
	  }      
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

  }

  
  wmesh_shape_restrict_t<double> shape_restrict;
  status = wmesh_shape_restrict_def(&shape_restrict,
				    element,
				    shape_family,
				    degree);
  WMESH_STATUS_CHECK(status);

  status = wmesh_shape_restrict_info(shape_restrict);
  WMESH_STATUS_CHECK(status);

  return WMESH_STATUS_SUCCESS;
}
