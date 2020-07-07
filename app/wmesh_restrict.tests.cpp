
#include "app.hpp"

#include <iostream>
#include "wmesh-blas.hpp"
#include <algorithm>
#include "bms.h"
#include "wmesh_t.hpp"
#include "wmesh-enums.hpp"
#include "wmesh_shape_restrict_t.hpp"


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
#if 0
template<typename T>
wmesh_int_t wmesh_cubature_hash(wmesh_int_t element_,
				wmesh_int_t family_,
				wmesh_int_t degree_)
{
  //
  // element 0-7
  // family  0-7
  // degree  free
  //

}
#endif  
#if 0
static void PrintBits(size_t const size, void const * const ptr)
  {
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;
    printf("##\n");
    for (i=size-1;i>=0;i--)
      {
	for (j=7;j>=0;j--)
	  {
	    byte = b[i] & (1<<j);
	    byte >>= j;
	    printf("%u", byte);
	  }
      }
    puts("");
  };
#endif
#include "wmesh_cubature_factory_t.hpp"
#include "wmesh_nodes_factory_t.hpp"
#include "wmesh_shape_eval_factory_t.hpp"

int main(int argc, char ** argv)
{
#if 0
  wmesh_shape_info_t shape_info;
  wmesh_shape_info_def(&shape_info, WMESH_SHAPE_FAMILY_LAGRANGE, 1);
  
  wmesh_nodes_info_t nodes_info;
  wmesh_nodes_info_def(&nodes_info, WMESH_NODES_FAMILY_GAUSSLOBATTO, 3);
  
  {
    const wmesh_int_t element = WMESH_ELEMENT_TETRAHEDRON;
    auto nodes = wmesh_nodes_factory_t<double>::nodes_instance(element,							       
							       nodes_info);
    
    auto shape_eval_nodes = wmesh_shape_eval_factory_t<double>::shape_eval_instance(element,
										    shape_info,
										    nodes);
    std::cout << shape_eval_nodes->m_f << std::endl;
    std::cout << shape_eval_nodes->m_diff[0] << std::endl;
    std::cout << shape_eval_nodes->m_diff[1] << std::endl;
    shape_eval_nodes = wmesh_shape_eval_factory_t<double>::shape_eval_instance(element,
									       shape_info,
									       nodes);
    std::cout << shape_eval_nodes->m_f << std::endl;
    std::cout << shape_eval_nodes->m_diff[0] << std::endl;
    std::cout << shape_eval_nodes->m_diff[1] << std::endl;

    auto cubature = wmesh_cubature_factory_t<double>::cubature_instance(element,
									WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE,
									2);
    
    shape_eval_nodes = wmesh_shape_eval_factory_t<double>::shape_eval_instance(element,
									       shape_info,
									       cubature);
    std::cout << shape_eval_nodes->m_f << std::endl;
    shape_eval_nodes = wmesh_shape_eval_factory_t<double>::shape_eval_instance(element,
									       shape_info,
									       cubature);
    std::cout << shape_eval_nodes->m_f << std::endl;
    shape_eval_nodes = wmesh_shape_eval_factory_t<double>::shape_eval_instance(element,
									       shape_info,
									       cubature);
    std::cout << shape_eval_nodes->m_f << std::endl;
    shape_eval_nodes = wmesh_shape_eval_factory_t<double>::shape_eval_instance(element,
									       shape_info,
									       cubature);
    std::cout << shape_eval_nodes->m_f << std::endl;
    
  }
  
  exit(1);
 
  for (wmesh_int_t i=WMESH_ELEMENT_EDGE;i<WMESH_ELEMENT_ALL;++i)
    {
      auto nodes = wmesh_nodes_factory_t<double>::nodes_instance(i,
								 WMESH_NODES_FAMILY_GAUSSLOBATTO,
								 3);
      std::cout << "nodes first order" << std::endl;
      std::cout << nodes->m_c << std::endl;
    }
#endif
  
#if 0
  for (wmesh_int_t i=WMESH_ELEMENT_EDGE;i<WMESH_ELEMENT_ALL;++i)
    {
      auto nodes = wmesh_nodes_factory_t<double>::nodes_instance(i,
								 WMESH_NODES_FAMILY_LAGRANGE,
								 2);
      std::cout << "nodes first order" << std::endl;
      std::cout << nodes->m_c << std::endl;
    }
  
  for (wmesh_int_t i=WMESH_ELEMENT_EDGE;i<WMESH_ELEMENT_ALL;++i)
    {
  
      auto cubature2 = wmesh_cubature_factory_t<double>::cubature_instance(i,
									   WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE,
									   2);
      
      std::cout << cubature2->m_c << std::endl;
    }
  
  auto cubature3 = wmesh_cubature_factory_t<double>::cubature_instance(WMESH_ELEMENT_TETRAHEDRON,
								       WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE,
								       17);
  exit(1);
#endif  

#if 0  
  
  for (wmesh_int_t i=0;i<8;++i)
    {
      for (wmesh_int_t j=0;j<8;++j)
	{
	  for (wmesh_int_t k=0;k<64;++k)
	    {
	      //    std::cout << "(i,j,k) = ("<<i<<","<<j<<","<<k<<")" << std::endl;
	      wmesh_int_t encod;
	      wmesh_int_t status = wmesh_cubature_uniqueid_encod(encod,i,j,k);
	      auto ret = mymap.insert(std::pair<wmesh_int_t,wmesh_mat_t<double>*>(encod,(wmesh_mat_t<double>*)malloc(sizeof(wmesh_mat_t<double>))));
	      if (ret.second == false)
		{
		  std::cerr << "already registered" << std::endl;
		}
	      //	      std::cout << encod<<std::endl;
	      wmesh_int_t e,f,d;
	      status = wmesh_cubature_uniqueid_decod(encod,e,f,d);
	      WMESH_CHECK(e==i);
	      WMESH_CHECK(f==j);
	      WMESH_CHECK(d==k);
	    }
	}
    }
  wmesh_int_t a = 2;
  wmesh_int_t b = 4;
  wmesh_int_t c = 43;


  wmesh_int_t h = 2047;
  PrintBits(sizeof(wmesh_int_t),&h);
  wmesh_int_t dd = h >> 3 * 2;
  PrintBits(sizeof(wmesh_int_t),&dd);
  wmesh_int_t d = (~((wmesh_int_t)0) << 6)&h;
  PrintBits(sizeof(wmesh_int_t),&d);
  wmesh_int_t r = ~(~((wmesh_int_t)0) << 6)&h;
  PrintBits(sizeof(wmesh_int_t),&r);

  
  wmesh_int_t mask_element = ~(~((wmesh_int_t)0) << 3);
  PrintBits(sizeof(wmesh_int_t),&mask_element);
  wmesh_int_t mask_family = ~(~((wmesh_int_t)0) << 3) << 3;
  PrintBits(sizeof(wmesh_int_t),&mask_family);
  wmesh_int_t mask_degree = (~((wmesh_int_t)0) << 6);
  PrintBits(sizeof(wmesh_int_t),&mask_degree);

  wmesh_int_t val = 1394521;
  PrintBits(sizeof(wmesh_int_t),&val);

  {
    wmesh_int_t tmp = val & mask_element;
    PrintBits(sizeof(wmesh_int_t),&tmp);
  }
  {
    wmesh_int_t tmp = (val & mask_family) >> 3;
    PrintBits(sizeof(wmesh_int_t),&tmp);
  }
  {
    wmesh_int_t tmp = (val & mask_degree) >> 6;
    PrintBits(sizeof(wmesh_int_t),&tmp);
  }
  exit(1);
#endif
  //  wmesh_cubature_factory<double>& cubature_factory = wmesh_cubature_factory<double>::getInstance();


  
  wmesh_cubature_t<double> cubature;
  {
  wmesh_status_t status =  wmesh_cubature_def(&cubature,
					      WMESH_ELEMENT_TRIANGLE,
					      WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE,
					      4);
  WMESH_STATUS_CHECK(status);
  } 

  
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

#if 0  
  wmesh_shape_restrict_t<double> shape_restrict;
  status = wmesh_shape_restrict_def(&shape_restrict,
				    element,
				    shape_family,
				    degree);
  WMESH_STATUS_CHECK(status);

  status = wmesh_shape_restrict_info(shape_restrict);
  WMESH_STATUS_CHECK(status);
#endif

#if 0
  const wmesh_mat_t<double>* g = wmesh_calculate_shape_restrict<double>(WMESH_ELEMENT_TETRAHEDRON,
									0,
									0,
									
									WMESH_SHAPE_FAMILY_LAGRANGE,
									2,
									WMESH_NODES_FAMILY_LAGRANGE,
									3);
  std::cout << *g << std::endl;
#endif
  
#if 0
  wmesh_mat_t<double> y;
  wmesh_mat_t<double>::alloc(&y,1,4);
  wmesh_mat_t<double> e;
  wmesh_mat_t<double>::alloc(&e,1,g->n);
  y.v[0]=0;
  y.v[1]=1;
  y.v[2]=0;
  y.v[3]=0;

  wmesh_mat_gemm(1.0,y,*g,0.0,e);
  std::cout << e << std::endl;
#endif  
  return WMESH_STATUS_SUCCESS;
}
