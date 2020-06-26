
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




#include <map>




template<typename T>
class wmesh_cubature_factory_t
{
private:  

  static   void uniqueid_encod(wmesh_int_t& 	val_,
			       wmesh_int_t 	element_,
			       wmesh_int_t 	family_,
			       wmesh_int_t 	degree_)
  {    
    val_ = element_ + (family_ << 3) + (degree_ << 6);
  }
  
  static   void uniqueid_decod(wmesh_int_t 	val_,
			       wmesh_int_t& 	element_,
			       wmesh_int_t& 	family_,
			       wmesh_int_t& 	degree_)
  {
    const wmesh_int_t s_mask_element = ~( ~((wmesh_int_t)0) << 3);
    const wmesh_int_t s_mask_family = ~( ~((wmesh_int_t)0) << 3) << 3;
    const wmesh_int_t s_mask_degree = ( ~((wmesh_int_t)0) << 6);
    element_ = val_ & s_mask_element;
    family_ = (val_ & s_mask_family) >> 3;
    degree_ = (val_ & s_mask_degree) >> 6;
  }

public:
  static wmesh_cubature_factory_t& instance()
  {
    static wmesh_cubature_factory_t s_instance;
    return s_instance;
  }
  
  static const wmesh_cubature_t<T>* cubature_instance(wmesh_int_t cubature_element_,
						      wmesh_int_t cubature_family_,
						      wmesh_int_t cubature_degree_)
  {    
    wmesh_int_t cubature_uniqueid;
    uniqueid_encod(cubature_uniqueid,
		   cubature_element_,
		   cubature_family_,
		   cubature_degree_);

    auto ret = instance().s_map.find(cubature_uniqueid);
    wmesh_cubature_t<T>*cubature = nullptr;
    if (ret == instance().s_map.end())
      {
	cubature = (wmesh_cubature_t<T>*)malloc(sizeof(wmesh_cubature_t<T>));
	wmesh_status_t status = wmesh_cubature_def(cubature,
				    cubature_element_,
				    cubature_family_,
				    cubature_degree_);
	if (status != WMESH_STATUS_SUCCESS)
	  {
	    std::cerr << "wmesh_cubature_def error" << std::endl;
	    exit(1);
	  }
	auto ret_insert = instance().s_map.insert(std::pair<wmesh_int_t,wmesh_cubature_t<T>*>(cubature_uniqueid, cubature));
	if (ret_insert.second == false)
	  {
	    std::cerr << "not found but already registered" << std::endl;
	    exit(1);
	  }
      }
    else
      {
	std::cerr << "already registered" << std::endl;
	cubature = ret->second;
      }
    return cubature;
  }
  
private:
  wmesh_cubature_factory_t(){};
  ~wmesh_cubature_factory_t(){};
  wmesh_cubature_factory_t(const wmesh_cubature_factory_t&)= delete;
  wmesh_cubature_factory_t& operator=(const wmesh_cubature_factory_t&)= delete;
  
private:
  static std::map<wmesh_int_t,wmesh_cubature_t<T>*> s_map;
};

template<typename T>
std::map<wmesh_int_t,wmesh_cubature_t<T>*> wmesh_cubature_factory_t<T>::s_map;


int main(int argc, char ** argv)
{
  auto cubature2 = wmesh_cubature_factory_t<double>::cubature_instance(WMESH_ELEMENT_TRIANGLE,
								       WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE,
								       4);
  
  std::cout << cubature2->m_c << std::endl;

  auto cubature3 = wmesh_cubature_factory_t<double>::cubature_instance(WMESH_ELEMENT_TETRAHEDRON,
								       WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE,
								       17);
  
  exit(1);
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
