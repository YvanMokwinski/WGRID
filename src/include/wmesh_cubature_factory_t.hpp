#pragma once

#include "wmesh-types.hpp"
#include "wmesh_cubature_t.hpp"
#include <map>
#include <iostream>


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
