#pragma once

#include "wmesh-types.hpp"
#include "wmesh_nodes_boundary_t.hpp"

#include <map>
#include <iostream>

template<typename T>
class wmesh_nodes_boundary_factory_t
{
public: 

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
  static wmesh_nodes_boundary_factory_t& instance()
  {
    static wmesh_nodes_boundary_factory_t s_instance;
    return s_instance;
  }
  
  static const wmesh_nodes_boundary_t<T>* nodes_boundary_instance(wmesh_int_t nodes_element_,
							 wmesh_int_t nodes_family_,
							 wmesh_int_t nodes_degree_)
  {    
    wmesh_int_t nodes_uniqueid;
    uniqueid_encod(nodes_uniqueid,
		   nodes_element_,
		   nodes_family_,
		   nodes_degree_);
    
    auto ret = instance().s_map.find(nodes_uniqueid);
    wmesh_nodes_boundary_t<T>*nodes = nullptr;
    if (ret == instance().s_map.end())
      {
	nodes = (wmesh_nodes_boundary_t<T>*)malloc(sizeof(wmesh_nodes_boundary_t<T>));
	wmesh_status_t status = wmesh_nodes_boundary_def(nodes,
							 nodes_element_,
							 nodes_family_,
							 nodes_degree_);
	if (status != WMESH_STATUS_SUCCESS)
	  {
	    std::cerr << "wmesh_nodes_def error" << std::endl;
	    exit(1);
	  }
	auto ret_insert = instance().s_map.insert(std::pair<wmesh_int_t,wmesh_nodes_boundary_t<T>*>(nodes_uniqueid, nodes));
	if (ret_insert.second == false)
	  {
	    std::cerr << "not found but already registered" << std::endl;
	    exit(1);
	  }
      }
    else
      {
	std::cerr << "already registered" << std::endl;
	nodes = ret->second;
      }
    return nodes;
  }

  
  static const wmesh_nodes_boundary_t<T>* nodes_instance(wmesh_int_t element_,
							 const wmesh_nodes_info_t& nodes_info_)
  {
    return nodes_instance(element_,
			  nodes_info_.m_family,
			  nodes_info_.m_degree);
  }
  
private:
  wmesh_nodes_boundary_factory_t(){};
  ~wmesh_nodes_boundary_factory_t(){};
  wmesh_nodes_boundary_factory_t(const wmesh_nodes_boundary_factory_t&)= delete;
  wmesh_nodes_boundary_factory_t& operator=(const wmesh_nodes_boundary_factory_t&)= delete;
  
private:
  static std::map<wmesh_int_t,wmesh_nodes_boundary_t<T>*> s_map;
};

template<typename T>
std::map<wmesh_int_t,wmesh_nodes_boundary_t<T>*> wmesh_nodes_boundary_factory_t<T>::s_map;
