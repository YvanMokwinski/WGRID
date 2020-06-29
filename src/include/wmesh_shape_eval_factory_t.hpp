#pragma once

#include "wmesh-types.hpp"
#include <map>
#include "wmesh_shape_eval_t.hpp"
#include "wmesh_nodes_t.hpp"
#include "wmesh_cubature_t.hpp"
#include <iostream>
template<typename T>
class wmesh_shape_eval_factory_t
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
  static wmesh_shape_eval_factory_t& instance()
  {
    static wmesh_shape_eval_factory_t s_instance;
    return s_instance;
  }
  
  using matdescr_t = typename std::pair<wmesh_int_t,const wmesh_mat_t<T>*>;
  using key_t = typename std::pair<wmesh_int_t,matdescr_t>;
  using val_t = wmesh_shape_eval_t<T>*;

  static const wmesh_shape_eval_t<T>* shape_eval_instance(wmesh_int_t			element_,
							  wmesh_int_t			shape_family_,
							  wmesh_int_t			shape_degree_,
							  wmesh_int_t			nodes_storage_,
							  const wmesh_mat_t<T> * 	nodes_)
  {
    wmesh_int_t shape_eval_uniqueid;
    uniqueid_encod(shape_eval_uniqueid,
		   element_,
		   shape_family_,
		   shape_degree_);
    
    auto ret = instance().s_map.find(key_t(shape_eval_uniqueid, matdescr_t(nodes_storage_,nodes_)));
    wmesh_shape_eval_t<T>*shape_eval = nullptr;
    if (ret == instance().s_map.end())
      {
	shape_eval = (wmesh_shape_eval_t<T>*)malloc(sizeof(wmesh_shape_eval_t<T>));
	wmesh_status_t status = wmesh_shape_eval_def(shape_eval,
						     element_,
						     shape_family_,
						     shape_degree_,
						     nodes_storage_,
						     nodes_);
	
	if (status != WMESH_STATUS_SUCCESS)
	  {
	    std::cerr << "wmesh_shape_eval_def error" << std::endl;
	    exit(1);
	  }
	auto ret_insert = instance().s_map.insert(std::pair<key_t,wmesh_shape_eval_t<T>*>( key_t(shape_eval_uniqueid, matdescr_t(nodes_storage_, nodes_)), shape_eval));
	if (ret_insert.second == false)
	  {
	    std::cerr << "not found but already registered" << std::endl;
	    exit(1);
	  }
      }
    else
      {
	std::cerr << "already registered" << std::endl;
	shape_eval = ret->second;
      }
    return shape_eval;
  }

  
  static const wmesh_shape_eval_t<T>* shape_eval_instance(wmesh_int_t element_,
							  const wmesh_shape_info_t&	shape_info_,
							  const wmesh_nodes_t<T> * 	nodes_)
  {
    
    return shape_eval_instance(element_,shape_info_.m_family,shape_info_.m_degree,nodes_->m_c_storage, &nodes_->m_c);
  }
  
  static const wmesh_shape_eval_t<T>* shape_eval_instance(wmesh_int_t element_,
							  const wmesh_shape_info_t&	shape_info_,
							  const wmesh_cubature_t<T> * 	cubature_)
  {
    
    return shape_eval_instance(element_,shape_info_.m_family,shape_info_.m_degree,cubature_->m_c_storage, &cubature_->m_c);
  }


  static const wmesh_shape_eval_t<T>* shape_eval_instance(wmesh_shape_t& shape_,
							  const wmesh_nodes_t<T> * 	nodes_)
  {
    return shape_eval_instance(shape_.m_element,shape_.m_family,shape_.m_degree,nodes_->m_c_storage, &nodes_->m_c);
  }
  
  static const wmesh_shape_eval_t<T>* shape_eval_instance(wmesh_shape_t&shape_,
							  const wmesh_cubature_t<T> * 	cubature_)
  {
    return shape_eval_instance(shape_.m_element,shape_.m_family,shape_.m_degree,cubature_->m_c_storage, &cubature_->m_c);
  }

  
private:
  wmesh_shape_eval_factory_t(){};
  ~wmesh_shape_eval_factory_t(){};
  wmesh_shape_eval_factory_t(const wmesh_shape_eval_factory_t&)= delete;
  wmesh_shape_eval_factory_t& operator=(const wmesh_shape_eval_factory_t&)= delete;
  
private:
  static std::map<key_t,val_t> s_map;
};

template<typename T>
std::map<typename wmesh_shape_eval_factory_t<T>::key_t,typename wmesh_shape_eval_factory_t<T>::val_t> wmesh_shape_eval_factory_t<T>::s_map;
