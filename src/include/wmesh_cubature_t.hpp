#pragma once

#include "wmesh-types.hpp"

template<typename T>
struct wmesh_cubature_t
{
private:
  wmesh_int_t 		m_element;
  wmesh_int_t 		m_family;
  wmesh_int_t 		m_degree;
  wmesh_int_t 		m_c_storage;
  wmesh_mat_t<T> 	m_c;
  wmesh_mat_t<T> 	m_w;

public:
  wmesh_cubature_t(wmesh_int_t 		element_,
		   wmesh_int_t 		family_,
		   wmesh_int_t 		degree_);

  inline wmesh_int_t 		get_element() 			const {return this->m_element;};
  inline wmesh_int_t 		get_family() 			const {return this->m_family;};
  inline wmesh_int_t 		get_degree() 			const {return this->m_degree;};
  inline wmesh_int_t 		get_coordinates_storage() 	const {return this->m_c_storage;};
  inline const wmesh_mat_t<T>& 	get_coordinates() 		const {return this->m_c;};
  inline const wmesh_mat_t<T>& 	get_weights() 			const {return this->m_w;};
  inline wmesh_int_t 		get_num_points()		const {return this->m_w.n;};
  
};

#if 0
template<typename T>
wmesh_status_t wmesh_cubature_def(wmesh_cubature_t<T>*__restrict__ self_,
				  wmesh_int_t 		element_,
				  wmesh_int_t 		family_,
				  wmesh_int_t 		degree_);
#endif
