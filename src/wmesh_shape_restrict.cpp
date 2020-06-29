

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh.hpp"
#include <chrono>

#include "bms.h"
#include "wmesh-utils.hpp"
#include "wmesh-blas.hpp"
#include "bms_templates.hpp"
#include "bms.hpp"
#include "wmesh_shape_restrict_t.hpp"
#include "wmesh_shape_eval_t.hpp"
#include "wmesh_shape_eval_factory_t.hpp"

#include "wmesh_nodes_factory_t.hpp"
#include "wmesh_nodes_boundary_factory_t.hpp"

template<typename T>
static std::ostream& operator<<(std::ostream&out_,
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


template<typename T>
const wmesh_mat_t<T>* wmesh_calculate_shape_restrict(wmesh_int_t 					element_,
						     wmesh_int_t 					ifacet_,
						     wmesh_int_t 					irot_,
					       
						     wmesh_int_t 					source_shape_family_,
						     wmesh_int_t 					source_shape_degree_,
						     
						     wmesh_int_t 					target_nodes_family_,
						     wmesh_int_t 					target_nodes_degree_)
{

  
  const wmesh_nodes_boundary_t<T>* nodes_boundary = wmesh_nodes_boundary_factory_t<T>::nodes_boundary_instance(element_,
													       target_nodes_family_,
													       target_nodes_degree_);
  const wmesh_shape_eval_t<T> * shape_eval = wmesh_shape_eval_factory_t<T>::shape_eval_instance(element_,
												source_shape_family_,
												source_shape_degree_,
												nodes_boundary->m_facets_nodes_storage,
												&nodes_boundary->m_facets_nodes[ifacet_][irot_]);	    
  return &shape_eval->m_f;
}

template
const wmesh_mat_t<double>* wmesh_calculate_shape_restrict<double>(wmesh_int_t 					element_,
								  wmesh_int_t 					ifacet_,
								  wmesh_int_t 					irot_,
								  
								  wmesh_int_t 					source_shape_family_,
								  wmesh_int_t 					source_shape_degree_,
								  
								  wmesh_int_t 					target_nodes_family_,
								  wmesh_int_t 					target_nodes_degree_);

template
const wmesh_mat_t<float>* wmesh_calculate_shape_restrict<float>(wmesh_int_t 					element_,
								wmesh_int_t 					ifacet_,
								wmesh_int_t 					irot_,
								
								wmesh_int_t 					source_shape_family_,
								wmesh_int_t 					source_shape_degree_,
								
								wmesh_int_t 					target_nodes_family_,
								wmesh_int_t 					target_nodes_degree_);
