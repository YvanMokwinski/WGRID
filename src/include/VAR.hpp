#pragma once
#include "wmesh-types.hpp"
template <typename impl_t>
struct VAR
{  

  inline wmesh_int_t 		ndofselm() const
  {
    return static_cast<const impl_t&>(*this).ndofselm();
  };
  
  inline wmesh_int_t 		ndofs() const
  {
    return static_cast<const impl_t&>(*this).ndofs();
  };
  
  inline void 		dofselm(wmesh_int_t 	id,
				wmesh_int_t 	icomp,
				double * __restrict__ 	dofs,
				wmesh_int_t 	inc) 	const
  {
    return static_cast<const impl_t&>(*this).dofselm(id,icomp,dofs,inc);
  };
  
  inline const wmesh_shape_t* 	shape() const
  {
    return static_cast<const impl_t&>(*this).shape();
  };
  
  inline void 		clear()
  {
    return static_cast<impl_t&>(*this).clear();
  };
  
  inline void 		setdofselm(wmesh_int_t 			id,
				   wmesh_int_t 			icomp,
				   const double * __restrict__	dofs,
				   wmesh_int_t 			inc)
  {
    return static_cast<impl_t&>(*this).setdofselm(id,icomp,dofs,inc);
  };
  
};
