#pragma once
#include "wmesh_t.hpp"
#include "VAR.hpp"
#include <iostream>
#include <math.h>
struct DG_VAR : public VAR<DG_VAR>
{

private:
  wmesh_mat_t<double>		m_values3;
  wmesh_shape_t 	m_shape;

  double * __restrict__	m_memory;
  bool 			m_owner;
  wmesh_t * 		m_mesh;

public:

  inline const wmesh_t* 	mesh()		const;
  inline const wmesh_mat_t<double>& 	matrix() 	const;
  inline const wmesh_shape_t* 		shape() 	const;
  inline wmesh_int_t 			ndofselm() 	const;
  inline wmesh_int_t 			ndofs() 	const;
  
  inline wmesh_t* 		mesh();
  inline wmesh_mat_t<double>& 	matrix();
  inline void 			clear();

  inline DG_VAR& operator *= 	(double alpha_);
  inline DG_VAR& operator += 	(const DG_VAR&that);
  inline DG_VAR& operator -= 	(const DG_VAR&that);
  //  inline DG_VAR& operator = 	(const temp_jacvar&that);;

  inline DG_VAR(wmesh_t * 	mesh_,
		wmesh_shape_t* 		shape);

  inline DG_VAR(wmesh_t * 	mesh_,
		wmesh_shape_t* 		shape,
		double*__restrict__ 		mem_,
		double*__restrict__ 		memory_,
		wmesh_int_t 		ld);
  
  inline ~DG_VAR();
  
  inline void 		dofselm		(wmesh_int_t id, wmesh_int_t icomp, double*__restrict__ dofs, wmesh_int_t inc) const;
  inline wmesh_int_t 	compute_nrms	(double nrms[]) const;
  inline void 		setdofselm	(wmesh_int_t id, wmesh_int_t icomp, const double*__restrict__ dofs, wmesh_int_t inc);
  
};

inline const wmesh_t* 		DG_VAR::mesh()	 const 	{ return this->m_mesh; };
inline wmesh_t* 		DG_VAR::mesh() 		{ return this->m_mesh; };
inline const wmesh_shape_t*	DG_VAR::shape()  const 	{ return &this->m_shape; };
inline wmesh_mat_t<double>& 		DG_VAR::matrix() 	{ return this->m_values3; };
inline const wmesh_mat_t<double>& 	DG_VAR::matrix() const 	{ return this->m_values3; };

inline DG_VAR& DG_VAR::operator *= (double alpha_)
{
  const wmesh_int_t m = this->m_values3.m;
  const wmesh_int_t n = this->m_values3.n;
  const wmesh_int_t ld = this->m_values3.ld;
  double * v = this->m_values3.v;
  for (wmesh_int_t j=0;j<n;++j)
    {
      for (wmesh_int_t i=0;i<m;++i)
	{
	  v[ld*j+i] *= alpha_;
	}      
    }
  return *this;
};

inline void DG_VAR::clear()
{
  const wmesh_int_t m = this->m_values3.m;
  const wmesh_int_t n = this->m_values3.n;
  const wmesh_int_t ld = this->m_values3.ld;
    double * v = this->m_values3.v;
  for (wmesh_int_t j=0;j<n;++j)
    {
      for (wmesh_int_t i=0;i<m;++i)
	{
	  v[ld*j+i] *= static_cast<double>(0);
	}      
    }

};

inline DG_VAR& DG_VAR::operator += (const DG_VAR&that)
{

  const wmesh_int_t m = this->m_values3.m;
  const wmesh_int_t n = this->m_values3.n;
  const wmesh_int_t ld = this->m_values3.ld;
  double * v = this->m_values3.v;
  if (m==that.m_values3.m)
    {
      fprintf(stderr,"error dimension (LINE %d FILE %s)\n",__LINE__,__FILE__);
      exit(1);
    }
  if (n==that.m_values3.n)
    {
      fprintf(stderr,"error dimension (LINE %d FILE %s)\n",__LINE__,__FILE__);
      exit(1);
    }

  for (wmesh_int_t j=0;j<n;++j)
    {
      for (wmesh_int_t i=0;i<m;++i)
	{
	  v[ld*j+i] += that.m_values3.v[that.m_values3.ld*j+i];
	}      
    }
  return *this;
};


inline DG_VAR& DG_VAR::operator -= (const DG_VAR&that)
{
  const wmesh_int_t m = this->m_values3.m;
  const wmesh_int_t n = this->m_values3.n;
  const wmesh_int_t ld = this->m_values3.ld;
    double * v = this->m_values3.v;
  if (m==that.m_values3.m)
    {
      fprintf(stderr,"error dimension (LINE %d FILE %s)\n",__LINE__,__FILE__);
      exit(1);
    }
  if (n==that.m_values3.n)
    {
      fprintf(stderr,"error dimension (LINE %d FILE %s)\n",__LINE__,__FILE__);
      exit(1);
    }
  for (wmesh_int_t j=0;j<n;++j)
    {
      for (wmesh_int_t i=0;i<m;++i)
	{
	  v[ld*j+i] -= that.m_values3.v[that.m_values3.ld*j+i];
	}      
    }
  return *this;
};


inline wmesh_int_t DG_VAR::ndofselm() 	const 	{ return m_shape.m_ndofs; }
inline wmesh_int_t DG_VAR::ndofs() 	const	{ return this->ndofselm() * m_mesh->m_num_triangles; }

inline DG_VAR::DG_VAR(wmesh_t * mesh_, wmesh_shape_t* shape) : m_mesh(mesh_)
{
  memcpy(&this->m_shape,shape,sizeof(wmesh_shape_t));
  this->m_memory = (double*__restrict__)calloc(m_mesh->m_num_triangles*shape->m_ndofs,sizeof(double));
  this->m_owner = true;
  wmesh_mat_t<double>::define(&this->m_values3,
			      this->m_shape.m_ndofs,
			      m_mesh->m_num_triangles,
			      this->m_memory,
			      this->m_shape.m_ndofs);
};
  
inline DG_VAR::DG_VAR(wmesh_t * mesh_,wmesh_shape_t* shape,double*__restrict__ mem_,double*__restrict__ memory_,wmesh_int_t ld) : m_mesh(mesh_)
{
  if (ld < shape->m_ndofs)
    {
      std::cerr << "Invalid ld(= " << ld << ") < nshape (= " << shape->m_ndofs << " )" <<  std::endl;
      exit(1);
    }
    
  memcpy(&this->m_shape,shape,sizeof(wmesh_shape_t));
  this->m_owner = false;
  this->m_memory = memory_;
  wmesh_mat_t<double>::define(&this->m_values3,this->m_shape.m_ndofs,this->m_mesh->m_num_triangles,this->m_memory,ld);
};

inline DG_VAR::~DG_VAR()
{
  if (this->m_owner && this->m_memory)
    {
      free(this->m_memory);
      this->m_memory = NULL;
    }
};

inline void DG_VAR::dofselm(wmesh_int_t id, wmesh_int_t icomp, double*__restrict__ dofs, wmesh_int_t inc) const
{
  wmesh_int_t n = ndofselm();
  for (wmesh_int_t i=0;i<n;++i)
    {
      dofs[i*inc] = m_values3.v[m_values3.ld * id + i];
    }
};
  
inline void DG_VAR::setdofselm(wmesh_int_t id, wmesh_int_t icomp, const double*__restrict__ dofs, wmesh_int_t inc) 
{
  wmesh_int_t n = ndofselm();
  for (wmesh_int_t i=0;i<n;++i)
    {
      m_values3.v[m_values3.ld * id + i] = dofs[i*inc];
    }
};
  
inline wmesh_int_t DG_VAR::compute_nrms(double nrms[]) const
{

  wmesh_int_t degree = this->m_shape.m_degree;
  //  I n1 = 1;
  for (wmesh_int_t ideg=0;ideg<=degree;++ideg)
    {
      double nrm=0.0;
      for (wmesh_int_t ielm=0;ielm<m_mesh->m_num_triangles;++ielm)
	{
	  double h1 = 0.0;
	  wmesh_int_t start = ((ideg) * (ideg+1) )/2;
	  wmesh_int_t bound = ((ideg+1) * (ideg+2) )/2;
	  for (wmesh_int_t j = start;j<bound;++j)
	    {
	      double x = m_values3.v[ielm*m_values3.ld+j];
	      h1 += x*x;
	    }

	  fprintf(stdout,"come and compute jacelm");
	  h1 *= 1.0;
	  // h1 *= m_mesh->jacelm[ielm];
	  exit(1);
	  
	  nrm+=h1;
	}
      nrms[1+ideg] = nrm;
    }

  nrms[0] = 0.0;
  for (wmesh_int_t ideg=0;ideg<=degree;++ideg)
    {
      nrms[0] += nrms[1+ideg];
    }

  nrms[0] = sqrt(nrms[0]);
  for (wmesh_int_t ideg=0;ideg<=degree;++ideg)
    {
      nrms[1+ideg] = sqrt(nrms[1+ideg]);
    }
    
  return 2+degree;
};
  
  

