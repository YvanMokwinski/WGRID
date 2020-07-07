#pragma once
#include "VAR.hpp"
#include "wmesh_t.hpp"

struct CG_VAR : public VAR<CG_VAR>
{
private: wmesh_mat_t<double> 	m_values3;
  double* 			m_memory3;
  
public:
  wmeshspace_t * 	m_meshspace;
  wmesh_int_t	m_ncomponents;
  wmesh_int_t	m_ndofsPerComponent;
  wmesh_shape_t	m_shape;
  bool 		m_owner;
  
  template <typename fct_t>
  void ih(fct_t f)
  {
    //
    // numnodes.
    //
    double uvw[3];
    wmesh_int_t nddlu;
    nddlu = m_meshspace->m_ndofs;
    wmesh_int_t nddlv	= nddlu;
    for (int i=0;i<nddlu;++i)
      {
	//
	fprintf(stdout," come here to correct, we need dof coordinates \n");
	exit(1);
	//
	double x=0,y=0;
#if 0
	double x= m_meshspace->coo[2*i+0];
	double y= m_meshspace->coo[2*i+1];
#endif
	f(x,y,uvw);
	m_memory3[i] = uvw[0];
	m_memory3[nddlu+i] = uvw[1];
      }
  };
  
  inline void	clear()
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
  
  inline const  wmesh_shape_t* shape() const
  {
    return &this->m_shape;
  };

  inline void dofselm(wmesh_int_t id, wmesh_int_t icomp, double * __restrict__ dofs, wmesh_int_t inc) const
  {
    fprintf(stdout," come here to correct, we need dof cnc \n");
    exit(1);
#if 0
    wmesh_int_t n = ndofselm();
    for (wmesh_int_t i=0;i<n;++i)
      {
	dofs[i*inc] = m_values3.v[m_values3.ld * icomp + m_meshspace->cnc[6*id+i]-1];
      }
#endif
  };
  
  inline void setdofselm(wmesh_int_t id, wmesh_int_t icomp, const double * __restrict__ dofs, wmesh_int_t inc) 
  {
    fprintf(stdout," come here to correct, we need dof cnc \n");
    exit(1);
#if 0
    wmesh_int_t n = ndofselm();
    for (wmesh_int_t i=0;i<n;++i)
      {
	m_values3.v[m_values3.ld * icomp + m_meshspace->cnc[6*id+i]-1] = dofs[i*inc];
      }
#endif
  };
  
  inline wmesh_int_t ndofselm() const 	{ return m_shape.m_ndofs; }
  inline wmesh_int_t ndofs() const	{ return this->m_ndofsPerComponent * this->m_ncomponents; }
  
  inline CG_VAR(wmeshspace_t * meshspace_, wmesh_int_t ncomponents, wmesh_shape_t* shape) : m_meshspace(meshspace_)
  {
    this->m_ndofsPerComponent = meshspace_->m_ndofs;
    this->m_ncomponents = ncomponents;
    memcpy(&this->m_shape,shape,sizeof(wmesh_shape_t));
    this->m_memory3 = (double * __restrict__)calloc(this->m_ndofsPerComponent * this->m_ncomponents,sizeof(double));
    this->m_owner = true;

    this->m_values3.v = this->m_memory3;
    this->m_values3.m = this->m_ndofsPerComponent;
    this->m_values3.n = this->m_ncomponents;
    this->m_values3.ld = this->m_ndofsPerComponent;

  };
  
  inline CG_VAR(wmeshspace_t * meshspace_,wmesh_int_t ncomponents, wmesh_shape_t* shape,double * __restrict__ mem_,double * __restrict__ memory_,wmesh_int_t ld) : m_meshspace(meshspace_)
  {
    this->m_ndofsPerComponent = meshspace_->m_ndofs;
    m_ncomponents = ncomponents;
    if (ld < this->m_ndofsPerComponent)
      {
	std::cerr << "Invalid ld(= " << ld << ") < ndofsPerComponent (= " << this->m_ndofsPerComponent << " )" <<  std::endl;
      }
    memcpy(&this->m_shape,shape,sizeof(wmesh_shape_t));
    this->m_owner = false;
    this->m_memory3 = memory_;
    this->m_values3.v = this->m_memory3;
    this->m_values3.m = this->m_ndofsPerComponent;
    this->m_values3.n = this->m_ncomponents;
    this->m_values3.ld = ld;

  };

  inline ~CG_VAR()
  {
    if (this->m_owner && this->m_memory3)
      {
	free(this->m_memory3);
	this->m_memory3 = NULL;
      }
  };
  
};
