#include "wmesh_t.hpp"
#include "wmesh_utils.hpp"
#include <iostream>
#include <string.h>
#include "wmesh-blas.hpp"


template<typename T>
wmesh_status_t wmesh_mat_gemm(wmesh_int_t 		storage_a_,
			      wmesh_int_t 		storage_b_,
			      wmesh_int_t 		storage_c_,
			      T 			alpha_,
			      const wmesh_mat_t<T>&	a_,
			      const wmesh_mat_t<T>&	b_,
			      T 			beta_,
			      wmesh_mat_t<T>&		c_)
{

  if (WMESH_STORAGE_INTERLEAVE == storage_c_)
    {
      wmesh_mat_gemm((WMESH_STORAGE_INTERLEAVE == storage_b_) ? "N" : "T",
		     (WMESH_STORAGE_INTERLEAVE == storage_a_) ? "N" : "T",
		     alpha_,
		     b_,
		     a_,
		     beta_,
		     c_);
    }
  else
    {
      wmesh_mat_gemm((WMESH_STORAGE_INTERLEAVE == storage_a_) ? "T" : "N",
		     (WMESH_STORAGE_INTERLEAVE == storage_b_) ? "T" : "N",
		     alpha_,
		     a_,
		     b_,
		     beta_,
		     c_);      
    }

  return WMESH_STATUS_SUCCESS;

}

template
wmesh_status_t wmesh_mat_gemm<float>(wmesh_int_t 			storage_a_,
				     wmesh_int_t 			storage_b_,
				     wmesh_int_t 			storage_c_,
				     float 				alpha_,
				     const wmesh_mat_t<float>&		a_,
				     const wmesh_mat_t<float>&		b_,
				     float 				beta_,
				     wmesh_mat_t<float>&		c_);


template
wmesh_status_t wmesh_mat_gemm<double>(wmesh_int_t 			storage_a_,
				     wmesh_int_t 			storage_b_,
				     wmesh_int_t 			storage_c_,
				     double 				alpha_,
				     const wmesh_mat_t<double>&		a_,
				     const wmesh_mat_t<double>&		b_,
				     double 				beta_,
				     wmesh_mat_t<double>&		c_);



wmesh_status_t wmesh_get_facet_dofs_ids(const wmesh_t * 	self_,
					wmesh_int_t 		element_type_,
					wmesh_int_t 		element_idx_,
					wmesh_int_t 		facet_type_,
					wmesh_int_t 		facet_idx_,
					wmesh_int_p 		ndofs_,
					wmesh_int_p 		dofs_,
					wmesh_int_t 		dofs_inc_);

#if 0
wmesh_status_t wmesh_get_dofs_ids(const wmesh_t * 	self_,
				  wmesh_int_t 		element_type_,
				  wmesh_int_t 		element_idx_,
				  wmesh_int_p 		dofs_,
				  wmesh_int_t 		dofs_inc_)
{
  for (wmesh_int_t i=0;i<self_->m_c2n.m_m[element_type_];++i)
    {
      dofs_[dofs_inc_*i] = self_->m_c2n.m_data[self_->m_c2n.m_ptr[element_type_] + self_->m_c2n.m_ld[element_type_] * element_idx_ + i];
    }
  return WMESH_STATUS_SUCCESS;
}
#endif


wmesh_status_t wmesh_get_facet_dofs_ids(const wmesh_t * 	self_,
					wmesh_int_t 		element_type_,
					wmesh_int_t 		element_idx_,
					wmesh_int_t 		facet_type_,
					wmesh_int_t 		facet_idx_,
					wmesh_int_p 		ndofs_,
					wmesh_int_p 		dofs_,
					wmesh_int_t 		dofs_inc_)
{
  const wmesh_int_t topodim = self_->m_topology_dimension;
  if (topodim==2)
    {
      const wmesh_int_t n = self_->m_s_e2n.m_m[element_type_];
      ndofs_[0] = n;
      for (wmesh_int_t i=0;i<n;++i)
	{
	  wmesh_int_t l = self_->m_s_e2n.m_data[self_->m_s_e2n.m_ptr[element_type_] + self_->m_s_e2n.m_ld[element_type_] * facet_idx_ + i];
	  dofs_[dofs_inc_*i] = self_->m_c2n.m_data[self_->m_c2n.m_ptr[element_type_] + self_->m_c2n.m_ld[element_type_] * element_idx_ + l];
	}
    }
  else   if (topodim==3)
    {
      if (facet_type_ == 0)
	{
	  const wmesh_int_t n = self_->m_s_t2n.m_m[element_type_];
	  ndofs_[0] = n;
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      wmesh_int_t l = self_->m_s_t2n.m_data[self_->m_s_t2n.m_ptr[element_type_] + self_->m_s_t2n.m_ld[element_type_] * facet_idx_ + i];
	      dofs_[dofs_inc_*i] = self_->m_c2n.m_data[self_->m_c2n.m_ptr[element_type_] + self_->m_c2n.m_ld[element_type_] * element_idx_ + l];
	    }
	}
      else if  (facet_type_ == 1)
	{
	  const wmesh_int_t n = self_->m_s_q2n.m_m[element_type_];
	  ndofs_[0] = n;
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      wmesh_int_t l = self_->m_s_q2n.m_data[self_->m_s_q2n.m_ptr[element_type_] + self_->m_s_q2n.m_ld[element_type_] * facet_idx_ + i];
	      dofs_[dofs_inc_*i] = self_->m_c2n.m_data[self_->m_c2n.m_ptr[element_type_] + self_->m_c2n.m_ld[element_type_] * element_idx_ + l];
	    }
	}
    }
  else
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
    }
  return WMESH_STATUS_SUCCESS;
}



template <typename T>
wmesh_status_t wmesh_get_cooelm(const wmesh_t * __restrict__	self_,
				wmesh_int_t			itype_,
				wmesh_int_t			ielm_,
				wmesh_int_t			cooelm_storage_,
				wmesh_int_t			cooelm_m_,
				wmesh_int_t			cooelm_n_,
				T * __restrict__			cooelm_,
				wmesh_int_t			cooelm_ld_)
{
  WMESH_CHECK_POINTER(self_);
  WMESH_CHECK_POINTER(cooelm_);
  switch(cooelm_storage_)
    {
    case WMESH_STORAGE_INTERLEAVE:
      {
	for (wmesh_int_t j=0;j<cooelm_n_;++j)
	  {
	    for (wmesh_int_t i=0;i<cooelm_m_;++i)
	      {
		cooelm_[cooelm_ld_*j+i]
		  = self_->m_coo[self_->m_coo_ld * ( self_->m_c2n.m_data[self_->m_c2n.m_ptr[itype_] + self_->m_c2n.m_ld[itype_] * ielm_ + j] - 1) + i];
	      }
	  }
	return WMESH_STATUS_SUCCESS;  
      }

    case WMESH_STORAGE_BLOCK:
      {
	for (wmesh_int_t j=0;j<cooelm_n_;++j)
	  {
	    for (wmesh_int_t i=0;i<cooelm_m_;++i)
	      {
		cooelm_[cooelm_ld_*i+j]
		  = self_->m_coo[self_->m_coo_ld * ( self_->m_c2n.m_data[self_->m_c2n.m_ptr[itype_] + self_->m_c2n.m_ld[itype_] * ielm_ + j] - 1) + i];
	      }
	  }
	return WMESH_STATUS_SUCCESS;  
      }

    }
  return WMESH_STATUS_INVALID_ENUM;  
}

template
wmesh_status_t wmesh_get_cooelm<float>(const wmesh_t * __restrict__	self_,
				wmesh_int_t			itype_,
				wmesh_int_t			ielm_,
				wmesh_int_t			cooelm_storage_,
				wmesh_int_t			cooelm_m_,
				wmesh_int_t			cooelm_n_,
				float * __restrict__			cooelm_,
				wmesh_int_t			cooelm_ld_);
template
wmesh_status_t wmesh_get_cooelm<double>(const wmesh_t * __restrict__	self_,
				wmesh_int_t			itype_,
				wmesh_int_t			ielm_,
				wmesh_int_t			cooelm_storage_,
				wmesh_int_t			cooelm_m_,
				wmesh_int_t			cooelm_n_,
				double * __restrict__			cooelm_,
				wmesh_int_t			cooelm_ld_);

extern "C"
{
 
  wmesh_status_t 	wmesh_get_topodim(const wmesh_t* 	self_,
					  wmesh_int_p 		topodim_)
  {
    topodim_[0] = self_->m_topology_dimension;
    return WMESH_STATUS_SUCCESS;  
  }
  
  wmesh_status_t wmesh_kill(wmesh_t* self_)
  {
    if (self_)
      {
	free(self_);
      }
    return WMESH_STATUS_SUCCESS;  
  };

    //!
  //! @brief Get the number of entities with a specific dimension.
  //!
  wmesh_status_t wmesh_factory	(wmesh_t** 		self__,
				 wmesh_int_t 		topology_dimension_,    

				 wmesh_int_t 		c2n_size_,
				 const_wmesh_int_p 	c2n_ptr_,
				 const_wmesh_int_p 	c2n_m_,
				 const_wmesh_int_p 	c2n_n_,
				 wmesh_int_p 		c2n_v_,
				 const_wmesh_int_p	c2n_ld_,

				 wmesh_int_t 		c_c_size_,
				 const_wmesh_int_p 	c_c_ptr_,
				 const_wmesh_int_p 	c_c_m_,
				 const_wmesh_int_p 	c_c_n_,
				 wmesh_int_p 		c_c_v_,
				 const_wmesh_int_p	c_c_ld_,
				 
				 wmesh_int_t 		bf2n_size_,
				 const_wmesh_int_p 	bf2n_ptr_,
				 const_wmesh_int_p 	bf2n_m_,
				 const_wmesh_int_p 	bf2n_n_,
				 wmesh_int_p 		bf2n_v_,
				 const_wmesh_int_p	bf2n_ld_,
				 
				 wmesh_int_t 		bf_c_size_,
				 const_wmesh_int_p 	bf_c_ptr_,
				 const_wmesh_int_p 	bf_c_m_,
				 const_wmesh_int_p 	bf_c_n_,
				 wmesh_int_p 		bf_c_v_,
				 const_wmesh_int_p	bf_c_ld_,
				 
				 wmesh_int_t		coo_m_,
				 wmesh_int_t		coo_n_,
				 double * 		coo_v_,
				 wmesh_int_t 		coo_ld_,
				 
				 wmesh_int_p 		n_c_v_,
				 wmesh_int_t		n_c_ld_)
  {
    
    self__[0] = (wmesh_t*)calloc(1,sizeof(wmesh_t));
    wmesh_t * self_ = self__[0];

    self_->m_topology_dimension = topology_dimension_;
    self_->m_num_nodes 	= coo_n_;
    self_->m_coo_n 	= coo_n_;
    self_->m_coo_m 	= coo_m_;
    self_->m_coo 	= coo_v_;
    self_->m_coo_ld 	= coo_ld_;
    
    wmesh_status_t status;

    
    status = wmesh_int_sparsemat_new(&self_->m_c2n,
				     c2n_size_,
				     c2n_ptr_,
				     c2n_m_,
				     c2n_n_,
				     c2n_v_,
				     c2n_ld_);
    WMESH_STATUS_CHECK(status);
    
    
    status = wmesh_int_sparsemat_new(&self_->m_c_c,
				     c_c_size_,
				     c_c_ptr_,
				     c_c_m_,
				     c_c_n_,
				     c_c_v_,
				     c_c_ld_);
    WMESH_STATUS_CHECK(status);
    
    
    wmesh_int_mat_def(&self_->m_n_c,
		      1,
		      coo_n_,
		      n_c_v_,
		      n_c_ld_);
    
    wmesh_int_t num_cells = 0;
    for (wmesh_int_t i=0;i<c2n_size_;++i) num_cells += c2n_n_[i];

    self_->m_num_cells  = num_cells;


    if (bf2n_size_>0)
      {
	wmesh_int_sparsemat_new(&self_->m_bf2n,
				c2n_size_,
				bf2n_ptr_,
				bf2n_m_,
				bf2n_n_,
				bf2n_v_,
				bf2n_ld_);

	
	wmesh_int_sparsemat_new(&self_->m_bf_c,
				c2n_size_,
				bf_c_ptr_,
				bf_c_m_,
				bf_c_n_,
				bf_c_v_,
				bf_c_ld_);
      }
    
    if (c2n_size_==1)
      {
	self_->m_num_edges 		= c2n_n_[0];
      }
    else if (c2n_size_==2)
      {
	self_->m_num_triangles 		= c2n_n_[0];
	self_->m_num_quadrilaterals 	= c2n_n_[1];
      }
    else if (c2n_size_==4)
      {
	self_->m_num_tetrahedra 	= c2n_n_[0];
	self_->m_num_pyramids 		= c2n_n_[1];
	self_->m_num_wedges 		= c2n_n_[2];
	self_->m_num_hexahedra 		= c2n_n_[3];
      }
    else
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
      }

    wmesh_int_t dimension = (c2n_size_==4) ? 3 : ( (c2n_size_==2) ? 2 : 1 );
    if (dimension > 1)
      {
	wmesh_status_t status = wmesh_build_s_e2n(&self_->m_s_e2n, dimension);
	WMESH_STATUS_CHECK(status);
	if (dimension > 2)
	  {
	    //
	    // 1.b/ Local triangles to nodes.
	    //
	    status  = wmesh_build_s_t2n(&self_->m_s_t2n);
	    WMESH_STATUS_CHECK(status);
	    
	    //
	    // 1.c/ Local quadrilaterals to nodes.
	    //
	    status = wmesh_build_s_q2n(&self_->m_s_q2n);
	    WMESH_STATUS_CHECK(status);	    
	  }
      }
    return WMESH_STATUS_SUCCESS;
  }

  //!
  //! @brief Get the number of entities with a specific dimension.
  //!
  wmesh_status_t wmesh_def	(wmesh_t** 		self__,
				 wmesh_int_t 		topology_dimension_, 
				 wmesh_int_t 		ntypes_,
				 const_wmesh_int_p 	c2n_ptr_,
				 const_wmesh_int_p 	c2n_m_,
				 const_wmesh_int_p 	c2n_n_,
				 wmesh_int_p 		c2n_v_,
				 const_wmesh_int_p	c2n_ld_,
				 wmesh_int_t		coo_m_,
				 wmesh_int_t		coo_n_,
				 double * 		coo_,
				 wmesh_int_t 		coo_ld_)
  {
    wmesh_status_t status;
    self__[0] = (wmesh_t*)calloc(1,sizeof(wmesh_t));
    wmesh_t * self_ = self__[0];
    self_->m_topology_dimension = topology_dimension_;
    self_->m_num_nodes	= coo_n_;
    self_->m_coo_n = coo_n_;
    self_->m_coo_m = coo_m_;
    self_->m_coo_ld 	= coo_ld_;    
    self_->m_coo 	= coo_;

    
    status = wmesh_int_sparsemat_new(&self_->m_c2n,
				     ntypes_,
				     c2n_ptr_,
				     c2n_m_,
				     c2n_n_,			    
				     c2n_v_,
				     c2n_ld_);
    WMESH_STATUS_CHECK(status);
    wmesh_int_t num_cells = 0;
    for (wmesh_int_t i=0;i<ntypes_;++i) num_cells += c2n_n_[i];
    
    status = wmesh_int_mat_def(&self_->m_n_c,
			       1,
			       self_->m_num_nodes,
			       (wmesh_int_p)calloc(self_->m_num_nodes,sizeof(wmesh_int_t)),
			       1);
    WMESH_STATUS_CHECK(status);
    

    wmesh_int_t c_c_m[4] = {1,1,1,1};
    wmesh_int_t c_c_ld[4] = {1,1,1,1};
    wmesh_int_t c_c_ptr[4+1];
    wmesh_int_t c_c_n[4];

    
    for (wmesh_int_t i=0;i<ntypes_;++i)
      {
	c_c_n[i] = c2n_n_[i];
      }
    
    c_c_ptr[0] = 0;
    for (wmesh_int_t i=0;i<ntypes_;++i)
      {
	c_c_ptr[i+1] = c_c_ptr[i] + c_c_ld[i] * c_c_n[i];
      }

    status = wmesh_int_sparsemat_new(&self_->m_c_c,
				     ntypes_,
				     c_c_ptr,
				     c_c_m,
				     c_c_n,
				     (wmesh_int_p)calloc(c_c_ptr[ntypes_],sizeof(wmesh_int_t)),
				     c_c_ld);
    

#if 0
    self_->m_flag_nodes = (wmesh_int_p)calloc(num_nodes_,sizeof(wmesh_int_t));
    self_->m_flag_cells = (wmesh_int_p)calloc(num_cells,sizeof(wmesh_int_t));
#endif
    self_->m_num_cells  = num_cells;

    //
    // 1/ Initialize the reference shape edges to nodes.
    //    
    //
    // 1.a/ Local edges to nodes.
    //
    wmesh_int_t dimension = (ntypes_==4) ? 3 : ( (ntypes_==2) ? 2 : 1 );
    if (dimension > 1)
      {
	wmesh_status_t status  = wmesh_build_s_e2n(&self_->m_s_e2n, dimension);
	WMESH_STATUS_CHECK(status);
	if (dimension > 2)
	  {
	    //
	    // 1.b/ Local triangles to nodes.
	    //
	    status  = wmesh_build_s_t2n(&self_->m_s_t2n);
	    WMESH_STATUS_CHECK(status);
	    
	    //
	    // 1.c/ Local quadrilaterals to nodes.
	    //
	    status = wmesh_build_s_q2n(&self_->m_s_q2n);
	    WMESH_STATUS_CHECK(status);	    
	  }
      }
    return WMESH_STATUS_SUCCESS;
  }

  wmesh_status_t wmesh_fprintf(const wmesh_t* 		self_,
				FILE * out_)
  {
    wmesh_int_t status = wmesh_int_sparsemat_fprintf(&self_->m_c2n,out_);
    WMESH_STATUS_CHECK(status);
    return WMESH_STATUS_SUCCESS;
  };

  
#if 0
  wmesh_status_t wmesh_reorder_cell(wmesh_t*self_,double * box)
  {
    double xyz[3];
    double x,y,z;
    const double * coo = self_->m_coo;
#if 0
    double xyz[3],box[6];

    double xmin,ymin,zmin;
    double xmax,ymax,zmax;
    
    xmin = xmax = coo[0];
    ymin = ymax = coo[1];
    zmin = zmax = coo[2];
    for (wmesh_int_t i=0;i<num_nodes;++i)
      {
	
	x = coo[3*i+0];
	y = coo[3*i+1];
	z = coo[3*i+2];

	xmin = (x < xmin) ? x : xmin;
	ymin = (y < ymin) ? y : ymin;
	zmin = (z < zmin) ? z : zmin;

	xmax = (x > xmax) ? x : xmax;
	ymax = (y > ymax) ? y : ymax;
	zmax = (z > zmax) ? z : zmax;
	
      }
    box[0] = xmin;
    box[1] = ymin;
    box[2] = zmin;

    box[3] = xmax;
    box[4] = ymax;
    box[5] = zmax;

#endif
    
    for (int l=0;l<4;++l)
      {
	auto c2n_m 	= self_->m_c2n.m_m[l];
	auto c2n_ld 	= self_->m_c2n.m_ld[l];
	auto c2n_n	= self_->m_c2n.m_n[l];
	auto c2n_v 	= self_->m_c2n.m_data + self_->m_c2n.m_ptr[l];
	
	unsigned long long int* p = (unsigned long long int*)malloc(2 * c2n_n * sizeof(unsigned long long int));
	for (wmesh_int_t j=0;j<c2n_n;++j)
	  {
	    x = y = z = 0.0;
	    for (wmesh_int_t i=0;i<c2n_m;++i)
	      {
		wmesh_int_t v = c2n_v[j*c2n_ld + i];
		x += coo[3*(v-1)+0];
		y += coo[3*(v-1)+1];
		z += coo[3*(v-1)+2];
	      }
	    xyz[0] = x/c2n_m;
	    xyz[1] = y/c2n_m;
	    xyz[2] = z/c2n_m;
	    p[2 * j + 0] = hilbert_coordinate(xyz,
					      box,
					      31);
	    p[2 * j + 1] = j;
	  }
	
	
	qsort(p,
	      c2n_n,
	      sizeof(unsigned long long int)*2,
	      hilbert_sort_predicate);
#if 0
	for (wmesh_int_t i=0;i<num_nodes;++i)
	  {
	    // 0     // 0
	    // 4     // 5
	    // 3     // 4
	    // 5     // 2
	    // 2     // 1
	    // 1     // 3
	    p[2* ( p[2*i+1] ) + 0] = i;	
	  }
#endif
	//
	// Now reorder the coordinates.
	//
	double * tmp = (double*)malloc(sizeof(double)*c2n_m*self_->m_c2n.m_n[l]);
	for (wmesh_int_t j=0;j<self_->m_c2n.m_n[l];++j)
	  {
	    wmesh_int_t to = p[2*j+1];
	    for (wmesh_int_t i=0;i<c2n_m;++i)
	      {		
		tmp[j*c2n_m + i] = c2n_v[to*c2n_ld + i];
	      }
	  }
	
	for (wmesh_int_t j=0;j<self_->m_c2n.m_n[l];++j)
	  {
	    for (wmesh_int_t i=0;i<c2n_m;++i)
	      {		
		c2n_v[j*c2n_ld + i] = tmp[j*c2n_m + i];
	      }
	  }

	free(tmp);	
	free(p);
      }
    return 0;
  }
#endif

  static int hilbert_sort_predicate(const void * a_,
				    const void * b_)
  {
    const unsigned long long int * a = (const unsigned long long int*)a_;
    const unsigned long long int * b = (const unsigned long long int*)b_;
    if (*a < *b)
      {
	return -1;
      }
    else if (*a > *b)    
      {
	return 1;
      }
    else
      {
	return 0;
      }
  }

#if 0
  
struct hilbert_reordering_t
{
  unsigned long long int cod;
  unsigned long long int id;
};


int hilbert_comp_t(const void * a_,const void * b_)
{
  const I * a = (const I*)a_;
  const I * b = (const I*)b_;
  if (*a < *b)
    {
      return -1;
    }
  else if (*a > *b)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

#if 0
  coordinates2d(mesh_->nelm, 
    bb,
    vertices,
    2,
    (unsigned long long int*)tmp,
    2);

  qsort(tmp,mesh_->nelm,2*sizeof(I),hilbert_comp_t);
  for (I ielm=0;ielm<mesh_->nelm;++ielm)
    {
      // perm_[tmp[2*ielm+1]] = ielm;
      perm_[ielm] = tmp[2*ielm+1];
    }
#endif

static void coordinates2d(long long int				N_, 
			  double				box_[],
			  double* 				vertex_,
			  long long int                         vertexoff_,
			  unsigned long long int *  		cod2_,
			  long long int				codoff_)
{
  static const unsigned long long m = 1LL<<62;
  static const int BitTab[2] 	= {1,2};
  static const int GeoCod[4]	= {1,2,0,3};
  static const int HilCod[4][4] = {{0,3,2,1}, {0,1,2,3}, {0,1,2,3}, {2,1,0,3}};
  static const double len 	= 4.611686018427387904e+18;
  double box[4],dbl;
  unsigned long long int IntCrd[2],cod;
  int rot[4];
  box[0] 	= box_[0];
  box[1] 	= box_[1];
  box[2] 	= len / (box_[2] - box_[0]);
  box[3] 	= len / (box_[3] - box_[1]);
  double loc[2];
  { int i;
    for(i=0; i<N_; i++)
      {
	loc[0] = vertex_[vertexoff_*i+0];
	loc[1] = vertex_[vertexoff_*i+1];
	/* Convert double precision coordinates to integers */
	dbl 	= (loc[0] - box[0]) * box[0+2];
	IntCrd[0] = (unsigned long long int)dbl;
	dbl 	= (loc[1] - box[1]) * box[1+2];
	IntCrd[1] = (unsigned long long int)dbl;
	/* Binary hilbert renumbering loop */
	cod = 0;
	rot[0] = GeoCod[0];
	rot[1] = GeoCod[1];
	rot[2] = GeoCod[2];
	rot[3] = GeoCod[3];
	{ int b;
	  for(b=0;b<31;b++)
	    {
	      int GeoWrd = 0;

	      if(IntCrd[0] & m)
		GeoWrd |= BitTab[0];
	      IntCrd[0] = IntCrd[0]<<1;

	      if(IntCrd[1] & m)
		GeoWrd |= BitTab[1];
	      IntCrd[1] = IntCrd[1]<<1;

	      const int NewWrd = rot[ GeoWrd ];
	      cod = cod<<2 | NewWrd;
	      rot[0] = HilCod[ NewWrd ][ rot[0] ];
	      rot[1] = HilCod[ NewWrd ][ rot[1] ];
	      rot[2] = HilCod[ NewWrd ][ rot[2] ];
	      rot[3] = HilCod[ NewWrd ][ rot[3] ];
	    } }
	cod2_[codoff_*i+0] 	= cod;
	cod2_[codoff_*i+1] 	= i;
      } }


  
}
#endif

  
  wmesh_status_t wmesh_reorder(wmesh_t*self_)
  {
#if 1
    wmesh_int_t num_nodes 	= self_->m_num_nodes;
    double * 	coo 	  	= self_->m_coo;
    unsigned long long int* p	= (unsigned long long int*)malloc(2 * num_nodes * sizeof(unsigned long long int));
    double xyz[3],box[6];

    double x,y,z;
    double xmin,ymin,zmin;
    double xmax,ymax,zmax;
    
    xmin = xmax = coo[0];
    ymin = ymax = coo[1];
    zmin = zmax = coo[2];
    for (wmesh_int_t i=0;i<num_nodes;++i)
      {
	
	x = coo[3*i+0];
	y = coo[3*i+1];
	z = coo[3*i+2];

	xmin = (x < xmin) ? x : xmin;
	ymin = (y < ymin) ? y : ymin;
	zmin = (z < zmin) ? z : zmin;

	xmax = (x > xmax) ? x : xmax;
	ymax = (y > ymax) ? y : ymax;
	zmax = (z > zmax) ? z : zmax;
	
      }
    box[0] = xmin;
    box[1] = ymin;
    box[2] = zmin;

    box[3] = xmax;
    box[4] = ymax;
    box[5] = zmax;
    printf("hey\n");

    for (wmesh_int_t i=0;i<num_nodes;++i)
      {
	xyz[0] = coo[3*i+0];
	xyz[1] = coo[3*i+1];
	xyz[2] = coo[3*i+2];
	
	p[2 * i + 0] = hilbert_coordinate(xyz,
					  box,
					  37);
	p[2 * i + 1] = i;
      }

    qsort(p,
	  num_nodes,
	  sizeof(unsigned long long int)*2,
	  hilbert_sort_predicate);
    
    for (wmesh_int_t i=0;i<num_nodes;++i)
      {
	// 0     // 0
	// 4     // 5
	// 3     // 4
	// 5     // 2
	// 2     // 1
	// 1     // 3
	p[2* ( p[2*i+1] ) + 0] = i;
	
      }

    for (int l=0;l<self_->m_c2n.m_size;++l)
      for (int j=0;j<self_->m_c2n.m_n[l];++j)
	for (int i=0;i<self_->m_c2n.m_m[l];++i)
	  {
	    wmesh_int_t v = self_->m_c2n.m_data[ self_->m_c2n.m_ptr[l] + j*self_->m_c2n.m_ld[l] + i];
	    self_->m_c2n.m_data[ self_->m_c2n.m_ptr[l] + j*self_->m_c2n.m_ld[l] + i] = p[2*(v-1) + 0] + 1;
	  }
    
    //
    // Now reorder the coordinates.
    //
    double * tmp = (double*)malloc(sizeof(double)*3*num_nodes);
    for (wmesh_int_t i=0;i<num_nodes;++i)
      {
	wmesh_int_t j = p[2*i+1];
	tmp[3*i+0] = coo[3*j+0];
	tmp[3*i+1] = coo[3*j+1];
	tmp[3*i+2] = coo[3*j+2];
      }
    
    for (wmesh_int_t i=0;i<num_nodes;++i)
      {
	coo[3*i+0] = tmp[3*i+0];
	coo[3*i+1] = tmp[3*i+1];
	coo[3*i+2] = tmp[3*i+2];
      }    
    free(tmp);
    //
    // Now do the renumbering of c2n .
    //
    free(p);
#endif
    return 0;
  }
  
#if 0
  
  wmesh_status_t wmesh_reorder_nodes(wmesh_t*		self_,
				     wmesh_int_n 	p_n_,
				     const_wmesh_int_p 	p_v_,
				     wmesh_int_n 	p_ld_)
  {
    wmesh_status_t status;
    wmesh_int_t num_nodes 	= self_->m_num_nodes;
    wmesh_int_t work_n 		= num_nodes * 3;
    void * work 		= malloc(sizeof(double) * num_nodes * 3);
    if (!work)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }
    status = wmesh_permutate(self_->m_coo,p_n_,p_v_,p_ld_,work_n,(double*)work);
    WMESH_STATUS_CHECK(status);
    status = wmesh_permutate(self_->m_flag_nodes,p_n_,p_v_,p_ld_,work_n,(wmesh_int_p)work);
    WMESH_STATUS_CHECK(status);
    free(p);
    return 0;
  }

  wmesh_status_t wmesh_reorder_cells(wmesh_t*		self_,
				     wmesh_int_n 	p_n_,
				     const_wmesh_int_p 	p_v_,
				     wmesh_int_n 	p_ld_)
  {
    wmesh_status_t status;
    wmesh_int_t num_nodes 	= self_->m_num_nodes;
    wmesh_int_t work_n 		= num_nodes * 3;
    void * work 		= malloc(sizeof(double) * num_nodes * 3);
    if (!work)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }
    status = wmesh_permutate(self_->m_coo,p_n_,p_v_,p_ld_,work_n,(double*)work);
    WMESH_STATUS_CHECK(status);
    status = wmesh_permutate(self_->m_flag_nodes,p_n_,p_v_,p_ld_,work_n,(wmesh_int_p)work);
    WMESH_STATUS_CHECK(status);
    free(p);
    return 0;
  }
#endif

  
  wmesh_status_t wmesh_write(const wmesh_t* 		self_,
			      const char * 		filename_)
  {
    WMESH_CHECK_POINTER(self_);
    WMESH_CHECK_POINTER(filename_);
    const char * extension = file_extension(filename_);
    if (!strcmp(extension,".mesh")||!strcmp(extension,".meshb"))
      {
	return wmesh_write_medit(self_,0 == strcmp(extension,".meshb"),filename_);
      }
    else if (!strcmp(extension,".vtk"))
      {
	return wmesh_write_vtk(self_,filename_);
      }
    else
      {
	std::cerr << "// wmesh_write, no extension found in filename '" << filename_<< "'" << std::endl;
	return WMESH_STATUS_INVALID_ARGUMENT;
      }
    return WMESH_STATUS_INVALID_ARGUMENT;
  };

  wmesh_status_t wmesh_read(wmesh_t ** 		self__,
			    const char * 	filename_)
  {
    WMESH_CHECK_POINTER(self__);
    WMESH_CHECK_POINTER(filename_);
    const char * extension = file_extension(filename_);
    wmesh_int_t status;
    if (!strcmp(extension,".mesh")||!strcmp(extension,".meshb"))
      {
	return wmesh_read_medit(self__,
				0 == strcmp(extension,".meshb"),
				filename_);
      }
    else if (!strcmp(extension,".vtk"))
      {
	std::cerr << "// wmesh_read, cannot read from .vtk file. " << std::endl;
	WMESH_STATUS_CHECK(WMESH_STATUS_NOT_IMPLEMENTED);
      }
    else
      {
	std::cerr << "// wmesh_read, no extension found in filename '" << filename_ << "'" << std::endl;
	return WMESH_STATUS_INVALID_ARGUMENT;
      }
    return status;
  };
  

};
