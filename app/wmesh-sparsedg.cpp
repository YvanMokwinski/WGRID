

#include "cmdline.hpp"
#include <iostream>
#include <math.h>
#include "wmesh-blas.hpp"
#include <algorithm>
#include "bms.h"
#include "wmesh.hpp"
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
wmesh_status_t wmesh_csrjacobian(wmesh_int_t 		c2n_size_,
				 const_wmesh_int_p 	c2n_ptr_,
				 const_wmesh_int_p 	c2n_m_,
				 const_wmesh_int_p 	c2n_n_,
				 const_wmesh_int_p 	c2n_v_,
				 const_wmesh_int_p 	c2n_ld_,
				 
				 wmesh_int_t 		c2c_size_,
				 const_wmesh_int_p 	c2c_ptr_,
				 const_wmesh_int_p 	c2c_m_,
				 const_wmesh_int_p 	c2c_n_,
				 const_wmesh_int_p 	c2c_v_,
				 const_wmesh_int_p 	c2c_ld_,
				 
				 const_wmesh_int_p 	row_size_blocks_,
				 const_wmesh_int_p 	col_size_blocks_,
				 
				 wmesh_int_p 		csr_m__,
				 wmesh_int_p 		csr_n__,
				 wmesh_int_p 		csr_nnz__,
				 wmesh_int_p*		csr_row_ptr__,
				 wmesh_int_p*		csr_col_ind__)
{
  
  WMESH_POINTER_CHECK(c2n_ptr_);
  WMESH_POINTER_CHECK(c2n_m_);
  WMESH_POINTER_CHECK(c2n_n_);
  WMESH_POINTER_CHECK(c2n_v_);
  WMESH_POINTER_CHECK(c2n_ld_);

  WMESH_POINTER_CHECK(c2c_ptr_);
  WMESH_POINTER_CHECK(c2c_m_);
  WMESH_POINTER_CHECK(c2c_n_);
  WMESH_POINTER_CHECK(c2c_v_);
  WMESH_POINTER_CHECK(c2c_ld_);

  WMESH_POINTER_CHECK(row_size_blocks_);
  WMESH_POINTER_CHECK(col_size_blocks_);

  WMESH_POINTER_CHECK(csr_m__);
  WMESH_POINTER_CHECK(csr_n__);
  WMESH_POINTER_CHECK(csr_nnz__);
  WMESH_POINTER_CHECK(csr_row_ptr__);
  WMESH_POINTER_CHECK(csr_row_ind__);


  //
  // Compute the number of rows and columns.
  //
  wmesh_int_t num_rows 		= 0;
  wmesh_int_t num_columns 	= 0;
  {
    for (wmesh_int_t cell_type = 0;cell_type < c2n_size_;++cell_type)
      {
	num_rows += c2n_n_[cell_type] * row_size_blocks_[cell_type];
      }
    for (wmesh_int_t cell_type = 0;cell_type < c2n_size_;++cell_type)
      {
	num_columns += c2n_n_[cell_type] * col_size_blocks_[cell_type];
      }
  }

  //
  // Allocate array csr_row_ptr.
  //
  csr_row_ptr__[0] = (wmesh_int_p)malloc(num_rows + 1,sizeof(wmesh_int_t) );
  if (!csr_row_ptr__[0])
    {
      WMESH_CHECK_STATUS(WMESH_CHECK_ERROR_MEMORY);
    }
  wmesh_int_p csr_row_ptr_ = csr_row_ptr__[0];

  //
  // Allocate array csr_col_ind.
  //
  wmesh_int_t global_shifts[5];
  global_shifts[0] = 0;
  for (wmesh_int_t cell_type = 0;cell_type < c2n_size_;++cell_type)
    {
      global_shifts[cell_type+1] = global_shifts[cell_type] + c2n_n_[cell_type];
    }

  wmesh_int_t global_dof_shifts[5];
  global_dof_shifts[0] = 0;
  for (wmesh_int_t cell_type = 0;cell_type < c2n_size_;++cell_type)
    {
      global_dof_shifts[cell_type+1] = global_dof_shifts[cell_type] + c2n_n_[cell_type] * row_size_blocks_[cell_type];
    }

  wmesh_int_t nnz = 0;  
  //
  // Get the number of interior faces.
  //
  csr_col_ind__[0] = (wmesh_int_p)malloc(nnz * sizeof(wmesh_int_t) );
  if (!csr_col_ind__[0])
    {
      WMESH_CHECK_STATUS(WMESH_CHECK_ERROR_MEMORY);
    }
  wmesh_int_p csr_col_ind_ = csr_col_ind__[0];



  
  
  
  wmesh_int_t at = 0;
  wmesh_int_t nei_ids[8];
  wmesh_int_t block_row_size;
  wmesh_int_t block_col_ind[8];

  wmesh_int_t row_idx = 0;
  
  //
  // For each cell type.
  //
  for (wmesh_int_t cell_type = 0;cell_type < c2n_size_;++cell_type)
    {
      wmesh_int_t row_idx_start = global_dof_shifts[cell_type];
      wmesh_int_t ib_size 	= row_size_blocks_[cell_type];      
      const wmesh_int_t c2n_n 	= c2n_n_[cell_type];
      const wmesh_int_t c2c_m 	= c2c_m_[cell_type];
      
      //
      // Traverse the adjacency graph.
      //
      for (wmesh_int_t cell_idx=0;cell_idx < c2n_n;++cell_idx)
	{
	  wmesh_int_t row_idx = row_idx_start + cell_idx * ib_size;
	  //
	  // Get the global cell index.
	  //
	  wmesh_int_t global_cell_idx = global_shifts[cell_type] + cell_idx;

	  //
	  // Load nei_ids
	  //
	  for (wmesh_int_t lidx=0;lidx<c2c_m;++lidx)
	    { 
	      nei_ids[lidx] = c2c_v_[c2c_ptr_[cell_type] + c2c_ld_[cell_type] * cell_idx + lidx];
	    }

	  //
	  // Set the sizes.
	  //
	  block_row_size = 0;
	  for (wmesh_int_t lidx=0;lidx<c2c_m;++lidx)
	    { 
	      if (nei_ids[lidx]>0)
		{
		  block_colidx1[block_row_size++] = nei_ids[lidx] - 1;
		}
	    }
	  block_colidx1[block_row_size++] = global_cell_idx;
	  
	  //
	  // Sort.
	  //
	  for (wmesh_int_t li=0;li<block_row_size;++li)
	    {
	      wmesh_int_t imin = li;
	      wmesh_int_t xmin = block_colidx1[li];
	      for (wmesh_int_t lj=li+1;lj<c2c_m;++lj)
		{
		  if (block_colidx1[lj] < xmin)
		    {
		      xmin =  block_colidx1[lj];
		      imin = lj;
		    }
		}
	      if (imin != li)
		{
		  wmesh_int_t tmp = block_colidx1[imin];
		  block_colidx1[imin] = block_colidx1[li];
		  block_colidx1[li] = tmp;
		}
	    }

	  wmesh_int_t jb_size 	= col_size_blocks_[nei_cell_type];	  
	  wmesh_int_t ib 	= global_cell_idx;
	  for (wmesh_int_t i=0;i<ib_size;++i)
	    {
	      for (wmesh_int_t k=0;k<block_row_size;++k)
		{
		  wmesh_int_t jb = block_colidx1[k];
		  for (wmesh_int_t j=0;j<jb_size;++j)
		    {
		      col_ind[row_ptr[at]++] = dof_col_shift_ + ;
		    }
		}
	    }
	}      
    }

  bsr_m__[0] = ncells_;
  bsr_n__[0] = ncells_;

}



wmesh_status_t wmesh_bsrjacobian(wmesh_int_t 	ncells_,
				     
				     wmesh_int_t 	c2n_size_,
				     const_wmesh_int_p 	c2n_ptr_,
				     const_wmesh_int_p 	c2n_m_,
				     const_wmesh_int_p 	c2n_n_,
				     const_wmesh_int_p 	c2n_v_,
				     const_wmesh_int_p 	c2n_ld_,
				     
				     wmesh_int_t 	c2c_size_,
				     const_wmesh_int_p 	c2c_ptr_,
				     const_wmesh_int_p 	c2c_m_,
				     const_wmesh_int_p 	c2c_n_,
				     const_wmesh_int_p 	c2c_v_,
				     const_wmesh_int_p 	c2c_ld_,

				     const_wmesh_int_p 	size_blocks_,
				     
				     wmesh_int_p 	bsr_mb__,
				     wmesh_int_p 	bsr_nb__,
				     wmesh_int_p 	bsr_nnzb__,
				     double ** 		bsr_val__,
				     wmesh_int_p*	bsr_row_ptr__,
				     wmesh_int_p*	bsr_col_ind__,
				     wmesh_int_p 	bsr_dim__)
{
  
  //  DG_JACOBIAN(I 	nelm_,
  //	      I 	nfaceinelm_,
  //	      I 	size_block_,
  //	      cst_pI 	adj_,
  //	      cst_pI 	adjoff_)
  
  bsr_mb__[0] = ncells_;
  bsr_nb__[0] = ncells_;

  //
  //
  //
  
  {
    this->m_size_block 			= size_block_;
    this->m_size_blockXsize_block 	= size_block_*size_block_;
    this->m_nelm			= nelm_;
    this->m_nfaceinelm			= nfaceinelm_;
    
    this->m_n 				= this->m_size_block * this->m_nelm;
    this->m_nc 				= this->m_size_block * this->m_size_block * (this->m_nfaceinelm+1) * this->m_nelm;
    this->m_begin 			= (pI)calloc(this->m_n+1,sizeof(I));
    this->m_index 			= (pI)malloc(this->m_nc*sizeof(I));
    this->m_values 			= (pR)malloc(this->m_nc*sizeof(R));
    
    for (I ielm=0;ielm<nelm_;++ielm)
      {
	//
	// Diagonal.
	//
	for (I k = 0;k <  this->m_size_block;++k)
	  {
	    m_begin[ielm * m_size_block + k + 1] += m_size_block;
	  }
#if 1
	//
	// Extra diagonal.
	//
	for (I localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    I nei = adj_[adjoff_[0]*ielm+localFaceIndex];
	    if (nei)
	      {
		for (I k = 0;k <  this->m_size_block;++k)
		  {
		    m_begin[ielm * m_size_block + k + 1] += m_size_block;
		  }
	      }
	  }
#endif
      }
    
    for (I i=2;i<=this->m_n;++i)
      {
	m_begin[i]+=m_begin[i-1];
      }
#if 0
    fprintf(stdout," begin[" ifmt "] = " ifmt "\n",i,m_begin[i]);
    fprintf(stdout," " ifmt " \n",m_begin[this->m_n]);
    exit(1);
#endif
    for (I ielm=0;ielm<nelm_;++ielm)
      {
#if 1
	//
	// Before diagonal.
	//
	for (I localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    I nei = adj_[adjoff_[0]*ielm+localFaceIndex];
	    if (nei)
	      {
		if (nei-1 < ielm)
		  {
		    for (I k = 0;k <  this->m_size_block;++k)
		      {
			for (I j = 0;j <  this->m_size_block;++j)
			  {					    
			    m_index[m_begin[ielm * m_size_block + k]] = (nei-1) * m_size_block + j;
			    m_begin[ielm * m_size_block + k]+=1;
			  }
		      }
		  }
	      }
	  }	
#endif	
	//
	// Diagonal.
	//
	for (I k = 0;k <  this->m_size_block;++k)
	  {
	    for (I j = 0;j <  this->m_size_block;++j)
	      {		
		m_index[m_begin[ielm * m_size_block + k]] = ielm * m_size_block + j;
		m_begin[ielm * m_size_block + k]+=1;
	      }
	  }
#if 1
	//
	// After diagonal.
	//
	for (I localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    I nei = adj_[adjoff_[0]*ielm+localFaceIndex];
	    if (nei)
	      {
		if (nei - 1 > ielm)
		  {
		    for (I k = 0;k <  this->m_size_block;++k)
		      {
			for (I j = 0;j <  this->m_size_block;++j)
			  {					    
			    m_index[m_begin[ielm * m_size_block + k]] = (nei-1) * m_size_block + j;
			    m_begin[ielm * m_size_block + k]+=1;
			  }
		      }
		  }
	      }
	  }
#endif	
      }

    for (I i=this->m_n;i>0;--i)
      {
	m_begin[i] = m_begin[i-1];
      }
    m_begin[0] = 0;

    
    for (I i=0;i<this->m_n;++i)
      {
	std::sort(&m_index[m_begin[i]],
		  &m_index[m_begin[i+1]],
		  [](const wmesh_int_t& a,const wmesh_int_t& b){return a < b;});
	
	//	qsort(&m_index[m_begin[i]],m_begin[i+1]-m_begin[i],sizeof(I),comp);
      }
    //    this->spy("dg.txt");

}
#endif

#if 0
#include "DG_VAR.hpp"
#include "DG_JACOBIAN.hpp"
  // #include "DG_BOUNDARY_CONDITION.hpp"
#endif  

  int main(int argc, char ** argv)
  {
#if 0
#if 0
#endif

    std::cout << "####### " << "####" << std::endl;    
    for (wmesh_int_t element = WMESH_ELEMENT_EDGE;element < WMESH_ELEMENT_ALL;++element)
      {
	wmesh_cubature_t<double> cubature;
	wmesh_status_t status =  wmesh_cubature_def(&cubature,
						    element,
						    WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE,
						    9);

	double fd = 0;
	for (wmesh_int_t i=0;i<cubature.m_w.n;++i)
	  {
	    fd += cubature.m_w.v[i];
	  }
	std::cout << "sum weights " << fd << std::endl;
      }
    
    
    { 
      wmesh_cubature_boundary_t<double> cubature_boundary;
      const wmesh_int_t num_facets = 3,elm = WMESH_ELEMENT_TRIANGLE;
      status =  wmesh_cubature_boundary_def(&cubature_boundary,
					     elm,
					    WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE,
					    9);
      WMESH_STATUS_CHECK(status);
      for (wmesh_int_t i=0;i<num_facets;++i)
	{
	  
	  std::cout << "facet " << i << std::endl;

#if 1
	  for (wmesh_int_t j=0;j<2*2;++j)
	    {
	      std::cout << "rot " << j << std::endl;
	      double fd = 0;
	      for (wmesh_int_t k=0;k<cubature_boundary.m_facets_cubature[i][j].m_w.n;++k)
		{
		  double det = (i==1) ? 1.0/sqrt(2.0) : 1.0;
		  fd += cubature_boundary.m_facets_cubature[i][j].m_w.v[k] * 1.0 * det;
		}
	      
	      std::cout << "sum weights over facet " << fd << std::endl;
	      //	      std::cout << cubature_boundary.m_facets_cubature[i][j].m_c << std::endl;
	    }
#endif
	}
    }
  
#if 0
  {
    wmesh_cubature_boundary_t<double> cubature_boundary;
    status =  wmesh_cubature_boundary_def(&cubature_boundary,
					  WMESH_ELEMENT_TRIANGLE,
					  WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE,
					  9);
    WMESH_STATUS_CHECK(status);
    for (wmesh_int_t i=0;i<3;++i)
      {
	std::cout << "edge " << i << std::endl;
	for (wmesh_int_t j=0;j<2*2;++j)
	  {
	    std::cout << "rot " << j << std::endl;
	    std::cout << cubature_boundary.m_facets_cubature[i][j].m_c << std::endl;
	  }
      }
  }

#endif
  
#if 0
  { 
    wmesh_cubature_boundary_t<double> cubature_boundary;
    status =  wmesh_cubature_boundary_def(&cubature_boundary,
					  WMESH_ELEMENT_TETRAHEDRON,
					  WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE,
					  3);
    WMESH_STATUS_CHECK(status);
    for (wmesh_int_t i=0;i<4;++i)
      {
	std::cout << "face " << i << std::endl;
	for (wmesh_int_t j=0;j<3*2;++j)
	  {
	    std::cout << "rot " << j << std::endl;
	    std::cout << cubature_boundary.m_facets_cubature[i][j].m_c << std::endl;
	  }
      }
  }
  { 
    wmesh_cubature_boundary_t<double> cubature_boundary;
    status =  wmesh_cubature_boundary_def(&cubature_boundary,
					  WMESH_ELEMENT_HEXAHEDRON,
					  WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE,
					  3);
    WMESH_STATUS_CHECK(status);
    for (wmesh_int_t i=0;i<6;++i)
      {
	std::cout << "face " << i << std::endl;
	for (wmesh_int_t j=0;j<4*2;++j)
	  {
	    std::cout << "rot " << j << std::endl;
	    std::cout << cubature_boundary.m_facets_cubature[i][j].m_c << std::endl;
	  }
      }
  }
#endif
  

  exit(1);
#endif
    wmesh_t* 		mesh 		= nullptr;
    wmesh_status_t 	status;
  
  //
  // Parameters.
  //
  WCOMMON::cmdline::str_t 	ofilename;
  const char * 			ifilename 	= nullptr;
  bool 				verbose 	= false;
  wmesh_int_t 			degree		= 0;

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
   
    //
    // Get the degree of the dg space.
    //
    if (false == cmd.option("-d", &degree))
      {
	fprintf(stderr,"missing output file, '-d' option.\n");
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
    
    if (cmd.get_nargs() == 1)
      {
	fprintf(stderr,"no file found.\n");
	return WMESH_STATUS_INVALID_ARGUMENT;
      }
    
    ifilename = cmd.get_arg(1);
  }
  
  //
  // Read the mesh.
  //
  status = wmesh_read(&mesh,
		      ifilename);
  WMESH_STATUS_CHECK(status);


  //
  // Analyze the mesh.
  //
  status = wmesh_analysis(mesh);
  
  WMESH_STATUS_CHECK(status);

  wmesh_shape_info_t shape_info_element;
  wmesh_shape_info_t shape_info_f;
  wmesh_shape_info_t shape_info_u;
  wmesh_shape_info_t shape_info_test;
  
  const wmesh_int_t u_degree 		= 2;
  const wmesh_int_t u_family 		= WMESH_SHAPE_FAMILY_LAGRANGE;
  const wmesh_int_t u_nodes_family 	= WMESH_NODES_FAMILY_LAGRANGE;

  const wmesh_int_t f_degree 		= degree;
  const wmesh_int_t f_family 		= WMESH_SHAPE_FAMILY_LAGRANGE;
  const wmesh_int_t f_nodes_family 	= WMESH_NODES_FAMILY_LAGRANGE;
  
  status = wmesh_shape_info_def(&shape_info_element,WMESH_SHAPE_FAMILY_LAGRANGE, 1);
  WMESH_STATUS_CHECK(status);
  
  status = wmesh_shape_info_def(&shape_info_u,u_family, u_degree);
  WMESH_STATUS_CHECK(status);
  
  status = wmesh_shape_info_def(&shape_info_f,f_family, f_degree);
  WMESH_STATUS_CHECK(status);
  
  status = wmesh_shape_info_def(&shape_info_test,f_family, f_degree);
  WMESH_STATUS_CHECK(status);
  
  //
  // Generate the finite element space for the velocity.
  // 
  wmeshspace_t * velocity_space;
  status = wmeshspace_def(&velocity_space,
			  u_nodes_family,
			  u_degree,
			  mesh);
  WMESH_STATUS_CHECK(status);
  
  wmesh_int_t csr_size 	= 0;
  wmesh_int_p csr_ptr 	= nullptr;
  wmesh_int_p csr_ind 	= nullptr;
  wmeshspacedg_sparse(mesh,
		      degree,
		      &csr_size,
		      &csr_ptr,
		      &csr_ind);

  double *		rhs 	= (double*) malloc(sizeof(double) * csr_size);
  double *		csr_val = (double*) malloc(sizeof(double) * csr_ptr[csr_size]);
  
#if 0
  
  //
  // Spy the symbolic matrix.
  //
  std::cout << "size " << csr_size << std::endl;
  std::cout << "nnz  " << csr_ptr[csr_size] << std::endl;

  FILE * f = fopen("out.txt","w");
  for (wmesh_int_t i=0;i<csr_size;++i)
    {
      for (wmesh_int_t s=csr_ptr[i];s<csr_ptr[i+1];++s)
	{
	  wmesh_int_t j = csr_ind[s];
	  fprintf(f," " WMESH_INT_FORMAT " " WMESH_INT_FORMAT "\n",csr_size - i, j);
	  //	  std::cout << csr_size - i << " " << j << std::endl;
	}
    }
  fclose(f);
#endif

  const wmesh_int_t topodim = mesh->m_topology_dimension;
  const wmesh_int_t 	velocity_storage 	= WMESH_STORAGE_INTERLEAVE;
  const wmesh_int_t	velocity_m 		= topodim;
  const wmesh_int_t	velocity_n 		= velocity_space->m_ndofs;
  const wmesh_int_t	velocity_ld 		= velocity_m;  
  double *		velocity 		= (double*)malloc(sizeof(double) * velocity_n * velocity_m);

  
  //
  // Interpolate the velocity.
  //
  {
    wmesh_int_t coo_dofs_m  	= topodim;
    wmesh_int_t coo_dofs_n  	= velocity_space->m_ndofs;
    wmesh_int_t coo_dofs_ld 	= coo_dofs_m;
    double * 	coo_dofs 	= (double*)malloc(sizeof(double) * coo_dofs_n * coo_dofs_m);    
    status 			=  wmeshspace_generate_coodofs(velocity_space,
							       WMESH_STORAGE_INTERLEAVE,
							       coo_dofs_m,
							       coo_dofs_n,
							       coo_dofs,
							       coo_dofs_ld);
    WMESH_STATUS_CHECK(status);
    double x[3];
    double u[3];
    for (wmesh_int_t j=0;j<coo_dofs_n;++j)
      {
	
	for (wmesh_int_t i=0;i<coo_dofs_m;++i)
	  {
	    x[i] = coo_dofs[coo_dofs_ld*j+i];
	  }

	for (wmesh_int_t i=0;i<coo_dofs_m;++i)
	  {	    
	    u[i] = 0.0 * x[i];
	  }	
	u[0] = 1.0;

	
	if (velocity_storage == WMESH_STORAGE_INTERLEAVE)
	  {
	    for (wmesh_int_t i=0;i<coo_dofs_m;++i)
	      {
		velocity[velocity_ld * j + i] = u[i];
	      }
	  }
	else 
	  {
	    for (wmesh_int_t i=0;i<coo_dofs_m;++i)
	      {
		velocity[velocity_ld * i + j] = u[i];
	      }
	  }	
      }
    
    free(coo_dofs);    
  }

  
  wmeshspacedg_t * spacedg;            
  status = wmeshspacedg_def(&spacedg,
			    f_nodes_family,
			    f_degree,
			    mesh);  
  WMESH_STATUS_CHECK(status);

  
  status = wmeshspacedg_advection(spacedg,
				  &shape_info_element,
				  &shape_info_u,
				  &shape_info_f,
				  &shape_info_test,
				  
				  velocity_space,
				  
				  velocity_storage,
				  velocity_m,
				  velocity_n,
				  velocity,
				  velocity_ld,
				  
				  csr_size,
				  csr_ptr,
				  csr_ind,
				  csr_val,				  
				  rhs);
  WMESH_STATUS_CHECK(status);

  //
  // Output.
  //  
  status =  bms_matrix_market_dense_dwrite(csr_size,
					   1,
					   rhs,
					   csr_size,
					   "%s_rhs.mtx",
					   ofilename);
  WMESH_STATUS_CHECK(status);
  
  status =  bms_matrix_market_csr_dwrite(csr_size,
					 csr_size,
					 csr_ptr[csr_size],
					 csr_ptr,
					 csr_ind,
					 csr_val,
					 "%s.mtx",
					 ofilename);
  WMESH_STATUS_CHECK(status);
  
  return WMESH_STATUS_SUCCESS;
}
