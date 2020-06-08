#pragma once
struct DG_DATA
{
  
  //
  //
  //
  wmesh_mat_t<double> m_flux{};

  //
  // submatrix.
  //
  wmesh_mat_t<double> m_flux_part_corr{};
  
  //
  // submatrix.
  //
  wmesh_mat_t<double> m_flux_part_sol{};
  
  wmesh_mat_t<double> uvw_ldofs;

  wmesh_mat_t<double> m_local_matrices;
  wmesh_mat_t<double> m_local_matrix_part_corr;  
  wmesh_mat_t<double> m_local_matrix_part_sol;

  wmesh_mat_t<double> vec_neicorr;
  wmesh_mat_t<double> vec_neisol;

  wmesh_mat_t<double> mat_belm;
  wmesh_mat_t<double> vec_nrmelm[nfaceinelm];
  wmesh_mat_t<double> vec_uface[nfaceinelm];

  wmesh_mat_t<double> m_local_rhs;
  wmesh_mat_t<double> hsol;

  double * m_udofs;
  double *  m_local_matrices_memory;
  double *  m_flux_memory{};
  double belm[dim*dim];
  double nrmelm[dim * nfaceinelm];
  double lcrhs[128];

  virtual ~DG_DATA()
  {
    if (m_flux_memory)
      {
	free(m_flux_memory);
	m_flux_memory = nullptr;
      }    
    if (m_udofs)
      {
	free(m_udofs);
	m_udofs = nullptr;
      }    
  };

  DG_DATA(wmesh_int_t teta_n_,
	  wmesh_int_t trial_n_,
	  wmesh_int_t test_n_,
	  wmesh_int_t teta_u_n_)
  {
    this->m_flux_memory = (double * )malloc(sizeof(double)*(trial_n_*test_n_+teta_n_*test_n_));
    wmesh_mat_t<double>::define(&this->m_flux_part_corr,test_n_,trial_n_,this->m_flux_memory,test_n_);
    wmesh_mat_t<double>::define(&this->m_flux_part_sol,test_n_,teta_n_,&this->m_flux_memory[trial_n_*test_n_],test_n_);    

    this->m_udofs = (double * )malloc(sizeof(double)*(teta_u_n_*dim));
    wmesh_mat_t<double>::define(&this->uvw_ldofs,teta_u_n_,dim,this->m_udofs,teta_u_n_);
    
    wmesh_mat_t<double>::define(&this->mat_belm,
		      2,2,belm,2);

    for (wmesh_int_t localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {
	wmesh_mat_t<double>::define(&this->vec_nrmelm[localFaceIndex],1,dim,&nrmelm[dim * localFaceIndex],1);
      }

    wmesh_int_t len  = trial_n_*test_n_+teta_n_*test_n_;
    m_local_matrices_memory = (double * )malloc(sizeof(double)*len);
    wmesh_mat_t<double>::define(&this->m_local_matrices,1,len,m_local_matrices_memory,1);
    wmesh_mat_t<double>::define(&this->m_local_matrix_part_corr,test_n_,trial_n_,  m_local_matrices_memory, test_n_);
    wmesh_mat_t<double>::define(&this->m_local_matrix_part_sol, test_n_,teta_n_,   m_local_matrices_memory + trial_n_ * test_n_, test_n_);
    wmesh_mat_t<double>::define(&this->m_local_rhs,1,trial_n_,lcrhs,1);
  };

};
