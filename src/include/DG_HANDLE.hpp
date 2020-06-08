#pragma once
#include "wmesh-types.hpp"
struct DG_HANDLE
{

  wmesh_mat_t<double> m_EVALU;
  wmesh_mat_t<double> m_UVWDOFS;
  wmesh_mat_t<double> m_BMAT;
  wmesh_mat_t<double> m_BRHS;

  wmesh_mat_t<double> m_mat_tmpbrhs_uelm;
  wmesh_mat_t<double> m_brhs_uvw;
  wmesh_mat_t<double> m_mat_tmpbrhs_uface[nfaceinelm];
  
};
