
#include "wmesh_nodes_family.h"
#include "bms.h"
#include <iostream>
#include <math.h>
#include "wmesh-status.h"
// #include "wfe_status.hpp"

#define MKL_ILP64 1
#include "mkl.h"
#include "wmesh.h"
#include "bms.hpp"
#define solve_params(_type) wmesh_int_t*n,_type*jacM,_type*wr,_type*wi,_type*vl,_type*vr,_type*work,wmesh_int_t*work_n_
template<typename T>
void solve(solve_params(T));

template<>
void solve<double>(solve_params(double))
{
  static constexpr const  char trN[1] = {'N'};
  static constexpr const char trV[1] = {'V'};
  wmesh_int_t ldvl=1,info;
  dgeev(trN,
	trV,
	n,
	jacM,
	n, 
	wr, 
	wi, 
	vl, 
	&ldvl, 
	vr, 
	n, 
	work, 
	work_n_,
	&info);
}

template<>
void solve<float>(solve_params(float))
{
  static constexpr const  char trN[1] = {'N'};
  static constexpr const char trV[1] = {'V'};
  wmesh_int_t ldvl=1,info;
  sgeev(trN,
	trV,
	n,
	jacM,
	n, 
	wr, 
	wi, 
	vl, 
	&ldvl, 
	vr, 
	n, 
	work, 
	work_n_,
	&info);
}

static inline wmesh_int_t flat3d(wmesh_int_t n_,wmesh_int_t l[])
{
  return l[0] + l[1]*n_ + l[2]*n_*n_;
}


static inline wmesh_int_t flat2d(wmesh_int_t n_,wmesh_int_t l[])
{
  return l[0] + l[1]*n_;
}



template<typename T>
static wmesh_status_t bms_simplex_1dto3d(wmesh_int_t 		d_,
					 wmesh_int_t 		n_,
					 const T*__restrict__ 	w_,
					 T *			rst_,
					 wmesh_int_t 		rst_ld_,
					 wmesh_int_p size_)
{
  static constexpr const T s_T0(0.0);
  static constexpr const T s_T1(1.0);
  static constexpr const T s_T2(2.0);
  static constexpr const T s_T3(3.0);
  static constexpr const T s_T4(4.0);  
  wmesh_int_t m 	= n_-1;
  wmesh_int_t next 	= 0;
  for (wmesh_int_t i=0;i<m+1;++i)
    {
      for (wmesh_int_t j=0;j<m+1-i;++j)
	{
	  wmesh_int_t l = m-i-j;	  
	  rst_[rst_ld_*next+0] = (d_*s_T1+s_T2*w_[i]-w_[j]-w_[l]) / s_T3;
	  rst_[rst_ld_*next+1] = (d_*s_T1+s_T2*w_[j]-w_[i]-w_[l]) / s_T3;
	  rst_[rst_ld_*next+2] = s_T0;
	  ++next;
	}
    }
  
  for (wmesh_int_t j=1;j<=m;++j)
    {
      for (wmesh_int_t k=2;k<=m+2-j;++k)
	{
	  wmesh_int_t l = m+3-j-k;
	  rst_[rst_ld_*next+0] = s_T0;
	  rst_[rst_ld_*next+1] = (d_*s_T1+s_T2*w_[j-1]-w_[k-1]-w_[l-1]) / s_T3;
	  rst_[rst_ld_*next+2] = (d_*s_T1+s_T2*w_[k-1]-w_[j-1]-w_[l-1]) / s_T3;
	  ++next;
	  
	}
    }

  for (wmesh_int_t i=2;i<=m;++i)
    {
      for (wmesh_int_t k=2;k<=m+2-i;++k)
	{
	  wmesh_int_t l = m+3-i-k;
	  rst_[rst_ld_*next+0] = (d_*s_T1+s_T2*w_[i-1]-w_[k-1]-w_[l-1]) / s_T3;
	  rst_[rst_ld_*next+1] = s_T0;
	  rst_[rst_ld_*next+2] = (d_*s_T1+s_T2*w_[k-1]-w_[i-1]-w_[l-1]) / s_T3;
	  ++next;
	}
    }
  
  for (wmesh_int_t i=2;i<=m;++i)
    {
      for (wmesh_int_t j=2;j<=m+1-i;++j)
	{
	  wmesh_int_t l = m+3-i-j;
	  T xi = (s_T1*d_+s_T2*w_[i-1]-w_[j-1]-w_[l-1]) / s_T3;
	  T eta = (s_T1*d_ +s_T2*w_[j-1]-w_[i-1]-w_[l-1]) / s_T3;
	  rst_[rst_ld_*next+0] = xi;
	  rst_[rst_ld_*next+1] = eta;
	  rst_[rst_ld_*next+2] = s_T1*d_-xi-eta;
	  ++next;
	}
    }

  for (wmesh_int_t i=2;i<=m;++i)
    {
      for (wmesh_int_t j=2;j<=m+1-i;++j)
	{
	  for (wmesh_int_t k=2;k<=m+2-i-j;++k)
	    {
	      wmesh_int_t l = m+4-i-j-k;
	      rst_[rst_ld_*next+0] = (d_*s_T1+s_T3*w_[i-1]-w_[j-1]-w_[k-1]-w_[l-1]) / s_T4;
	      rst_[rst_ld_*next+1] = (d_*s_T1+s_T3*w_[j-1]-w_[i-1]-w_[k-1]-w_[l-1]) / s_T4;
	      rst_[rst_ld_*next+2] = (d_*s_T1+s_T3*w_[k-1]-w_[i-1]-w_[j-1]-w_[l-1]) / s_T4;
	      ++next;
	    }
	}
    }
  size_[0] = next;
  return WMESH_STATUS_SUCCESS;  
}


template<typename T>
static wmesh_status_t bms_simplex_1dto2d(const wmesh_int_t 	n_,
					 const T*		w_,
					 wmesh_int_t		incw_,
					 wmesh_int_t		rs_storage_,
					 T*			rs_,
					 wmesh_int_t		rs_ld_) noexcept
{
  static constexpr const T zero(0.0);
  static constexpr const T one(1.0);
  static constexpr const T two(2.0);
  static constexpr const T three(3.0);
  static constexpr const T six(6.0);
	  
  const unsigned int numPointsOnEdge 	= (n_ > 2) ? n_ - 2  : 0;
  auto _degree = n_-1;
  unsigned int startedge0 = 3;
  unsigned int startedge1 = 3 + numPointsOnEdge + numPointsOnEdge-1;  
  unsigned int startedge2 = 3 + numPointsOnEdge*2 + numPointsOnEdge-1;
  unsigned int startInterior = 3*numPointsOnEdge+3;
  rs_[0*rs_ld_+0] 	= zero;
  rs_[0*rs_ld_+1] 	= zero;
  rs_[1*rs_ld_+0] 	= one;
  rs_[1*rs_ld_+1] 	= zero;
  rs_[2*rs_ld_+0] 	= zero;
  rs_[2*rs_ld_+1] 	= one;

  //
  // Third edge
  //
  for (unsigned int j=1;j<_degree;++j)
    {
      const auto wj 	= w_[j*incw_];
      const auto wk 	= w_[(_degree-j)*incw_];
      rs_[startedge2*rs_ld_+0] 	= zero;
      rs_[startedge2*rs_ld_+1] 	= (three + two * wj - wk) / six;
      --startedge2;
    }
	  
  //
  // First edge
  //
  for (unsigned int i=1;i<_degree;++i)
    {
      const auto wi 	= w_[i*incw_];
      const auto wk 	= w_[(_degree - i)*incw_];
      rs_[startedge0*rs_ld_+0] 	= (three + two * wi - wk) / six;
      rs_[startedge0*rs_ld_+1] 	= zero;
      ++startedge0;
    }
	  
  startedge0=3;
  // 
  // Second edge
  //
  for (unsigned int i=1;i<_degree;++i)
    {
      rs_[startedge1*rs_ld_+0] 	= rs_[startedge0++*rs_ld_+0]; 
      rs_[startedge1*rs_ld_+1] 	= rs_[++startedge2*rs_ld_+1];
      --startedge1;
    }
      
  //
  // Interior
  //
  for (unsigned int i=1;i<_degree;++i)
    {
      const auto wi = w_[i];
      for (unsigned int j=1;j<_degree-i;++j)
	{
	  const auto wj 		= w_[j*incw_];
	  const auto wk 		= w_[(_degree-i-j)*incw_];
	  rs_[startInterior*rs_ld_+0] 	= ( two* (one + wi) - (wj + wk) ) / six;
	  rs_[startInterior*rs_ld_+1] 	= ( two* (one + wj) - (wi + wk) ) / six;
	  startInterior++;
	}
    }
  return WMESH_STATUS_SUCCESS;
};






template<typename T>
wmesh_status_t wmesh_nodes_legendre(wmesh_int_t 		nspl_,
				    wmesh_int_t 		storagep_,
				    T*__restrict__ 		p_,
				    wmesh_int_t 		ldp_, 
				    T*__restrict__ 		w_,
				    wmesh_int_t 		incw_, 
				    wmesh_int_t 		work_n_,
				    T*__restrict__		work_)
{
  static constexpr T rhalf = static_cast<T>(0.5);
  static constexpr T r0 = static_cast<T>(0.0);
  static constexpr T r1 = static_cast<T>(1.0);
  static constexpr T r2 = static_cast<T>(2.0);
  static constexpr T r4 = static_cast<T>(4.0);
  
  //  if (WFE::storage_t::is_invalid(storagep_)) return WMESH_STATUS_INVALID_ENUM;
  
  if (nspl_ < 0)    	return WMESH_STATUS_INVALID_SIZE;
  if (!p_)      	return WMESH_STATUS_INVALID_POINTER;
  if (ldp_  < 1) 	return WMESH_STATUS_INVALID_SIZE;
  if (!w_)      	return WMESH_STATUS_INVALID_POINTER;
  if (incw_ < 1) 	return WMESH_STATUS_INVALID_SIZE;
  
  if (work_n_ < (nspl_+2)*(6 + nspl_) + 4*(nspl_+1) ) return WMESH_STATUS_INVALID_SIZE;

  bool p_is_block = (storagep_ == WMESH_STORAGE_BLOCK);
  if (p_is_block && ldp_ < nspl_) return WMESH_STATUS_INVALID_SIZE;  

  wmesh_int_t n    	= nspl_+1;
  wmesh_int_t nn 	= nspl_+2;
  T 
    * __restrict__ vl = NULL,
    * __restrict__ b    = &work_[nn],
    * __restrict__ jacM = &work_[2*nn],
    * __restrict__ wr   = &work_[2*nn+nn*nn],
    * __restrict__ wi   = &work_[2*nn+nn*nn+nn],
    * __restrict__ vr   = &work_[2*nn+nn*nn+2*nn],
    * __restrict__ work  = &work_[2*nn+2*nn*nn+2*nn];   
  for (wmesh_int_t k=0;k<n;++k)
    {
      b[k]=(T)0.0;
    }
  
  b[0] = r2;
  b[1] = r0;
  for (wmesh_int_t k=2;k<n;++k)
    {
      b[k] = r1 / (r4 - r1/( ((T)(k-1))  * ((T)(k-1)) ) );
    }
  
  for (wmesh_int_t i=0;i<n;++i)
    {
      for (wmesh_int_t j=0;j<n;++j)
	{
	  jacM[i+j*n]=r0;
	}
    }
  
  for (wmesh_int_t i=0;i<n;++i)
    {
      jacM[i+i*n]=(T)0.0;
    }
  for (wmesh_int_t i=0;i<n-1;++i)
    {
      jacM[i+i*n+1]=sqrt(b[1+i]);
    }
  for (wmesh_int_t i=1;i<n;++i)
    {
      jacM[i+i*n-1]=sqrt(b[i]);
    }

  //
  // Eigen values / vectors decompositi
  //
  solve(&n,jacM,wr,wi,vl,vr,work,&work_n_);
 
#if 0
  {
    static constexpr const  char trN[1] = {'N'};
    static constexpr const char trV[1] = {'V'};
    wmesh_int_t ldvl=1,info;
    dgeev(trN,
	  trV,
	  &n,
	  jacM,
	  &n, 
	  wr, 
	  wi, 
	  vl, 
	  &ldvl, 
	  vr, 
	  &n, 
	  work, 
	  &work_n_,
	  &info);
  }
#endif  
  if (p_is_block)
    {
      for (wmesh_int_t i=0;i<n-1;++i) 
	{
	  p_[i] = wr[i];
	}
    }
  else
    {
      for (wmesh_int_t i=0;i<n-1;++i) 
	{
	  p_[i * ldp_] = wr[i];
	}	  
    }
  
  for (wmesh_int_t i=0;i<n-1;++i) 
    {
      w_[i * incw_] = vr[1+i*n]*vr[1+i*n] * r2;
    }

  for (int i = 0; i < (n-1)-1; i++)
    {
      int min_idx = i;
      for (int j = i+1; j < (n-1); j++)
	{
	  if (p_[j] < p_[min_idx])
	    {
	      min_idx = j;
	    }
	}
      
      { auto tmp = p_[min_idx]; p_[min_idx] = p_[i]; p_[i] = tmp; }
      { auto tmp = w_[min_idx]; w_[min_idx] = w_[i]; w_[i] = tmp; }
      
    }
  
  //
  // Force symmetry.
  //
  if (nspl_ % 2 > 0)
    {
      p_[ (nspl_/2) *ldp_] = 0.0;
    }

  if (p_is_block)
    {
      for (int i=0;i<nspl_/2;++i)
	{
	  T tmp = (p_[i ] - p_[(nspl_-1-i) ])* rhalf;
	  p_[(nspl_-1-i)] = -tmp;
	  p_[i] = tmp;
	}
    }
  else
    {
      for (int i=0;i<nspl_/2;++i)
	{
	  T tmp = (p_[i *ldp_] - p_[(nspl_-1-i) *ldp_])*0.5;
	  p_[(nspl_-1-i) *ldp_] = -tmp;
	  p_[i *ldp_] = tmp;
	}
    }

  return WMESH_STATUS_SUCCESS;
}



template<typename T>
wmesh_status_t
bms_nodes_pyramid(wmesh_int_t		degree_,

		  wmesh_int_t		c_storage_,
		  wmesh_int_t		c_m_,
		  wmesh_int_t		c_n_,
		  T*__restrict__ 	c_v_,
		  wmesh_int_t		c_ld_,

		  wmesh_int_t 		level2_,
		  T 			shift_,
		  T 			scal_,
		  wmesh_int_p 		perm_,
		  wmesh_int_t		work_n_,
		  T* __restrict__ 	work_)  
{
  
  static const wmesh_int_t s_rst_m  = 3;
  static const wmesh_int_t s_rst_n  = 5;
  static const wmesh_int_t s_rst_ld = s_rst_m; 
  static const wmesh_int_t s_rst_v[] =  { 0,0,0,
					  1,0,0,
					  1,1,0,
					  0,1,0,
					  0,0,1};
  
  static const wmesh_int_t s_t2n[] = {0,1,4,
				      1,2,4, 
				      2,3,4, 
				      3,0,4};

  static const wmesh_int_t s_q2n[] = {0,3,2,1};
  
  static const wmesh_int_t s_e2n[] = {0,1,
				      1,2,
				      2,3,
				      3,0,
				      0,4,
				      1,4,
				      2,4,
				      3,4};

  static const wmesh_int_t s_t2n_diag[] = {0,2,4};
  static const wmesh_int_t s_tet2n[] = {1,2,0,4,
					3,0,2,4};

  wmesh_int_t 		r_n_ 	= (degree_+1);

  

  wmesh_int_t r_inc_ = 1;
  //  const wmesh_int_t   	n 	= (degree_>0)? degree_-1 : 0;
  wmesh_status_t status;
  T r_v_[512],w[512];
  T one = 1.0;
  T two	= 2.0;
  T three= 3.0;
  T six	= 6.0;
  T eight= 8.0;
  
  T z_v_[512];
  z_v_[0]  = -1.0;
  status =  wmesh_nodes_legendre(degree_-1,
				 WMESH_STORAGE_INTERLEAVE,
				 z_v_+1,
				 1, 
				 w,
				 1,
				 work_n_,
				 work_);
  z_v_[degree_]  = 1.0; 
  WMESH_STATUS_CHECK(status);
  
    wmesh_int_t ijk[3];
    wmesh_int_t degree = degree_;
  wmesh_int_t level = 0;
  {
    
    //      shift_ = (z_v_[level] + one) / two;
    //      scal_ = (z_v_[degree_ - level] + one) / two - shift_;
    
    std::cout << "degree " << degree << std::endl;

    {
      r_v_[0]  = -1.0;
      status =  wmesh_nodes_legendre(degree-1,
				     WMESH_STORAGE_INTERLEAVE,
				     r_v_+1,
				     1, 
				       w,
				       1,
				       work_n_,
				       work_);
	r_v_[degree_]  = 1.0; 
	WMESH_STATUS_CHECK(status);
      }
      
      //
      // nodes
      //
      for (wmesh_int_t j=0;j<s_rst_n;++j)
	{
	  
	  for (wmesh_int_t i=0;i<s_rst_m;++i)
	    {
	      ijk[i] = level + s_rst_v[s_rst_ld*j+i] * degree;
	    }
	  
	  wmesh_int_t idx = perm_[flat3d(r_n_,ijk)]-1;
	  for (wmesh_int_t i=0;i<s_rst_m;++i)
	    {
	      c_v_[c_ld_*idx+i] = shift_ + s_rst_v[s_rst_ld*j+i] * scal_;
	    }
	}
      
      
      //
      // edges
      //
      for (wmesh_int_t ie=0;ie<8;++ie)
	{
	  wmesh_int_t node_0 = s_e2n[ie*2+0];
	  wmesh_int_t node_1 = s_e2n[ie*2+1];
	  std::cout << "node " << node_0  << " " << node_1 << std::endl;
	  for (wmesh_int_t i=1;i<degree;++i)
	    {
	      T l0 	= one - (r_v_[ r_inc_* i ] + one) / two;
	      T l1 	= (r_v_[ r_inc_* i ] + one) / two;
	      
	      wmesh_int_t s0 = degree - i;
	      wmesh_int_t s1 = i;
	      
	      for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
		{
		  ijk[idim] = level + s0 * s_rst_v[s_rst_ld*node_0 + idim] + s1 * s_rst_v[s_rst_ld*node_1 + idim];
		}
	      
	      wmesh_int_t idx = perm_[flat3d(r_n_,ijk)]-1;
	      
	      for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
		{
		  c_v_[c_ld_*idx+idim] = shift_ + (l0 * s_rst_v[s_rst_ld*node_0 + idim] + l1 * s_rst_v[s_rst_ld*node_1 + idim]) * scal_;
		  //	      std::cout << "idim  " << idim << " = " << c_v_[c_ld_*idx+idim] << std::endl;
		}	  
	    }
	}
      
      for (wmesh_int_t it=0;it<4;++it)
	{
	  wmesh_int_t node_0 = s_t2n[it*3+0];
	  wmesh_int_t node_1 = s_t2n[it*3+1];
	  wmesh_int_t node_2 = s_t2n[it*3+2];
	  
	  for (wmesh_int_t i=1;i<degree_;++i)
	    {
	      const auto ri = r_v_[i];
	      for (wmesh_int_t j=1;j<degree_-i;++j)
		{
		  const auto rj	= r_v_[j*r_inc_];
		  const auto rl	= r_v_[(degree_ - i - j)*r_inc_];

		  T l1  = ( two * (one + ri) - (rj + rl) ) / six;
		  T l2  = ( two * (one + rj) - (ri + rl) ) / six;
		  T l0  = one - l1 - l2;

		  wmesh_int_t s1  = i;
		  wmesh_int_t s2  = j;
		  wmesh_int_t s0  = degree_-s1 - s2;
		  for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
		    ijk[idim] = level +
		      s0 * s_rst_v[s_rst_ld*node_0 + idim] +
		      s1 * s_rst_v[s_rst_ld*node_1 + idim] +
		      s2 * s_rst_v[s_rst_ld*node_2 + idim];
	      
		  wmesh_int_t idx = perm_[flat3d(r_n_,ijk)]-1;	  
		  for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
		    c_v_[idx * c_ld_ + idim] =
		      shift_ + 
		      ( l0 * s_rst_v[s_rst_ld*node_0 + idim] +
			l1 * s_rst_v[s_rst_ld*node_1 + idim] +
			l2 * s_rst_v[s_rst_ld*node_2 + idim] ) * scal_;
		}
	    }
	}
  

      for (wmesh_int_t iq=0;iq<1;++iq)
	{
	  wmesh_int_t node_0 = s_q2n[iq*4+0];
	  wmesh_int_t node_1 = s_q2n[iq*4+1];
	  wmesh_int_t node_2 = s_q2n[iq*4+2];
	  wmesh_int_t node_3 = s_q2n[iq*4+3];
	
	  for (wmesh_int_t i=1;i<degree_;++i)
	    {
	      const auto ri = (r_v_[i*r_inc_]+one)/two;
	      for (wmesh_int_t j=1;j<degree_;++j)
		{
		  const auto rj = (r_v_[j*r_inc_]+one) / two;
	       
		  T l0  = (one - ri)*(one - rj);
		  T l1   = ri*(one - rj);
		  T l2  = ri*rj;
		  T l3  = (one - ri)*rj;
#if 0	     
		  wmesh_int_t s0  = (degree_ - i)*(degree_ - j);
		  wmesh_int_t s1  = i*(degree_ - j);
		  wmesh_int_t s2  = i * j ;
		  wmesh_int_t s3  = (degree_ - i) * j;
#endif
		  ijk[0] = j;
		  ijk[1] = i;
		  ijk[2] = 0;
#if 0
		  for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
		    ijk[idim] = level + 
		      s0 * s_rst_v[s_rst_ld*node_0 + idim] +
		      s1 * s_rst_v[s_rst_ld*node_1 + idim] +
		      s2 * s_rst_v[s_rst_ld*node_2 + idim] + 
		      s3 * s_rst_v[s_rst_ld*node_3 + idim];
#endif	      
		  std::cout << "ddddddddddddddddd " << ijk[0] << " " << ijk[1] << " " << ijk[2] << std::endl;
		  wmesh_int_t idx = perm_[flat3d(r_n_,ijk)]-1;
	      
		  for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
		    c_v_[idx * c_ld_ + idim] =
		      shift_ + (l0 * s_rst_v[s_rst_ld*node_0 + idim] +
				l1 * s_rst_v[s_rst_ld*node_1 + idim] +
				l2 * s_rst_v[s_rst_ld*node_2 + idim] + 
				l3 * s_rst_v[s_rst_ld*node_3 + idim]) * scal_;		
		}
	    }
	}	
    }


  {      
    wmesh_int_t node_0 = s_t2n_diag[3*0+0];
    wmesh_int_t node_1 = s_t2n_diag[3*0+1];
    wmesh_int_t node_2 = s_t2n_diag[3*0+2];
    
    for (wmesh_int_t i=1;i<degree_;++i)
      {
	const auto ri = r_v_[i];
	for (wmesh_int_t j=1;j<degree_-i;++j)
	  {
	    const auto rj	= r_v_[j*r_inc_];
	    const auto rl	= r_v_[(degree_ - i - j)*r_inc_];
	    
	    T l1  = ( two * (one + ri) - (rj + rl) ) / six;
	    T l2  = ( two * (one + rj) - (ri + rl) ) / six;
	    T l0  = one - l1 - l2;
	      
	    wmesh_int_t s1  = i;
	    wmesh_int_t s2  = j;
	    wmesh_int_t s0  = degree_-s1 - s2;
	    for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
	      ijk[idim] = level +
		s0 * s_rst_v[s_rst_ld*node_0 + idim] +
		s1 * s_rst_v[s_rst_ld*node_1 + idim] +
		s2 * s_rst_v[s_rst_ld*node_2 + idim];
	      
	    wmesh_int_t idx = perm_[flat3d(r_n_,ijk)]-1;	  
	    for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
	      c_v_[idx * c_ld_ + idim] =
		shift_ + 
		( l0 * s_rst_v[s_rst_ld*node_0 + idim] +
		  l1 * s_rst_v[s_rst_ld*node_1 + idim] +
		  l2 * s_rst_v[s_rst_ld*node_2 + idim] ) * scal_;
	  }
      }
  }

  
  for (wmesh_int_t itet=0;itet<2;++itet)
    {
      
      wmesh_int_t node_0 = s_tet2n[4*itet+0];
      wmesh_int_t node_1 = s_tet2n[4*itet+1];
      wmesh_int_t node_2 = s_tet2n[4*itet+2];
      wmesh_int_t node_3 = s_tet2n[4*itet+3];      
      for (wmesh_int_t i=1;i<degree_;++i)
	{
	  const auto ri = r_v_[i*r_inc_];
	  for (wmesh_int_t j=1;j<degree_-i;++j)
	    {
	      const auto rj	= r_v_[j*r_inc_];	      
	      for (wmesh_int_t k=1;k<degree_-i-j;++k)
		{
		  const auto rk	= r_v_[k*r_inc_];
		  
		  const auto rl	= r_v_[(degree_ - i - j - k)*r_inc_];
		  
		  T l1  = ( two + three * ri - (rj + rk + rl) ) / eight;
		  T l2  = ( two + three * rj - (ri + rk + rl) ) / eight;
		  T l3  = ( two + three * rk - (ri + rj + rl) ) / eight;
		  T l0  = one - l1 - l2 - l3;

		  wmesh_int_t s1  = i;
		  wmesh_int_t s2  = j;
		  wmesh_int_t s3  = k;
		  wmesh_int_t s0  = degree_ - s1 - s2 - s3;
		  for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
		    ijk[idim] = level +
		      s0 * s_rst_v[s_rst_ld*node_0 + idim] +
		      s1 * s_rst_v[s_rst_ld*node_1 + idim] +
		      s2 * s_rst_v[s_rst_ld*node_2 + idim] +
		      s3 * s_rst_v[s_rst_ld*node_3 + idim];

		  wmesh_int_t idx = perm_[flat3d(r_n_,ijk)]-1;	  
		  for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
		    {
		      c_v_[idx * c_ld_ + idim] =
			shift_ + 
			( l0 * s_rst_v[s_rst_ld*node_0 + idim] +
			  l1 * s_rst_v[s_rst_ld*node_1 + idim] +
			  l2 * s_rst_v[s_rst_ld*node_2 + idim] +
			  l3 * s_rst_v[s_rst_ld*node_3 + idim]) * scal_;
		      std:: cout << " " << idim <<  " " << c_v_[idx * c_ld_ + idim] << std::endl;

		    }

		}
	    }
	}
      
    }
      

  
  
  return WMESH_STATUS_SUCCESS;
}

template<typename T>
wmesh_status_t
bms_nodes(wmesh_int_t 		element_,
	  wmesh_int_t 		family_,
	  wmesh_int_t		degree_,

	  wmesh_int_t		b_storage_,
	  wmesh_int_t		b_m_,
	  wmesh_int_t		b_n_,
	  const_wmesh_int_p 	b_v_,
	  wmesh_int_t		b_ld_,

	  wmesh_int_t		c_storage_,
	  wmesh_int_t		c_m_,
	  wmesh_int_t		c_n_,
	  T*__restrict__ 	c_v_,
	  wmesh_int_t		c_ld_,

	  wmesh_int_t		work_n_,
	  T* __restrict__ 	work_)  
{

  
#ifndef NDEBUG
  std::cout << "// WMESH : bms_nodes : element = " << element_ << std::endl;
  std::cout << "// WMESH : bms_nodes : family  = " << family_ << std::endl;
  std::cout << "// WMESH : bms_nodes : degree  = " << degree_ << std::endl;
#endif
  
  switch(family_)
    {
    case WMESH_NODES_FAMILY_LAGRANGE:
      {
	
#if 0
	status = bms_ordering(element_,
			      degree_,
			      c_storage_,
			      c_m_,
			      c_n_,
			      c_v_,
			      c_ld_);
#endif
	
	if (degree_ > 0)
	  {
	    for (wmesh_int_t j=0;j<c_n_;++j)
	      {
		for (wmesh_int_t i=0;i<c_m_;++i)
		  {
		    c_v_[c_ld_*j+i] = ((T)b_v_[b_ld_*j+i]) / ((T)degree_);
		  }
	      }
	  }
	
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_NODES_FAMILY_GAUSSLOBATTO:
      {
	wmesh_status_t status;

	if (element_ == WMESH_ELEMENT_EDGE)
	  {
	    T w[64], work[1024];
	    c_v_[0]  = -1.0;
	    status =  wmesh_nodes_legendre(degree_-1,
					   c_storage_,
					   c_v_ + 1 * c_ld_,
					   c_ld_, 
					   w,
					   1,
					   1024,
					   work);
	    c_v_[degree_*c_ld_]  = 1.0;	    
	    WMESH_STATUS_CHECK(status);
	    return WMESH_STATUS_SUCCESS;
	  }
	
	T * p1d 	= 	work_;
	wmesh_int_t n1d = 	(degree_+1);
	work_n_  	-= 	n1d;
	
	status  = bms_nodes(WMESH_ELEMENT_EDGE,
			    WMESH_NODES_FAMILY_GAUSSLOBATTO,
			    degree_,
			    b_storage_,
			    b_m_,
			    b_n_,
			    b_v_,
			    b_ld_,
			    WMESH_STORAGE_INTERLEAVE,
			    1,
			    n1d,
			    p1d,
			    1,
			    work_n_,
			    work_ + n1d);

	for (wmesh_int_t i=0;i<n1d;++i)
	  {
	    std::cout << "# " << p1d[i] << std::endl;
	  }
	
#if 0
	std::cout << "---- " << std::endl;
	for (wmesh_int_t j=0;j<b_n_;++j)
	  {
	    for (wmesh_int_t i=0;i<b_m_;++i)
	      {
		std::cout << " " << b_v_[b_ld_*j+i];
	      }
	    std::cout << std::endl;
	  }
       	std::cout << "---- " << std::endl;
#endif
	if (element_ == WMESH_ELEMENT_PYRAMID)
	  {
	    wmesh_int_p perm = (wmesh_int_p)malloc(n1d*n1d*n1d*sizeof(wmesh_int_t));
	    
	    status =  bms_ordering_flat(degree_,
					b_storage_,
					b_m_,
					b_n_,
					b_v_,
					b_ld_,
					perm);
	    WMESH_STATUS_CHECK(status);

	    wmesh_int_t level = 0;
	    T shift = 0.0;
	    T scal  = 1.0;
	    status  = bms_nodes_pyramid(degree_,			      			      
					c_storage_,
					c_m_,
				       c_n_,
				       c_v_,
				       c_ld_,				       
				       level,
				       shift,
				       scal,
				       perm,
				       work_n_,
				       work_);
	    WMESH_STATUS_CHECK(status);
	  }
	else
	  {

	    wmesh_int_t iwork_n = n1d*n1d*n1d;
	    wmesh_int_p iwork = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*iwork_n);
	    status = bms_transform(element_,
				   n1d,
				   p1d,
				   1,
			       
				   b_storage_,
				   b_m_,
				   b_n_,
				   b_v_,
				   b_ld_,
			       
				   c_storage_,
				   c_m_,
				   c_n_,
				   c_v_,					 
				   c_ld_,
			       
				   iwork_n,
				   iwork);	
	    free(iwork);
	    WMESH_STATUS_CHECK(status);

	  }
#if 0
	for (wmesh_int_t j=0;j<c_n_;++j)
	  {
	    for (wmesh_int_t i=0;i<c_m_;++i)
	      {
		std::cout << " " << c_v_[c_ld_*j+i];
	      }
	    std::cout << std::endl;
	  }
#endif
	return WMESH_STATUS_SUCCESS;
      }
    case WMESH_NODES_FAMILY_BEZIER:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);	
      }

    }
  
  return WMESH_STATUS_INVALID_ENUM;  
}

extern "C"
{

  wmesh_status_t
  bms_nodes_buffer_size(wmesh_int_t 	element_,
			wmesh_int_t 	family_,
			wmesh_int_t	degree_,
			wmesh_int_p	work_n_)
  {
    WMESH_CHECK_POINTER(work_n_);
    work_n_[0] = (degree_+1) * (degree_+1) * (degree_+1);
    return WMESH_STATUS_SUCCESS;
  }
  
  wmesh_status_t
  bms_snodes(wmesh_int_t 	element_,
	     wmesh_int_t 	family_,
	     wmesh_int_t	degree_,
	     
	     wmesh_int_t	b_storage_,
	     wmesh_int_t	b_m_,
	     wmesh_int_t	b_n_,
	     const_wmesh_int_p 	b_v_,
	     wmesh_int_t	b_ld_,

	     wmesh_int_t	c_storage_,
	     wmesh_int_t	c_m_,
	     wmesh_int_t	c_n_,
	     float* 		c_v_,
	     wmesh_int_t	c_ld_,
	     
	     wmesh_int_t	work_n_,
	     float* 		work_)
  {
    return bms_nodes(element_,
		     family_,
		     degree_,
		     b_storage_,
		     b_m_,
		     b_n_,
		     b_v_,
		     b_ld_,
		     c_storage_,
		     c_m_,
		     c_n_,
		     c_v_,
		     c_ld_,
		     work_n_,
		     work_);    
  }

  wmesh_status_t
  bms_dnodes(wmesh_int_t 	element_,
	     wmesh_int_t 	family_,
	     wmesh_int_t	degree_,
	     wmesh_int_t	b_storage_,
	     wmesh_int_t	b_m_,
	     wmesh_int_t	b_n_,
	     const_wmesh_int_p 	b_v_,
	     wmesh_int_t	b_ld_,
	     wmesh_int_t	c_storage_,
	     wmesh_int_t	c_m_,
	     wmesh_int_t	c_n_,
	     double* 		c_v_,
	     wmesh_int_t	c_ld_,
	     wmesh_int_t	work_n_,
	     double* 		work_)
  {    
    return bms_nodes(element_,
		     family_,
		     degree_,
		     b_storage_,
		     b_m_,
		     b_n_,
		     b_v_,
		     b_ld_,
		     c_storage_,
		     c_m_,
		     c_n_,
		     c_v_,
		     c_ld_,
		     work_n_,
		     work_);    
  }

}
