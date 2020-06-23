#include "bms.h"
#include "wmesh-status.h"
#include "wmesh-enums.h"
#include <iostream>

static inline wmesh_int_t bms_ordering_flat3d(wmesh_int_t n_,wmesh_int_t l[])
{
  return l[0] + l[1]*n_ + l[2]*n_*n_;
}


static inline wmesh_int_t bms_ordering_flat2d(wmesh_int_t n_,wmesh_int_t l[])
{
  return l[0] + l[1]*n_;
}


template<typename T>
wmesh_status_t
bms_transform_pyramid(wmesh_int_t 		r_n_,
		      const T*__restrict__	r_v_,
		      wmesh_int_t 		r_inc_,
		  
		      wmesh_int_t		c_storage_,
		      wmesh_int_t		c_m_,
		      wmesh_int_t		c_n_,
		      T*__restrict__ 		c_v_,
		      wmesh_int_t		c_ld_,

		      const_wmesh_int_p		p_v_,
		      wmesh_int_t 		p_inc_)  
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

  static constexpr const T one(1);
  static constexpr const T two(2);
  static constexpr const T three(3);
  static constexpr const T six(6);
  static constexpr const T eight(8);  

		  

  wmesh_int_t ijk[3];
  wmesh_int_t degree = r_n_ - 1;
  
  //
  // nodes
  //
  for (wmesh_int_t j=0;j<s_rst_n;++j)
    {
      
      for (wmesh_int_t i=0;i<s_rst_m;++i)
	{
	  ijk[i] = s_rst_v[s_rst_ld*j+i] * degree;
	}
	  
      wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat3d(r_n_,ijk)]-1;
      for (wmesh_int_t i=0;i<s_rst_m;++i)
	{
	  c_v_[c_ld_*idx+i] = s_rst_v[s_rst_ld*j+i];
	}
    }
      
      
  //
  // edges
  //
  for (wmesh_int_t ie=0;ie<8;++ie)
    {
      wmesh_int_t node_0 = s_e2n[ie*2+0];
      wmesh_int_t node_1 = s_e2n[ie*2+1];
      // std::cout << "node " << node_0  << " " << node_1 << std::endl;
      for (wmesh_int_t i=1;i<degree;++i)
	{
	  T l0 	= one - (r_v_[ r_inc_* i ] + one) / two;
	  T l1 	= (r_v_[ r_inc_* i ] + one) / two;
	      
	  wmesh_int_t s0 = degree - i;
	  wmesh_int_t s1 = i;
	      
	  for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
	    {
	      ijk[idim] = s0 * s_rst_v[s_rst_ld*node_0 + idim] + s1 * s_rst_v[s_rst_ld*node_1 + idim];
	    }
	      
	  wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat3d(r_n_,ijk)]-1;
	      
	  for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
	    {
	      c_v_[c_ld_*idx+idim] = (l0 * s_rst_v[s_rst_ld*node_0 + idim] + l1 * s_rst_v[s_rst_ld*node_1 + idim]);
	      //	      std::cout << "idim  " << idim << " = " << c_v_[c_ld_*idx+idim] << std::endl;
	    }	  
	}
    }

  //
  // boundary triangles
  //
  for (wmesh_int_t it=0;it<4;++it)
    {
      wmesh_int_t node_0 = s_t2n[it*3+0];
      wmesh_int_t node_1 = s_t2n[it*3+1];
      wmesh_int_t node_2 = s_t2n[it*3+2];
	  
      for (wmesh_int_t i=1;i<degree;++i)
	{
	  const auto ri = r_v_[i];
	  for (wmesh_int_t j=1;j<degree-i;++j)
	    {
	      const auto rj	= r_v_[j*r_inc_];
	      const auto rl	= r_v_[(degree - i - j)*r_inc_];

	      T l1  = ( two * (one + ri) - (rj + rl) ) / six;
	      T l2  = ( two * (one + rj) - (ri + rl) ) / six;
	      T l0  = one - l1 - l2;

	      wmesh_int_t s1  = i;
	      wmesh_int_t s2  = j;
	      wmesh_int_t s0  = degree-s1 - s2;
	      for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
		ijk[idim] = 
		  s0 * s_rst_v[s_rst_ld*node_0 + idim] +
		  s1 * s_rst_v[s_rst_ld*node_1 + idim] +
		  s2 * s_rst_v[s_rst_ld*node_2 + idim];
	      
	      wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat3d(r_n_,ijk)]-1;	  
	      for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
		c_v_[idx * c_ld_ + idim] =
		  
		  ( l0 * s_rst_v[s_rst_ld*node_0 + idim] +
		    l1 * s_rst_v[s_rst_ld*node_1 + idim] +
		    l2 * s_rst_v[s_rst_ld*node_2 + idim] );
	    }
	}
    }
  

  for (wmesh_int_t iq=0;iq<1;++iq)
    {
      wmesh_int_t node_0 = s_q2n[iq*4+0];
      wmesh_int_t node_1 = s_q2n[iq*4+1];
      wmesh_int_t node_2 = s_q2n[iq*4+2];
      wmesh_int_t node_3 = s_q2n[iq*4+3];
	
      for (wmesh_int_t i=1;i<degree;++i)
	{
	  const auto ri = (r_v_[i*r_inc_]+one)/two;
	  for (wmesh_int_t j=1;j<degree;++j)
	    {
	      const auto rj = (r_v_[j*r_inc_]+one) / two;
	       
	      T l0  = (one - ri)*(one - rj);
	      T l1   = ri*(one - rj);
	      T l2  = ri*rj;
	      T l3  = (one - ri)*rj;
#if 0	     
	      wmesh_int_t s0  = (degree - i)*(degree - j);
	      wmesh_int_t s1  = i*(degree - j);
	      wmesh_int_t s2  = i * j ;
	      wmesh_int_t s3  = (degree - i) * j;
#endif

	      ijk[0] = j;
	      ijk[1] = i;
	      ijk[2] = 0;

#if 0
	      for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
		ijk[idim] = 
		  s0 * s_rst_v[s_rst_ld*node_0 + idim] +
		  s1 * s_rst_v[s_rst_ld*node_1 + idim] +
		  s2 * s_rst_v[s_rst_ld*node_2 + idim] + 
		  s3 * s_rst_v[s_rst_ld*node_3 + idim];
#endif	      
	      //	      std::cout << "ddddddddddddddddd " << ijk[0] << " " << ijk[1] << " " << ijk[2] << std::endl;
	      wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat3d(r_n_,ijk)]-1;
	      
	      for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
		c_v_[idx * c_ld_ + idim] =
		  (l0 * s_rst_v[s_rst_ld*node_0 + idim] +
			    l1 * s_rst_v[s_rst_ld*node_1 + idim] +
			    l2 * s_rst_v[s_rst_ld*node_2 + idim] + 
			    l3 * s_rst_v[s_rst_ld*node_3 + idim]);		
	    }
	}
    }	



  {      
    wmesh_int_t node_0 = s_t2n_diag[3*0+0];
    wmesh_int_t node_1 = s_t2n_diag[3*0+1];
    wmesh_int_t node_2 = s_t2n_diag[3*0+2];
    
    for (wmesh_int_t i=1;i<degree;++i)
      {
	const auto ri = r_v_[i];
	for (wmesh_int_t j=1;j<degree-i;++j)
	  {
	    const auto rj	= r_v_[j*r_inc_];
	    const auto rl	= r_v_[(degree - i - j)*r_inc_];
	    
	    T l1  = ( two * (one + ri) - (rj + rl) ) / six;
	    T l2  = ( two * (one + rj) - (ri + rl) ) / six;
	    T l0  = one - l1 - l2;
	      
	    wmesh_int_t s1  = i;
	    wmesh_int_t s2  = j;
	    wmesh_int_t s0  = degree-s1 - s2;
	    for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
	      ijk[idim] = 
		s0 * s_rst_v[s_rst_ld*node_0 + idim] +
		s1 * s_rst_v[s_rst_ld*node_1 + idim] +
		s2 * s_rst_v[s_rst_ld*node_2 + idim];
	      
	    wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat3d(r_n_,ijk)]-1;	  
	    for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
	      c_v_[idx * c_ld_ + idim] =
		
		( l0 * s_rst_v[s_rst_ld*node_0 + idim] +
		  l1 * s_rst_v[s_rst_ld*node_1 + idim] +
		  l2 * s_rst_v[s_rst_ld*node_2 + idim] );
	  }
      }
  }

  
  for (wmesh_int_t itet=0;itet<2;++itet)
    {
      
      wmesh_int_t node_0 = s_tet2n[4*itet+0];
      wmesh_int_t node_1 = s_tet2n[4*itet+1];
      wmesh_int_t node_2 = s_tet2n[4*itet+2];
      wmesh_int_t node_3 = s_tet2n[4*itet+3];      
      for (wmesh_int_t i=1;i<degree;++i)
	{
	  const auto ri = r_v_[i*r_inc_];
	  for (wmesh_int_t j=1;j<degree-i;++j)
	    {
	      const auto rj	= r_v_[j*r_inc_];	      
	      for (wmesh_int_t k=1;k<degree-i-j;++k)
		{
		  const auto rk	= r_v_[k*r_inc_];
		  
		  const auto rl	= r_v_[(degree - i - j - k)*r_inc_];
		  
		  T l1  = ( two + three * ri - (rj + rk + rl) ) / eight;
		  T l2  = ( two + three * rj - (ri + rk + rl) ) / eight;
		  T l3  = ( two + three * rk - (ri + rj + rl) ) / eight;
		  T l0  = one - l1 - l2 - l3;

		  wmesh_int_t s1  = i;
		  wmesh_int_t s2  = j;
		  wmesh_int_t s3  = k;
		  wmesh_int_t s0  = degree - s1 - s2 - s3;
		  for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
		    ijk[idim] = 
		      s0 * s_rst_v[s_rst_ld*node_0 + idim] +
		      s1 * s_rst_v[s_rst_ld*node_1 + idim] +
		      s2 * s_rst_v[s_rst_ld*node_2 + idim] +
		      s3 * s_rst_v[s_rst_ld*node_3 + idim];

		  wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat3d(r_n_,ijk)]-1;	  
		  for (wmesh_int_t idim=0;idim<s_rst_m;++idim)
		    {
		      c_v_[idx * c_ld_ + idim] =
			
			( l0 * s_rst_v[s_rst_ld*node_0 + idim] +
			  l1 * s_rst_v[s_rst_ld*node_1 + idim] +
			  l2 * s_rst_v[s_rst_ld*node_2 + idim] +
			  l3 * s_rst_v[s_rst_ld*node_3 + idim]);
		      //		      std:: cout << " " << idim <<  " " << c_v_[idx * c_ld_ + idim] << std::endl;

		    }

		}
	    }
	}
      
    }
  
  return WMESH_STATUS_SUCCESS;
}


template<typename T>
static wmesh_status_t bms_transform_tetrahedron(wmesh_int_t 		r_n_,
						const T*__restrict__ 	r_v_,
						wmesh_int_t 		r_inc_,
						wmesh_int_t		c_storage_,
						wmesh_int_t		c_m_,
						wmesh_int_t		c_n_,
						T*			c_v_,
						wmesh_int_t		c_ld_,
						const_wmesh_int_p	p_v_,
						wmesh_int_t		p_inc_)
{
  static constexpr const T s_T0(0);
  static constexpr const T s_T1(1);
  static constexpr const T s_T2(2);
  static constexpr const T s_T3(3);


  T d_ = 2;
  wmesh_int_t n_ 	= r_n_;
  wmesh_int_t m 	= r_n_-1;
  wmesh_int_t ijk[3];
  for (wmesh_int_t i=0;i<n_;++i)
    {
      for (wmesh_int_t j=0;j<n_ - i;++j)
	{
	  wmesh_int_t l = m - i - j;
	  ijk[0] = i;
	  ijk[1] = j;
	  ijk[2] = 0;
	  wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat3d(n_,ijk)]-1;
	  if(c_storage_ == WMESH_STORAGE_INTERLEAVE)
	    {
	      
	      //  c_v_[c_ld_ * idx + 0] = ( two * (one + ri) - (rj + rl) ) / six;
	      //  c_v_[c_ld_ * idx + 1] = ( two * (one + rj) - (ri + rl) ) / six;
			
	      c_v_[c_ld_*idx+0] = (d_*s_T1+s_T2*r_v_[i*r_inc_]-r_v_[j*r_inc_]-r_v_[l*r_inc_]) / s_T3 / d_;
	      c_v_[c_ld_*idx+1] = (d_*s_T1+s_T2*r_v_[j*r_inc_]-r_v_[i*r_inc_]-r_v_[l*r_inc_]) / s_T3 / d_;
	      c_v_[c_ld_*idx+2] = s_T0;
	    }
	  else
	    {
	      c_v_[c_ld_*0+idx] = (d_*s_T1+s_T2*r_v_[i*r_inc_]-r_v_[j*r_inc_]-r_v_[l*r_inc_]) / s_T3 / d_;
	      c_v_[c_ld_*1+idx] = (d_*s_T1+s_T2*r_v_[j*r_inc_]-r_v_[i*r_inc_]-r_v_[l*r_inc_]) / s_T3 / d_;
	      c_v_[c_ld_*2+idx] = s_T0;
	    }
	}
    }
  for (wmesh_int_t j=0;j < m;++j)
    {
      for (wmesh_int_t k=1;k < n_ - j;++k)
	{
	  ijk[0] = 0;
	  ijk[1] = j;
	  ijk[2] = k;
	  wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat3d(n_,ijk)]-1;

	  wmesh_int_t l = m - j - k;
	  if(c_storage_ == WMESH_STORAGE_INTERLEAVE)
	    {
	      c_v_[c_ld_*idx+0] = s_T0;
	      c_v_[c_ld_*idx+1] = (d_*s_T1+s_T2*r_v_[r_inc_*j]-r_v_[r_inc_*k]-r_v_[r_inc_*l]) / s_T3 / d_;
	      c_v_[c_ld_*idx+2] = (d_*s_T1+s_T2*r_v_[r_inc_*k]-r_v_[r_inc_*j]-r_v_[r_inc_*l]) / s_T3 / d_;
	    }
	  else
	    {
	      c_v_[c_ld_*0+idx] = s_T0;
	      c_v_[c_ld_*1+idx] = (d_*s_T1+s_T2*r_v_[r_inc_*j]-r_v_[r_inc_*k]-r_v_[r_inc_*l]) / s_T3 / d_;
	      c_v_[c_ld_*2+idx] = (d_*s_T1+s_T2*r_v_[r_inc_*k]-r_v_[r_inc_*j]-r_v_[r_inc_*l]) / s_T3 / d_;
	    }
	}
    }

  
  for (wmesh_int_t i=1;i< m;++i)
    {
      for (wmesh_int_t k=1;k < n_ - i;++k)
	{
	  ijk[0] = i;
	  ijk[1] = 0;
	  ijk[2] = k;
	  wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat3d(n_,ijk)]-1;


	  wmesh_int_t l = m - i - k;
	  if(c_storage_ == WMESH_STORAGE_INTERLEAVE)
	    {
	      c_v_[c_ld_*idx+0] = (d_*s_T1+s_T2*r_v_[r_inc_*i]-r_v_[r_inc_*k]-r_v_[r_inc_*l]) / s_T3 / d_;
	      c_v_[c_ld_*idx+1] = s_T0;
	      c_v_[c_ld_*idx+2] = (d_*s_T1+s_T2*r_v_[r_inc_*k]-r_v_[r_inc_*i]-r_v_[r_inc_*l]) / s_T3 / d_;
	    }
	  else
	    {
	      c_v_[c_ld_*0+idx] = (d_*s_T1+s_T2*r_v_[r_inc_*i]-r_v_[r_inc_*k]-r_v_[r_inc_*l]) / s_T3 / d_;
	      c_v_[c_ld_*1+idx] = s_T0;
	      c_v_[c_ld_*2+idx] = (d_*s_T1+s_T2*r_v_[r_inc_*k]-r_v_[r_inc_*i]-r_v_[r_inc_*l]) / s_T3 / d_;
	    }
	}
    }
  
  for (wmesh_int_t i = 1; i < m;++i)
    {
      for (wmesh_int_t j = 1;j < n_ - i;++j)
	{
	  ijk[0] = i;
	  ijk[1] = j;
	  ijk[2] = m - i - j;
	  wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat3d(n_,ijk)]-1;
	  
	  wmesh_int_t l = m - i - j;	  
	  T xi = (s_T1*d_+s_T2*r_v_[r_inc_*i]-r_v_[r_inc_*j]-r_v_[r_inc_*l]) / s_T3;
	  T eta = (s_T1*d_ +s_T2*r_v_[r_inc_*j]-r_v_[r_inc_*i]-r_v_[r_inc_*l]) / s_T3;
	  if(c_storage_ == WMESH_STORAGE_INTERLEAVE)
	    {
	      c_v_[c_ld_*idx+0] = xi / d_;
	      c_v_[c_ld_*idx+1] = eta / d_;
	      c_v_[c_ld_*idx+2] = (s_T1*d_-xi-eta) / d_;
	    }
	  else
	    {
	      c_v_[c_ld_*0+idx] = xi / d_;
	      c_v_[c_ld_*0+idx] = eta / d_;
	      c_v_[c_ld_*0+idx] = (s_T1*d_-xi-eta) / d_;
	    }
	}
    }
  
  static constexpr const T s_T4(4.0);
  for (wmesh_int_t i=1;i< m;++i)
    {
      for (wmesh_int_t j=1;j < m - i;++j)
	{
	  for (wmesh_int_t k=1;k < m - i - j;++k)
	    {
	      
	      ijk[0] = i;
	      ijk[1] = j;
	      ijk[2] = k;
	      
	      wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat3d(n_,ijk)]-1;
	      wmesh_int_t l = m - i - j - k;
	      if(c_storage_ == WMESH_STORAGE_INTERLEAVE)
		{
		  c_v_[c_ld_*idx+0] = (d_*s_T1+s_T3*r_v_[r_inc_*i]-r_v_[r_inc_*j]-r_v_[r_inc_*k]-r_v_[r_inc_*l]) / s_T4 / d_;
		  c_v_[c_ld_*idx+1] = (d_*s_T1+s_T3*r_v_[r_inc_*j]-r_v_[r_inc_*i]-r_v_[r_inc_*k]-r_v_[r_inc_*l]) / s_T4 / d_;
		  c_v_[c_ld_*idx+2] = (d_*s_T1+s_T3*r_v_[r_inc_*k]-r_v_[r_inc_*i]-r_v_[r_inc_*j]-r_v_[r_inc_*l]) / s_T4 / d_;
		}
	      else
		{
		  c_v_[c_ld_*0+idx] = (d_*s_T1+s_T3*r_v_[r_inc_*i]-r_v_[r_inc_*j]-r_v_[r_inc_*k]-r_v_[r_inc_*l]) / s_T4 / d_;
		  c_v_[c_ld_*1+idx] = (d_*s_T1+s_T3*r_v_[r_inc_*j]-r_v_[r_inc_*i]-r_v_[r_inc_*k]-r_v_[r_inc_*l]) / s_T4 / d_;
		  c_v_[c_ld_*2+idx] = (d_*s_T1+s_T3*r_v_[r_inc_*k]-r_v_[r_inc_*i]-r_v_[r_inc_*j]-r_v_[r_inc_*l]) / s_T4 / d_;
		}
	    }
	}
    }

  return WMESH_STATUS_SUCCESS;  
}


template<typename T>
wmesh_status_t bms_transform(wmesh_int_t 		element_,
			     wmesh_int_t 		r_n_,					 
			     const T*__restrict__	r_v_,
			     wmesh_int_t		r_inc_,
			     
			     wmesh_int_t		c_storage_,
			     wmesh_int_t		c_m_,
			     wmesh_int_t		c_n_,
			     T * __restrict__		c_v_,					 
			     wmesh_int_t		c_ld_,
			     
			     const_wmesh_int_p		p_v_,
			     wmesh_int_t		p_inc_) 
{
  
  static constexpr const T one(1.0);
  static constexpr const T two(2.0);
  static constexpr const T six(6.0);

  WMESH_CHECK_POSITIVE(r_n_);
  WMESH_CHECK_POINTER(r_v_);
  WMESH_CHECK(r_inc_ >= 1);
 
  wmesh_int_t c_dim 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_m_ : c_n_;
  
  WMESH_CHECK( (c_dim == 2) || (c_dim ==3) );   
  WMESH_CHECK(c_ld_ >= c_m_);
  
  wmesh_int_t degree = r_n_ - 1;
  wmesh_int_t ijk[3];  
  if (element_ == WMESH_ELEMENT_TRIANGLE)
    {
      switch(c_storage_)
	{
	case WMESH_STORAGE_INTERLEAVE:
	  {	
	    for (wmesh_int_t i=0;i<r_n_;++i)
	      {
		ijk[0] = i;
		const auto ri = r_v_[i];
		for (wmesh_int_t j=0;j<r_n_-i;++j)
		  {
		    ijk[1] = j;
		    wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat2d(r_n_,ijk)] - 1;
		    
		    const auto rj	= r_v_[j*r_inc_];
		    const auto rl	= r_v_[(degree - i - j)*r_inc_];	  
		    c_v_[c_ld_ * idx + 0] = ( two * (one + ri) - (rj + rl) ) / six;
		    c_v_[c_ld_ * idx + 1] = ( two * (one + rj) - (ri + rl) ) / six;
		  }
	      }
	    break;
	  }
      
	case WMESH_STORAGE_BLOCK:
	  {
	    for (wmesh_int_t i=0;i<=degree;++i)
	      {
		ijk[0] = i;
		const auto ri = r_v_[i];
		for (wmesh_int_t j=0;j<=degree-i;++j)
		  {
		    ijk[1] = j;
		    WMESH_CHECK(p_v_[p_inc_*bms_ordering_flat2d(r_n_,ijk)] > 0);
		    wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat2d(r_n_,ijk)] - 1;
		    
		    const auto rj	= r_v_[j*r_inc_];
		    const auto rl	= r_v_[(degree-i-j)*r_inc_];	  
		    c_v_[c_ld_ * 0 + idx] = ( two * (one + ri) - (rj + rl) ) / six;
		    c_v_[c_ld_ * 1 + idx] = ( two * (one + rj) - (ri + rl) ) / six;
		  }
	      }
	
	    break;
	  }
      
	default:
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
	  }
      
	}
    }
  else if (element_ == WMESH_ELEMENT_QUADRILATERAL || element_ == WMESH_ELEMENT_HEXAHEDRON)
    {
      switch(c_storage_)
	{
	case WMESH_STORAGE_INTERLEAVE:
	  {	
	    if (2 == c_dim)
	      {
		
		for (wmesh_int_t i=0;i<=degree;++i)
		  {
		    ijk[0] = i;
		    const auto ri = (r_v_[i]+one)/two;
		    for (wmesh_int_t j=0;j<=degree;++j)
		      {
			const auto rj	= (r_v_[j*r_inc_]+one)/two;
			ijk[1] = j;
#ifndef NDEBUG			    
			WMESH_CHECK(p_v_[p_inc_*bms_ordering_flat2d(r_n_,ijk)] > 0);
#endif
			wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat2d(r_n_,ijk)] - 1;
			c_v_[c_ld_ * idx + 0] = ri;
			c_v_[c_ld_ * idx + 1] = rj;
		      }
		  }

	      }
	    else
	      {
#ifndef NDEBUG			    
		WMESH_CHECK(3 == c_dim);
#endif
		//		std::cout << "hello degree="  << degree << std::endl;
		for (wmesh_int_t i=0;i<r_n_;++i)
		  {
		    ijk[0] = i;
		    const auto ri = (r_v_[i*r_inc_]+one)/two;
		    for (wmesh_int_t j=0;j<r_n_;++j)
		      {
			ijk[1] = j;
			const auto rj = (r_v_[j*r_inc_]+one)/two;
			for (wmesh_int_t k=0;k<r_n_;++k)
			  {
			    ijk[2] = k;
#ifndef NDEBUG			    
			    WMESH_CHECK(p_v_[p_inc_*bms_ordering_flat3d(r_n_,ijk)] > 0);
#endif
			    wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat3d(r_n_,ijk)] - 1;
			    // std::cout << " >  idx = " << idx + 1 << ", c_ld_ = " << c_ld_ << ", flat = "<< bms_ordering_flat3d(r_n_,ijk) << std::endl;
			    const auto rk = (r_v_[k*r_inc_]+one)/two;			
			    c_v_[c_ld_*idx+0] = ri;
			    c_v_[c_ld_*idx+1] = rj;
			    c_v_[c_ld_*idx+2] = rk;
			    //			    std::cout << "coo " << ri << " " << rj << " " << rk << ", ijk " << ijk[0] << " " << ijk[1] << " " << ijk[2] << std::endl;
			  }
		      }
		  }
	      }
	    
	    break;
	  }
	  
	case WMESH_STORAGE_BLOCK:
	  {
	    if (2 == c_dim)
	      {	    
		for (wmesh_int_t i=0;i<=degree;++i)
		  {
		    ijk[0] = i;
		    const auto ri = (r_v_[i]+one)/two;
		    for (wmesh_int_t j=0;j<=degree;++j)
		      {
			const auto rj	= (r_v_[j*r_inc_]+one)/two;
			ijk[1] = j;
#ifndef NDEBUG			    
			WMESH_CHECK(p_v_[p_inc_*bms_ordering_flat2d(r_n_,ijk)] > 0);
#endif
			wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat2d(r_n_,ijk)] - 1;
			c_v_[c_ld_ * 0 + idx] = ri;
			c_v_[c_ld_ * 1 + idx] = rj;
		      }
		  }

	      }
	    else
	      {
#ifndef NDEBUG			    		
		WMESH_CHECK(3 == c_dim);
#endif
		for (wmesh_int_t i=0;i<=degree;++i)
		  {
		    ijk[0] = i;
		    const auto ri = (r_v_[i*r_inc_]+one)/two;
		    for (wmesh_int_t j=0;j<=degree;++j)
		      {
			ijk[1] = j;
			const auto rj = (r_v_[j*r_inc_]+one)/two;
			for (wmesh_int_t k=0;k<=degree;++k)
			  {
			    ijk[2] = k;
#ifndef NDEBUG			    		
			    WMESH_CHECK(p_v_[p_inc_*bms_ordering_flat3d(r_n_,ijk)] > 0);
#endif
			    wmesh_int_t idx = p_v_[p_inc_*bms_ordering_flat3d(r_n_,ijk)] - 1;			
			    const auto rk = (r_v_[k*r_inc_]+one)/two;			
			    c_v_[c_ld_*0+idx] = ri;
			    c_v_[c_ld_*1+idx] = rj;
			    c_v_[c_ld_*2+idx] = rk;
			  }
		      }
		  }
	      }
	    break;
	  }
      
	default:
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
	  }
      
	}

      
    }
  else if (element_ == WMESH_ELEMENT_WEDGE)
    {
      switch(c_storage_)
	{
	case WMESH_STORAGE_INTERLEAVE:
	  {	
#ifndef NDEBUG			    		
	    WMESH_CHECK(3 == c_dim);
#endif
	    for (wmesh_int_t i=0;i<=degree;++i)
	      {
		ijk[0] = i;
		const auto ri = r_v_[i];
		for (wmesh_int_t j=0;j<=degree-i;++j)
		  {
		    ijk[1] = j;
		    const auto rj = r_v_[j*r_inc_];
		    for (wmesh_int_t k=0;k<=degree;++k)
		      {
			ijk[2] = k;
			const auto rk		= r_v_[k*r_inc_];	  
#ifndef NDEBUG			    		
			WMESH_CHECK(p_v_[p_inc_*bms_ordering_flat3d(r_n_,ijk)] > 0);
#endif
			wmesh_int_t idx 	= p_v_[p_inc_*bms_ordering_flat3d(r_n_,ijk)] - 1;		    
			
			const auto rl		= r_v_[(degree-i-j)*r_inc_];	  
			c_v_[c_ld_ * idx + 0] 	= ( two * (one + ri) - (rj + rl) ) / six;
			c_v_[c_ld_ * idx + 1] 	= ( two * (one + rj) - (ri + rl) ) / six;
			c_v_[c_ld_ * idx + 2] 	= (rk + one)/two;
		      }
		  }
	      }
	    break;
	  }
	  
	case WMESH_STORAGE_BLOCK:
	  {
	    for (wmesh_int_t i=0;i<=degree;++i)
	      {
		ijk[0] = i;
		const auto ri = r_v_[i];
		for (wmesh_int_t j=0;j<=degree-i;++j)
		  {
		    ijk[1] = j;
		    const auto rj = r_v_[j*r_inc_];
		    for (wmesh_int_t k=0;k<=degree;++k)
		      {
			ijk[2] = k;
			const auto rk		= r_v_[k*r_inc_];
#ifndef NDEBUG			    		
			WMESH_CHECK(p_v_[p_inc_*bms_ordering_flat3d(r_n_,ijk)] > 0);
#endif
			wmesh_int_t idx 	= p_v_[p_inc_*bms_ordering_flat3d(r_n_,ijk)] - 1;		    

			const auto rl		= r_v_[(degree-i-j)*r_inc_];	  
			c_v_[c_ld_ * 0 + idx] 	= ( two * (one + ri) - (rj + rl) ) / six;
			c_v_[c_ld_ * 1 + idx] 	= ( two * (one + rj) - (ri + rl) ) / six;
			c_v_[c_ld_ * 2 + idx] 	= (rk + one)/two;
		      }
		  }
	      }
	    break;
	  }
      
	default:
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
	  }
      
	}

    }
  else if (element_ == WMESH_ELEMENT_PYRAMID)
    {
      wmesh_status_t status =  bms_transform_pyramid(r_n_,
						     r_v_,
						     r_inc_,
						     c_storage_,
						     c_m_,
						     c_n_,
						     c_v_,
						     c_ld_,
						     p_v_,
						     p_inc_);
      WMESH_STATUS_CHECK(status);      
    }
  else if (element_ == WMESH_ELEMENT_TETRAHEDRON)
    {
      wmesh_status_t status =  bms_transform_tetrahedron(r_n_,
							 r_v_,
							 r_inc_,
							 c_storage_,
							 c_m_,
							 c_n_,
							 c_v_,
							 c_ld_,
							 p_v_,
							 p_inc_);
      WMESH_STATUS_CHECK(status);      
    }
  return WMESH_STATUS_SUCCESS;
};


template
wmesh_status_t bms_transform<float>(wmesh_int_t 		element_,
				    wmesh_int_t 		r_n_,					 
				    const float*__restrict__	r_,
				    wmesh_int_t			r_inc_,
				    
				    wmesh_int_t			c_storage_,
				    wmesh_int_t			c_m_,
				    wmesh_int_t			c_n_,
				    float * __restrict__	c_v_,					 
				    wmesh_int_t			c_ld_,
				    
				    const_wmesh_int_p 		p_v_,
				    wmesh_int_t			p_inc_) ;

template
wmesh_status_t bms_transform<double>(wmesh_int_t 		element_,
				     wmesh_int_t 		r_n_,					 
				    const double*__restrict__	r_,
				    wmesh_int_t			r_inc_,
				    
				    wmesh_int_t			c_storage_,
				    wmesh_int_t			c_m_,
				    wmesh_int_t			c_n_,
				    double * __restrict__	c_v_,					 
				    wmesh_int_t			c_ld_,
				    
				    const_wmesh_int_p 		p_v_,
				    wmesh_int_t			p_inc_);


extern "C"
{
  wmesh_status_t bms_ordering_flat(wmesh_int_t		degree_,
				   wmesh_int_t		b_storage_,
				   wmesh_int_t		b_m_,
				   wmesh_int_t		b_n_,
				   const_wmesh_int_p	b_v_,
				   wmesh_int_t		b_ld_,
				   wmesh_int_p		p_v_,
				   wmesh_int_t		p_inc_)
  {
    wmesh_int_t r_n_ 		= degree_ + 1;
    wmesh_int_t b_dim 		= (b_storage_ == WMESH_STORAGE_INTERLEAVE) ? b_m_ : b_n_;
    wmesh_int_t b_count 	= (b_storage_ == WMESH_STORAGE_INTERLEAVE) ? b_n_ : b_m_;  
    WMESH_CHECK( (b_dim == 2) || (b_dim ==3) );    
    WMESH_CHECK(b_ld_ >= b_m_);
  
    for (wmesh_int_t i=0;i<b_count;++i)
      {
	p_v_[p_inc_*i] = 0;
      }
  
    wmesh_int_t ijk[3];
    switch(b_storage_)
      {
      case WMESH_STORAGE_INTERLEAVE:
	{
	  for (wmesh_int_t j=0;j<b_count;++j)
	    {
	      for (wmesh_int_t i=0;i<b_dim;++i)
		{
		  ijk[i] = b_v_[b_ld_ * j + i];
		}

	      wmesh_int_t flat  = (b_dim==2) ? bms_ordering_flat2d(r_n_,ijk) : bms_ordering_flat3d(r_n_,ijk);
	      p_v_[p_inc_*flat] = j + 1;
	    }
	  break;
	}
      case WMESH_STORAGE_BLOCK:
	{
	  for (wmesh_int_t j=0;j<b_count;++j)
	    {
	      for (wmesh_int_t i=0;i<b_dim;++i)
		{
		  ijk[i] = b_v_[b_ld_ * i + j];
		}
	      wmesh_int_t flat  = (b_dim==2) ? bms_ordering_flat2d(r_n_,ijk) : bms_ordering_flat3d(r_n_,ijk);
	      p_v_[p_inc_*flat] = j + 1;
	    }
	  break;
	}
      default:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
	}
      }
    return WMESH_STATUS_SUCCESS;
  }

}
