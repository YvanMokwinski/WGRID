#include "bms.h"
#include "wmesh-status.h"
#include "wmesh-enums.h"
#include <iostream>

static inline wmesh_int_t flat3d(wmesh_int_t n_,wmesh_int_t l[])
{
  return l[0] + l[1]*n_ + l[2]*n_*n_;
}


static inline wmesh_int_t flat2d(wmesh_int_t n_,wmesh_int_t l[])
{
  return l[0] + l[1]*n_;
}

#if 0


template<typename T>
static wmesh_status_t bms_transform_pyramid(wmesh_int_t 		r_size_,
					    const_wmesh_int_p 		r_ptr_,
					    const_wmesh_int_p 		r_m_,
					    const_wmesh_int_p 		r_n_,
					    const T*__restrict__ 	r_v_,
					    wmesh_int_t 		r_ld_,
					    
					    wmesh_int_t			c_storage_,
					    wmesh_int_t			c_m_,
					    wmesh_int_t			c_n_,
					    T*				c_v_,
					    wmesh_int_t			c_ld_,
					    wmesh_int_p 		perm_)

{ // 5
  // 4 1
  // 3  2
  // 2   3
  // 1    4
  // 0      5
  for (wmesh_int_t isize=0;isize<r_size_;++isize)
    {
      
    }
  
  wmesh_int_t ijk[3];
  for (wmesh_int_t k=0;k<r_n_;++k)
    {
      const auto rk = (r_v_[k*r_inc_] + one) / two;
      for (wmesh_int_t j=0;j<r_n_-k;j++)
	{
	  const auto rj = (r_v_[j*r_inc_] + one) / two;
	  for (wmesh_int_t i=0;i<r_n_-k;++i)
	    {		
	      const auto ri = (r_v_[i*r_inc_] + one) / two;
	      ijk[0] = i;
	      ijk[1] = j;
	      ijk[2] = k;
	      wmesh_int_t idx = perm_[flat3d(r_n_,ijk)]-1;	      
	      c_v_[c_ld_*idx+0] = ri;
	      c_v_[c_ld_*idx+1] = rj;
	      c_v_[c_ld_*idx+2] = rk;
	    } 
	}
    }
}
#endif 

template<typename T>
static wmesh_status_t bms_transform_tetrahedron(wmesh_int_t 		r_n_,
						const T*__restrict__ 	r_v_,
						wmesh_int_t 		r_inc_,
						wmesh_int_t		c_storage_,
						wmesh_int_t		c_m_,
						wmesh_int_t		c_n_,
						T*			c_v_,
						wmesh_int_t		c_ld_,
						wmesh_int_p 		perm_)
{
  static constexpr const T s_T0(0);
  static constexpr const T s_T1(1);
  static constexpr const T s_T2(2);
  static constexpr const T s_T3(3);

  //  static constexpr const T d_(2.0);
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
	  wmesh_int_t idx = perm_[flat3d(n_,ijk)]-1;
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
	  wmesh_int_t idx = perm_[flat3d(n_,ijk)]-1;

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
	  wmesh_int_t idx = perm_[flat3d(n_,ijk)]-1;


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
	  wmesh_int_t idx = perm_[flat3d(n_,ijk)]-1;
	  
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
	      
	      wmesh_int_t idx = perm_[flat3d(n_,ijk)]-1;
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
			     
			     wmesh_int_t		b_storage_,
			     wmesh_int_t		b_m_,
			     wmesh_int_t		b_n_,
			     const_wmesh_int_p		b_v_,
			     wmesh_int_t		b_ld_,
			     
			     wmesh_int_t		c_storage_,
			     wmesh_int_t		c_m_,
			     wmesh_int_t		c_n_,
			     T * __restrict__		c_v_,					 
			     wmesh_int_t		c_ld_,
			     
			     wmesh_int_t 		work_n_,
			     wmesh_int_p		work_) 
{

  static constexpr const T one(1.0);
  static constexpr const T two(2.0);
  static constexpr const T six(6.0);

  WMESH_CHECK_POSITIVE(r_n_);
  WMESH_CHECK_POINTER(r_v_);
  WMESH_CHECK(r_inc_ >= 1);
  
  WMESH_CHECK_STORAGE(b_storage_);
  
  wmesh_int_t b_dim 	= (b_storage_ == WMESH_STORAGE_INTERLEAVE) ? b_m_ : b_n_;
  wmesh_int_t c_dim 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_m_ : c_n_;
  wmesh_int_t b_count 	= (b_storage_ == WMESH_STORAGE_INTERLEAVE) ? b_n_ : b_m_;
  wmesh_int_t c_count 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_n_ : c_m_;
  
  WMESH_CHECK( (b_dim == 2) || (b_dim ==3) );  
  WMESH_CHECK( b_dim == c_dim );  
  WMESH_CHECK( b_count == c_count );  
  
  WMESH_CHECK_POINTER(b_ld_ >= b_m_);
  WMESH_CHECK_POINTER(c_ld_ >= c_m_);
  
  wmesh_int_t required_work_n = (b_dim == 2) ? r_n_ * r_n_ : r_n_ * r_n_ * r_n_;
  if (work_n_ < required_work_n)
    {
      std::cerr << "// WMESH::ERROR work_n_ = " << work_n_
		<< ", required_work_n_ = " << required_work_n
		<< ", r_n_ = " << r_n_
		<< std::endl;
      WMESH_CHECK(work_n_ >= required_work_n);
    }
  
  WMESH_CHECK_POINTER(work_);

  for (wmesh_int_t i=0;i<required_work_n;++i)
    {
      work_[i] = 0;
    }


  wmesh_int_t degree = r_n_ - 1;
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

	    wmesh_int_t flat  = (b_dim==2) ? flat2d(r_n_,ijk) : flat3d(r_n_,ijk);
	    work_[flat] = j + 1;
	    //	    std::cout << "perm["  << flat << "] = "  << j+1 <<  ",  "  << ijk[0] << " "  << ijk[1] << " "  << ijk[2] << std::endl;
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
	    wmesh_int_t flat  = (b_dim==2) ? flat2d(r_n_,ijk) : flat3d(r_n_,ijk);
	    work_[flat] = j + 1;
	  }
	break;
      }
    default:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
      }
    }
  
  
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
		    wmesh_int_t idx = work_[flat2d(r_n_,ijk)] - 1;
		    
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
		    WMESH_CHECK(work_[flat2d(r_n_,ijk)] > 0);
		    wmesh_int_t idx = work_[flat2d(r_n_,ijk)] - 1;
		    
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
	    if (2 == b_dim)
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
			WMESH_CHECK(work_[flat2d(r_n_,ijk)] > 0);
#endif
			wmesh_int_t idx = work_[flat2d(r_n_,ijk)] - 1;
			c_v_[c_ld_ * idx + 0] = ri;
			c_v_[c_ld_ * idx + 1] = rj;
		      }
		  }

	      }
	    else
	      {
#ifndef NDEBUG			    
		WMESH_CHECK(3 == b_dim);
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
			    WMESH_CHECK(work_[flat3d(r_n_,ijk)] > 0);
#endif
			    wmesh_int_t idx = work_[flat3d(r_n_,ijk)] - 1;
			    // std::cout << " >  idx = " << idx + 1 << ", c_ld_ = " << c_ld_ << ", flat = "<< flat3d(r_n_,ijk) << std::endl;
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
	    if (2 == b_dim)
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
			WMESH_CHECK(work_[flat2d(r_n_,ijk)] > 0);
#endif
			wmesh_int_t idx = work_[flat2d(r_n_,ijk)] - 1;
			c_v_[c_ld_ * 0 + idx] = ri;
			c_v_[c_ld_ * 1 + idx] = rj;
		      }
		  }

	      }
	    else
	      {
#ifndef NDEBUG			    		
		WMESH_CHECK(3 == b_dim);
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
			    WMESH_CHECK(work_[flat3d(r_n_,ijk)] > 0);
#endif
			    wmesh_int_t idx = work_[flat3d(r_n_,ijk)] - 1;			
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
	    WMESH_CHECK(3 == b_dim);
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
			WMESH_CHECK(work_[flat3d(r_n_,ijk)] > 0);
#endif
			wmesh_int_t idx 	= work_[flat3d(r_n_,ijk)] - 1;		    
			
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
			WMESH_CHECK(work_[flat3d(r_n_,ijk)] > 0);
#endif
			wmesh_int_t idx 	= work_[flat3d(r_n_,ijk)] - 1;		    

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
#if 0
      wmesh_status_t status =  bms_transform_pyramid(r_n_,
						     r_v_,
						     r_inc_,
						     c_storage_,
						     c_m_,
						     c_n_,
						     c_v_,
						     c_ld_,
						     work_);
      WMESH_STATUS_CHECK(status);
#endif
      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
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
							 work_);
      WMESH_STATUS_CHECK(status);      
    }
  return WMESH_STATUS_SUCCESS;
};


template
wmesh_status_t bms_transform<float>(wmesh_int_t 		element_,
				    wmesh_int_t 		r_n_,					 
				    const float*__restrict__	r_,
				    wmesh_int_t			r_inc_,
				    
				    wmesh_int_t			b_storage_,
				    wmesh_int_t			b_m_,
				    wmesh_int_t			b_n_,
				    const_wmesh_int_p		b_v_,
				    wmesh_int_t			b_ld_,
				    
				    wmesh_int_t			c_storage_,
				    wmesh_int_t			c_m_,
				    wmesh_int_t			c_n_,
				    float * __restrict__	c_v_,					 
				    wmesh_int_t			c_ld_,
				    
				    wmesh_int_t 		work_n_,
				    wmesh_int_p			work_) ;


template
wmesh_status_t bms_transform<double>(wmesh_int_t 		element_,
				    wmesh_int_t 		r_n_,					 
				    const double*__restrict__	r_,
				    wmesh_int_t			r_inc_,
				    
				    wmesh_int_t			b_storage_,
				    wmesh_int_t			b_m_,
				    wmesh_int_t			b_n_,
				    const_wmesh_int_p		b_v_,
				    wmesh_int_t			b_ld_,
				    
				    wmesh_int_t			c_storage_,
				    wmesh_int_t			c_m_,
				    wmesh_int_t			c_n_,
				    double * __restrict__	c_v_,					 
				    wmesh_int_t			c_ld_,
				    
				    wmesh_int_t 		work_n_,
				    wmesh_int_p			work_) ;


template
wmesh_status_t bms_transform<wmesh_int_t>(wmesh_int_t 		element_,
				    wmesh_int_t 		r_n_,					 
				    const wmesh_int_t*__restrict__	r_,
				    wmesh_int_t			r_inc_,
				    
				    wmesh_int_t			b_storage_,
				    wmesh_int_t			b_m_,
				    wmesh_int_t			b_n_,
				    const_wmesh_int_p		b_v_,
				    wmesh_int_t			b_ld_,
				    
				    wmesh_int_t			c_storage_,
				    wmesh_int_t			c_m_,
				    wmesh_int_t			c_n_,
				    wmesh_int_t * __restrict__	c_v_,					 
				    wmesh_int_t			c_ld_,
				    
				    wmesh_int_t 		work_n_,
				    wmesh_int_p			work_) ;


extern "C"
{
  wmesh_status_t bms_ordering_flat(wmesh_int_t		degree_,
				   wmesh_int_t		b_storage_,
				   wmesh_int_t		b_m_,
				   wmesh_int_t		b_n_,
				   const_wmesh_int_p	b_v_,
				   wmesh_int_t		b_ld_,
				   wmesh_int_p		flat_)
  {
    wmesh_int_t r_n_ = degree_ + 1;
    wmesh_int_t b_dim 	= (b_storage_ == WMESH_STORAGE_INTERLEAVE) ? b_m_ : b_n_;
    wmesh_int_t b_count 	= (b_storage_ == WMESH_STORAGE_INTERLEAVE) ? b_n_ : b_m_;  
    WMESH_CHECK( (b_dim == 2) || (b_dim ==3) );    
    WMESH_CHECK_POINTER(b_ld_ >= b_m_);
  
    for (wmesh_int_t i=0;i<b_count;++i)
      {
	flat_[i] = 0;
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

	      wmesh_int_t flat  = (b_dim==2) ? flat2d(r_n_,ijk) : flat3d(r_n_,ijk);
	      flat_[flat] = j + 1;
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
	      wmesh_int_t flat  = (b_dim==2) ? flat2d(r_n_,ijk) : flat3d(r_n_,ijk);
	      flat_[flat] = j + 1;
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
