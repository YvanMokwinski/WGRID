#pragma once
#include "bms_traits.hpp"
#include "wmesh-blas.hpp"
#include "wmesh-math.hpp"
#include <iostream>
template <typename T> constexpr T Factorial(wmesh_int_t i)
{
  return (i==0) ? T(1.0) : T(i)*Factorial<T>(i-1);
};

template <typename T> constexpr T Gamma(wmesh_int_t i)
{
  return Factorial<T>(i-1);
};

template <typename T> constexpr T Pow2(wmesh_int_t i)
{
  return (i==0) ? T(1.0) : T(2.0) * Pow2<T>(i-1);
};

#if 0
template <typename T>
wmesh_status_t bms_jacobip(wmesh_int_t 			alpha_,
			   wmesh_int_t 			beta_,
			   wmesh_int_t 			N_,
			   wmesh_int_t 			x_n_,
			   const T * __restrict__  	x_,
			   wmesh_int_t  		x_ld_,
			   T *  __restrict__ 		y_,
			   wmesh_int_t  		y_ld_,
			   wmesh_int_t 			work_n_,			   
			   T *  __restrict__ 		work_)
{
  if (N_ == 0)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  y_[y_ld_*i]=static_cast<T>(1.0);
	}
    }
  else if (N_ == 1)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  y_[y_ld_*i]=x_[x_ld_*i];
	}
    }
  else if (N_ == 2)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  y_[y_ld_*i]=(x_[x_ld_*i]*x_[x_ld_*i]*3.0-1.0)/2.0;
	}
    }
  return 0;
  
}

#endif
#if 0



template <typename T>
wmesh_status_t bms_jacobip(wmesh_int_t 			alpha_,
			   wmesh_int_t 			beta_,
			   wmesh_int_t 			N_,
			   wmesh_int_t 			x_n_,
			   const T * __restrict__  	x_,
			   wmesh_int_t  		x_ld_,
			   T *  __restrict__ 		y_,
			   wmesh_int_t  		y_ld_,
			   wmesh_int_t 			work_n_,			   
			   T *  __restrict__ 		work_)
{
  WMESH_CHECK(alpha_ >= 0);
  WMESH_CHECK(beta_  >= 0);
  WMESH_CHECK(N_     >= 0);
  WMESH_CHECK_POSITIVE(x_n_);
  WMESH_CHECK(work_n_  >= 2*x_n_);
  WMESH_CHECK_POINTER(x_);
  WMESH_CHECK_POINTER(y_);
  WMESH_CHECK_POINTER(work_);

  static constexpr T
    r2 = T(2.0),
    r3 = T(3.0),
    r1 = T(1.0),
    r0 = T(0.0);

  for (wmesh_int_t i=0;i<x_n_;++i)
    {
      y_[i * y_ld_] = r1;
    }
  
  wmesh_int_t n1=1;
  T * __restrict__ pii 	= work_;
  T * __restrict__  pi 	= work_ + x_n_;
  if (N_>0)
    {
      xcopy(&x_n_,y_,&y_ld_,pii,&n1);
      for (wmesh_int_t i=0;i<x_n_;++i)
	{      
	  y_[i*y_ld_] = (alpha_ + 1) + ((alpha_+beta_+2)*(x_[i *x_ld_]-1.0)) / 2.0; 
	}
      
      if (N_>1)
	{
	  //
	  // pi = y
	  //
	  xcopy(&x_n_,y_,&y_ld_,pi,&n1);
	  for (int i=2; i<= N_; ++i)
	    {
	      h1 = r2*i+ab;
	      T ri = static_cast<T>(i);
	      anew = r2/(h1+r2)*wmesh_math<T>::xsqrt( (ri+1.0)*(ri+ab1)*(ri+a1)*(ri+b1)/(h1+r1)/(h1+r3));
	      bnew = -(alpha_*alpha_-beta_*beta_) / ( h1*(h1+r2) );
	      
	      for (wmesh_int_t j=0;j<x_n_;++j)
		{
		  y_[j*y_ld_] = ( (x_[j*x_ld_]-bnew) * pi[i] - (2.0*i +  ) *pii[i] ) / (  2.0*i*(i+alpha_+beta_)*(2.0*i+alpha_+beta_-2)  );
		}
	      for (wmesh_int_t j=0;j<x_n_;++j)
		{
		  y_[j*y_ld_] *= r1/anew;
		}
	      
	      //	      y = (x-bnew) * pi - aold*pii;
	      //	      BLAS_daxpy(&x_n_,&s,y_,&n1);
	      //	      y *= r1/anew;
	      //	      T s = r1 / anew;
	      //	      BLAS_dscal(&x_n_,&s,y_,&n1);
	      aold = anew;
	      xcopy(&x_n_,pi,&n1,pii,&n1);
	      xcopy(&x_n_,y_,&y_ld_,pi,&n1);
	      //  pii = pi;
	      //  pi = y;
	    }	  
	}     
    }

  return WMESH_STATUS_SUCCESS;
};
#endif



#if 1


template <typename T>
wmesh_status_t bms_jacobip00(wmesh_int_t 			N_,
			     wmesh_int_t 			x_n_,
			     const T * __restrict__  	x_,
			     wmesh_int_t  		x_ld_,
			     T *  __restrict__ 		y_,
			     wmesh_int_t  		y_ld_,
			     wmesh_int_t 			work_n_,			   
			     T *  __restrict__ 		work_)
{
  if (N_==0)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  y_[i * y_ld_] = 1.0;
	}            
    }
  else if (N_==1)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  y_[i * y_ld_] = x_[i*x_ld_];
	}            
    }
  else if (N_==2)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = (t*t*3.0-1.0)*0.5;
	}            
    }
  else if (N_==3)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = t*(5.0*t*t-3.0)/2.0;
	}
    }
  else if (N_==4)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = (35.0*t*t*t*t-30.0*t*t+3.0)/8.0;
	}
    }
  else
    {
      std::cerr << "error" << std::endl;
      return 1;
    }
  return 0;
}


template <typename T>
wmesh_status_t bms_jacobip11(wmesh_int_t 			N_,
			     wmesh_int_t 			x_n_,
			     const T * __restrict__  	x_,
			     wmesh_int_t  		x_ld_,
			     T *  __restrict__ 		y_,
			     wmesh_int_t  		y_ld_,
			     wmesh_int_t 			work_n_,			   
			     T *  __restrict__ 		work_)
{
  if (N_==0)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  y_[i * y_ld_] = 0.0;
	}            
    }
  else if (N_==1)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  y_[i * y_ld_] = 1.0;
	}            
    }
  else if (N_==2)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = 3.0*x_[i*x_ld_];
	}            
    }
  else if (N_==3)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = (5.0*t*t-1.0)*3.0/2.0;
	}
    }
  else if (N_==4)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = 5.0*t*(7.0*t*t-3.0)/2.0;
	}
    }
  else
    {
      std::cerr << "error" << std::endl;
      return 1;
    }
  return 0;

}
#else

template <typename T>
wmesh_status_t bms_jacobip00(wmesh_int_t 			N_,
			     wmesh_int_t 			x_n_,
			     const T * __restrict__  	x_,
			     wmesh_int_t  		x_ld_,
			     T *  __restrict__ 		y_,
			     wmesh_int_t  		y_ld_,
			     wmesh_int_t 			work_n_,			   
			     T *  __restrict__ 		work_)
{
  if (N_==0)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  y_[i * y_ld_] = 1.0;
	}            
    }
  else if (N_==1)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  y_[i * y_ld_] = x_[i*x_ld_];
	}            
    }
  else if (N_==2)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = t*t;
	}            
    }
  else if (N_==3)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = t*t*t;
	}
    }
  else if (N_==4)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = t*t*t*t;
	}
    }
  else
    {
      std::cerr << "error" << std::endl;
      return 1;
    }
  return 0;
}


template <typename T>
wmesh_status_t bms_jacobip11(wmesh_int_t 			N_,
			     wmesh_int_t 			x_n_,
			     const T * __restrict__  	x_,
			     wmesh_int_t  		x_ld_,
			     T *  __restrict__ 		y_,
			     wmesh_int_t  		y_ld_,
			     wmesh_int_t 			work_n_,			   
			     T *  __restrict__ 		work_)
{
  if (N_==0)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  y_[i * y_ld_] = 0.0;
	}            
    }
  else if (N_==1)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  y_[i * y_ld_] = 1.0;
	}            
    }
  else if (N_==2)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = 2.0*t;
	}            
    }
  else if (N_==3)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = 3.0*t*t;
	}
    }
  else if (N_==4)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = 4.0*t*t*t;
	}
    }
  else
    {
      std::cerr << "error" << std::endl;
      return 1;
    }
  return 0;

}
#endif


template<typename T>
wmesh_status_t bms_jacobip11old(wmesh_int_t 			N_,
			     wmesh_int_t 			x_n_,
			     const T * __restrict__  	x_,
			     wmesh_int_t  		x_ld_,
			     T *  __restrict__ 		y_,
			     wmesh_int_t  		y_ld_,
			     wmesh_int_t 			work_n_,			   
			     T *  __restrict__ 		work_)
{
  if (N_==0)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  y_[i * y_ld_] = 1.0;
	}            
    }
  else if (N_==1)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  y_[i * y_ld_] = 2.0*x_[i*x_ld_];
	}            
    }
  else if (N_==2)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = 3.0*(5.0*t*t-1.0)*0.25;
	}            
    }
  else if (N_==3)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = t*(7.0*t*t-3.0);
	}
    }
  else if (N_==4)
    {
      for (wmesh_int_t i=0;i<x_n_;++i)
	{
	  T t = x_[i*x_ld_];
	  y_[i * y_ld_] = 5.0*(21.0*t*t*t*t-14.0*t*t+1.0)/8.0;
	}
    }
  else
    {
      std::cerr << "error" << std::endl;
      return 1;
    }
  return 0;

}

template <typename T>
wmesh_status_t bms_jacobip(wmesh_int_t 			alpha_,
			   wmesh_int_t 			beta_,
			   wmesh_int_t 			N_,
			   wmesh_int_t 			x_n_,
			   const T * __restrict__  	x_,
			   wmesh_int_t  		x_ld_,
			   T *  __restrict__ 		y_,
			   wmesh_int_t  		y_ld_,
			   wmesh_int_t 			work_n_,			   
			   T *  __restrict__ 		work_)
{
  WMESH_CHECK(alpha_ >= 0);
  WMESH_CHECK(beta_  >= 0);
  WMESH_CHECK(N_     >= 0);
  WMESH_CHECK_POSITIVE(x_n_);
  WMESH_CHECK(work_n_  >= 2*x_n_);
  WMESH_CHECK_POINTER(x_);
  WMESH_CHECK_POINTER(y_);
  WMESH_CHECK_POINTER(work_);
#if 0
  if ( (alpha_==0)&&(beta_==0))
    {
      return bms_jacobip00(N_,
			 x_n_,
			 x_,
			 x_ld_,
			 y_,
			 y_ld_,
			 work_n_,			   
			 work_);
    }
  if ( (alpha_==1)&&(beta_==1))
    {
      return bms_jacobip11(N_,
			 x_n_,
			 x_,
			 x_ld_,
			 y_,
			 y_ld_,
			 work_n_,			   
			 work_);
    }

  std::cerr << "pas ici " << std::endl;
  exit(1);
#endif

  
  static constexpr T
    r2 = T(2.0),
    r3 = T(3.0),
    r1 = T(1.0),
    r0 = T(0.0);
  
  T aold   = r0,
    anew   = r0,
    bnew   = r0,
    h1     = r0,
    gamma1 = r0;
  
  const T
    ab  = alpha_ + beta_,
    ab1 = alpha_+ beta_ + 1,
    a1  = alpha_ + 1,
    b1  = beta_ + 1;
#if 0
  std::cout <<"alpha " << alpha_ << std::endl;
  std::cout <<"beta " << beta_ << std::endl;
  std::cout <<"N_ " << N_ << std::endl;
  std::cout <<"ab " << ab << std::endl;
  std::cout <<"a1 " << a1 << std::endl;
  std::cout <<"b1 " << b1 << std::endl;
#endif
  
  // Initial values P_0(x) and P_1(x)
  const T gamma0 = Pow2<T>(ab1)*Gamma<T>(a1)*Gamma<T>(b1)/Factorial<T>(ab1);
  const T y0 = r1 / wmesh_math<T>::xsqrt(gamma0);


  for (wmesh_int_t i=0;i<x_n_;++i)
    {
      y_[i * y_ld_] = 1.0;
    }      
  wmesh_int_t n1=1;
  T * __restrict__ pii 	= work_;
  T * __restrict__  pi 	= work_ + x_n_;
  if (N_>0)
    {
      xcopy(&x_n_,y_,&y_ld_,pii,&n1);
      for (wmesh_int_t i=0;i<x_n_;++i)
	{      
	  y_[i*y_ld_] = (1.0+alpha_) + (2.0+alpha_+beta_)*(x_[i *x_ld_]-1.0)/2.0;
	}
      if (N_>1)
	{	  
	  
	  //
	  // pi = y
	  //
	  xcopy(&x_n_,y_,&y_ld_,pi,&n1);
#if 1
	  for (wmesh_int_t i=2; i<=N_; ++i)
	    {
	      T ri = i;
	      T c0 = 2.0*ri*(ri+alpha_+beta_)*(2.0*ri+alpha_+beta_-2.0);
	      T c1 = -2.0*(-1.0+ri+alpha_)*(-1.0+ri+beta_)*(2.0*ri+beta_+alpha_);
	      for (wmesh_int_t j=0;j<x_n_;++j)
		{
		  T c2 = (2.0*ri+alpha_+beta_-1.0)  * (   (2.0*ri+alpha_+beta_-2.0)*(2.0*ri+alpha_+beta_)*x_[j*x_ld_] +alpha_*alpha_-beta_*beta_ );
		  y_[j*y_ld_] = ( c1 * pii[j] + c2 * pi[j]) / c0;
		  // fprintf(stdout,"%e %e c0=%e, alpha_ %ld beta %ld\n",y_[j*y_ld_],(3.0*x_[x_ld_*j]*x_[x_ld_*j]-1.0)*0.5,c0,alpha_,beta_);
		  //		  y_[j*y_ld_] = ( ( (alpha_*alpha_-beta_*beta_)* + x_[i*x_ld_]*(2.0*i+alpha_+beta_-2.0) )*pi[j]*pii[j])/;
		}
	      //   exit(1);
	      xcopy(&x_n_,pi,&n1,pii,&n1);
	      xcopy(&x_n_,y_,&y_ld_,pi,&n1);
	    }
#endif
	  //	  exit(1);
	}     
    }

  return WMESH_STATUS_SUCCESS;
};




template <typename T>
wmesh_status_t bms_jacobipold(wmesh_int_t 			alpha_,
			   wmesh_int_t 			beta_,
			   wmesh_int_t 			N_,
			   wmesh_int_t 			x_n_,
			   const T * __restrict__  	x_,
			   wmesh_int_t  		x_ld_,
			   T *  __restrict__ 		y_,
			   wmesh_int_t  		y_ld_,
			   wmesh_int_t 			work_n_,			   
			   T *  __restrict__ 		work_)
{
  WMESH_CHECK(alpha_ >= 0);
  WMESH_CHECK(beta_  >= 0);
  WMESH_CHECK(N_     >= 0);
  WMESH_CHECK_POSITIVE(x_n_);
  WMESH_CHECK(work_n_  >= 2*x_n_);
  WMESH_CHECK_POINTER(x_);
  WMESH_CHECK_POINTER(y_);
  WMESH_CHECK_POINTER(work_);

  static constexpr T
    r2 = T(2.0),
    r3 = T(3.0),
    r1 = T(1.0),
    r0 = T(0.0);
  
  T aold   = r0,
    anew   = r0,
    bnew   = r0,
    h1     = r0,
    gamma1 = r0;
  
  const T
    ab  = alpha_ + beta_,
    ab1 = alpha_+ beta_ + 1,
    a1  = alpha_ + 1,
    b1  = beta_ + 1;
#if 0
  std::cout <<"alpha " << alpha_ << std::endl;
  std::cout <<"beta " << beta_ << std::endl;
  std::cout <<"N_ " << N_ << std::endl;
  std::cout <<"ab " << ab << std::endl;
  std::cout <<"a1 " << a1 << std::endl;
  std::cout <<"b1 " << b1 << std::endl;
#endif
  
  // Initial values P_0(x) and P_1(x)
  const T gamma0 = Pow2<T>(ab1)*Gamma<T>(a1)*Gamma<T>(b1)/Factorial<T>(ab1);
  const T y0 = r1 / wmesh_math<T>::xsqrt(gamma0);
#if 0
  std::cout <<"gamma(a1) " << Gamma<T>(a1) << std::endl;
  std::cout <<"gamma(a1) " << Gamma<T>(a1) << std::endl;
  std::cout <<"factorial(b1) " << Factorial<T>(ab1) << std::endl;
  std::cout <<"r1 " << r1 << std::endl;
  std::cout <<"gammmaaaa " << gamma0 << std::endl;
  std::cout <<"y0 " << y0 << std::endl;
#endif  
  for (wmesh_int_t i=0;i<x_n_;++i)
    {
      y_[i * y_ld_] = y0;
    }
  
  wmesh_int_t n1=1;
  T * __restrict__ pii 	= work_;
  T * __restrict__  pi 	= work_ + x_n_;
  if (N_>0)
    {
      xcopy(&x_n_,y_,&y_ld_,pii,&n1);
      //      auto pii = y;
      gamma1 = (a1)*(b1)/(ab+3.0)*gamma0;
      for (wmesh_int_t i=0;i<x_n_;++i)
	{      
	  y_[i*y_ld_] = ((ab+r2)*x_[i *x_ld_]/r2 + (alpha_-beta_)/r2) / wmesh_math<T>::xsqrt(gamma1);
	}
      if (N_>1)
	{

	  //
	  // pi = y
	  //
	  xcopy(&x_n_,y_,&y_ld_,pi,&n1);

	  // Repeat value in recurrence.
	  aold = r2 / (r2+ab) * wmesh_math<T>::xsqrt((a1)*(b1)/(ab+r3));	  
	  // Forward recurrence using the symmetry of the recurrence.
	  for (int i=1; i<=(N_-1); ++i)
	    {
	      h1 = r2*i+ab;
	      T ri = static_cast<T>(i);
	      anew = r2/(h1+r2)*wmesh_math<T>::xsqrt( (ri+1.0)*(ri+ab1)*(ri+a1)*(ri+b1)/(h1+r1)/(h1+r3));
	      bnew = -(alpha_*alpha_-beta_*beta_) / ( h1*(h1+r2) );
	      
	      for (wmesh_int_t j=0;j<x_n_;++j)
		{
		  y_[j*y_ld_] = (x_[j*x_ld_]-bnew) * pi[j] - aold*pii[j];
		}
	      for (wmesh_int_t j=0;j<x_n_;++j)
		{
		  y_[j*y_ld_] *= r1/anew;
		}
	      
	      //	      y = (x-bnew) * pi - aold*pii;
	      //	      BLAS_daxpy(&x_n_,&s,y_,&n1);
	      //	      y *= r1/anew;
	      //	      T s = r1 / anew;
	      //	      BLAS_dscal(&x_n_,&s,y_,&n1);
	      aold = anew;
	      xcopy(&x_n_,pi,&n1,pii,&n1);
	      xcopy(&x_n_,y_,&y_ld_,pi,&n1);
	      //  pii = pi;
	      //  pi = y;
	    }	  
	}     
    }

  return WMESH_STATUS_SUCCESS;
};


template<wmesh_int_t TOPODIM_>
inline wmesh_status_t bms_template_topodim2numtypes(wmesh_int_p 	ntypes_)
{
  WMESH_CHECK_POINTER(ntypes_);
  ntypes_[0] = bms_traits_topodim<TOPODIM_>::s_ntypes;
  return WMESH_STATUS_SUCCESS;
}

template<wmesh_int_t TOPODIM_>
wmesh_status_t bms_template_topodim2elements(wmesh_int_p 	num_elements_,
					     wmesh_int_p 	elements_)
{
  num_elements_[0] = bms_traits_topodim<TOPODIM_>::s_ntypes;
  for (wmesh_int_t i=0;i<bms_traits_topodim<TOPODIM_>::s_ntypes;++i)
    {
      elements_[i] = bms_traits_topodim<TOPODIM_>::s_ntypes + i;
    }
  return WMESH_STATUS_SUCCESS;
}

template<wmesh_int_t ELEMENT_>
inline wmesh_status_t bms_template_element2topodim(wmesh_int_p 	topodim_)
{
  WMESH_CHECK_POINTER(topodim_);
  topodim_[0] = bms_traits_element<ELEMENT_>::s_topodim;
  return WMESH_STATUS_SUCCESS;
}

template<wmesh_int_t TOPODIM_>
inline wmesh_status_t bms_template_elements_num_facets(wmesh_int_p 	num_facets_);

template<>
inline wmesh_status_t bms_template_elements_num_facets<WMESH_TOPODIM_VOLUME>(wmesh_int_p 	num_facets_)
{
  WMESH_CHECK_POINTER(num_facets_);
  num_facets_[0] = 4;
  num_facets_[1] = 5;
  num_facets_[2] = 5;
  num_facets_[3] = 6;
  return WMESH_STATUS_SUCCESS;
}

template<>
inline wmesh_status_t bms_template_elements_num_facets<WMESH_TOPODIM_FACE>(wmesh_int_p 	num_facets_)
{
  WMESH_CHECK_POINTER(num_facets_);
  num_facets_[0] = 3;
  num_facets_[1] = 4;
  return WMESH_STATUS_SUCCESS;
}


template<>
inline wmesh_status_t bms_template_elements_num_facets<WMESH_TOPODIM_EDGE>(wmesh_int_p 	num_facets_)
{
  WMESH_CHECK_POINTER(num_facets_);
  num_facets_[0] = 2;
  return WMESH_STATUS_SUCCESS;
}

template<>
inline wmesh_status_t bms_template_elements_num_facets<WMESH_TOPODIM_NODE>(wmesh_int_p 	num_facets_)
{
  WMESH_CHECK_POINTER(num_facets_);
  num_facets_[0] = 0;
  return WMESH_STATUS_SUCCESS;
}


template<wmesh_int_t ELEMENT_>
inline wmesh_status_t bms_template_element_facets(wmesh_int_p num_facets_,wmesh_int_p facets_);

template<>
inline wmesh_status_t bms_template_element_facets<WMESH_ELEMENT_EDGE>(wmesh_int_p num_facets_,wmesh_int_p facets_)
{
  WMESH_CHECK_POINTER(facets_);
  num_facets_[0] = 2;
  facets_[0] = WMESH_ELEMENT_NODE;
  facets_[1] = WMESH_ELEMENT_NODE;
  return WMESH_STATUS_SUCCESS;
}


template<>
inline wmesh_status_t bms_template_element_facets<WMESH_ELEMENT_TRIANGLE>(wmesh_int_p num_facets_,wmesh_int_p facets_)
{
  WMESH_CHECK_POINTER(facets_);
  num_facets_[0] = 3;
  facets_[0] = WMESH_ELEMENT_EDGE;
  facets_[1] = WMESH_ELEMENT_EDGE;
  facets_[2] = WMESH_ELEMENT_EDGE;
  return WMESH_STATUS_SUCCESS;
}

template<>
inline wmesh_status_t bms_template_element_facets<WMESH_ELEMENT_QUADRILATERAL>(wmesh_int_p num_facets_,wmesh_int_p facets_)
{
  WMESH_CHECK_POINTER(facets_);
  num_facets_[0] = 4;
  facets_[0] = WMESH_ELEMENT_EDGE;
  facets_[1] = WMESH_ELEMENT_EDGE;
  facets_[2] = WMESH_ELEMENT_EDGE;
  facets_[3] = WMESH_ELEMENT_EDGE;
  return WMESH_STATUS_SUCCESS;
}



template<>
inline wmesh_status_t bms_template_element_facets<WMESH_ELEMENT_TETRAHEDRON>(wmesh_int_p num_facets_,wmesh_int_p facets_)
{
  WMESH_CHECK_POINTER(facets_);
  num_facets_[0] = 4;
  facets_[0] = WMESH_ELEMENT_TRIANGLE;
  facets_[1] = WMESH_ELEMENT_TRIANGLE;
  facets_[2] = WMESH_ELEMENT_TRIANGLE;
  facets_[3] = WMESH_ELEMENT_TRIANGLE;
  return WMESH_STATUS_SUCCESS;
}

template<>
inline wmesh_status_t bms_template_element_facets<WMESH_ELEMENT_PYRAMID>(wmesh_int_p num_facets_,wmesh_int_p facets_)
{
  WMESH_CHECK_POINTER(facets_);
  num_facets_[0] = 5;

  facets_[0] = WMESH_ELEMENT_TRIANGLE;
  facets_[1] = WMESH_ELEMENT_TRIANGLE;
  facets_[2] = WMESH_ELEMENT_TRIANGLE;
  facets_[3] = WMESH_ELEMENT_TRIANGLE;
  facets_[4] = WMESH_ELEMENT_QUADRILATERAL;
  return WMESH_STATUS_SUCCESS;
}

template<>
inline wmesh_status_t bms_template_element_facets<WMESH_ELEMENT_WEDGE>(wmesh_int_p num_facets_,wmesh_int_p facets_)
{
  WMESH_CHECK_POINTER(facets_);
  num_facets_[0] = 5;

  facets_[0] = WMESH_ELEMENT_TRIANGLE;
  facets_[1] = WMESH_ELEMENT_TRIANGLE;
  facets_[2] = WMESH_ELEMENT_QUADRILATERAL;
  facets_[3] = WMESH_ELEMENT_QUADRILATERAL;
  facets_[4] = WMESH_ELEMENT_QUADRILATERAL;
  return WMESH_STATUS_SUCCESS;
}

template<>
inline wmesh_status_t bms_template_element_facets<WMESH_ELEMENT_HEXAHEDRON>(wmesh_int_p num_facets_,wmesh_int_p facets_)
{
  WMESH_CHECK_POINTER(facets_);
  num_facets_[0] = 6;
  facets_[0] = WMESH_ELEMENT_QUADRILATERAL;
  facets_[1] = WMESH_ELEMENT_QUADRILATERAL;
  facets_[2] = WMESH_ELEMENT_QUADRILATERAL;
  facets_[3] = WMESH_ELEMENT_QUADRILATERAL;
  facets_[4] = WMESH_ELEMENT_QUADRILATERAL;
  facets_[5] = WMESH_ELEMENT_QUADRILATERAL;
  return WMESH_STATUS_SUCCESS;
}






template<wmesh_int_t ENTITY_>
inline wmesh_status_t bms_template_elements_num_entities(wmesh_int_t 		num_elements_,
							 const_wmesh_int_p 	elements_,
							 wmesh_int_p 		num_entities_)
{  
  static wmesh_int_t s_num_entities[8][8] = { {1,2,3,4,4,5,6,8},
					      {0,1,3,4,6,8,9,12},
					      {0,0,1,0,4,4,2,0},
					      {0,0,0,1,0,1,4,6},
					      {0,0,0,0,1,0,0,0},
					      {0,0,0,0,0,1,0,0},
					      {0,0,0,0,0,0,1,0},
					      {0,0,0,0,0,0,0,0} };
  for (wmesh_int_t i=0;i<num_elements_;++i)
    {
      num_entities_[i] = s_num_entities[ENTITY_][elements_[i]];
    }
  return WMESH_STATUS_SUCCESS;
}

template <wmesh_int_t ELEMENT> inline wmesh_int_t bms_template_ndofs(wmesh_int_t d_){ return -1; }

template <> inline wmesh_int_t bms_template_ndofs<WMESH_ELEMENT_EDGE>(wmesh_int_t d_) { return d_+1; }
template <> inline wmesh_int_t bms_template_ndofs<WMESH_ELEMENT_TRIANGLE>(wmesh_int_t d_) { return ((d_+1)*(d_+2)) / 2; }
template <> inline wmesh_int_t bms_template_ndofs<WMESH_ELEMENT_QUADRILATERAL>(wmesh_int_t d_) { return ((d_+1)*(d_+1)); }
template <> inline wmesh_int_t bms_template_ndofs<WMESH_ELEMENT_TETRAHEDRON>(wmesh_int_t d_) { return ((d_+1)*(d_+2)*(d_+3)) / 6; }
template <> inline wmesh_int_t bms_template_ndofs<WMESH_ELEMENT_PYRAMID>(wmesh_int_t d_) { return ((d_+2)*(d_+1)*(2*d_+3)) / 6; }
template <> inline wmesh_int_t bms_template_ndofs<WMESH_ELEMENT_WEDGE>(wmesh_int_t d_) { return ((d_+1)*(d_+1)*(d_+2)) / 2; }
template <> inline wmesh_int_t bms_template_ndofs<WMESH_ELEMENT_HEXAHEDRON>(wmesh_int_t d_) { return (d_+1)*(d_+1)*(d_+1); }


template <wmesh_int_t ELEMENT> inline wmesh_int_t bms_template_ndofs_interior(wmesh_int_t d_){ return -1;}

template <> inline wmesh_int_t bms_template_ndofs_interior<WMESH_ELEMENT_EDGE>(wmesh_int_t d_) { return (d_>0) ?d_-1 : 1; }
template <> inline wmesh_int_t bms_template_ndofs_interior<WMESH_ELEMENT_TRIANGLE>(wmesh_int_t d_) { return (d_>0) ?((d_-1)*(d_-2)) / 2 : 1; }
template <> inline wmesh_int_t bms_template_ndofs_interior<WMESH_ELEMENT_QUADRILATERAL>(wmesh_int_t d_) { return (d_>0) ?((d_-1)*(d_-1)) : 1; }
template <> inline wmesh_int_t bms_template_ndofs_interior<WMESH_ELEMENT_TETRAHEDRON>(wmesh_int_t d_) { return (d_>0) ?((d_-1)*(d_-2)*(d_-3)) / 6 : 1; }
template <> inline wmesh_int_t bms_template_ndofs_interior<WMESH_ELEMENT_PYRAMID>(wmesh_int_t d_) { return (d_>0) ?((d_-2)*(d_-1)*(2*d_-3)) / 6 : 1; }
template <> inline wmesh_int_t bms_template_ndofs_interior<WMESH_ELEMENT_WEDGE>(wmesh_int_t d_) { return (d_>0) ?((d_-1)*(d_-1)*(d_-2)) / 2 : 1; }
template <> inline wmesh_int_t bms_template_ndofs_interior<WMESH_ELEMENT_HEXAHEDRON>(wmesh_int_t d_) { return (d_>0) ? (d_-1)*(d_-1)*(d_-1) : 1; }






template<typename T>
wmesh_status_t bms_template_shape(wmesh_int_t 		element_,
				  wmesh_int_t 		family_,
				  wmesh_int_t 		degree_,

				  const_wmesh_int_p	diff_,
				  
				  wmesh_int_t 		c_storage_,					  
				  wmesh_int_t 		c_m_,
				  wmesh_int_t 		c_n_,
				  const T * 		c_,
				  wmesh_int_t 		c_ld_,
				  
				  wmesh_int_t 		b_storage_,
				  wmesh_int_t 		b_m_,
				  wmesh_int_t 		b_n_,
				  T* 			b_,
				  wmesh_int_t 		b_ld_,
				      
				  wmesh_int_t		iw_n_,
				  wmesh_int_p		iw_,
				  wmesh_int_t		rw_n_,
				  T * rw_);



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

	  wmesh_int_t		iwork_n_,
	  wmesh_int_p 		iwork_,
	  wmesh_int_t		rwork_n_,
	  T* __restrict__ 	rwork_);





template<typename T>
wmesh_status_t bms_mirrored_local_coordinates_edge(wmesh_int_t 			signed_rotation_,
						   wmesh_int_t 			c_storage_,
						   wmesh_int_t 			c_m_,
						   wmesh_int_t 			c_n_,
						   const T * __restrict__ 	c_v_,
						   wmesh_int_t 			c_ld_,
						   wmesh_int_t 			x_storage_,
						   wmesh_int_t 			x_m_,
						   wmesh_int_t 			x_n_,
						   T * __restrict__ 		x_v_,
						   wmesh_int_t 			x_ld_);

template<typename T>
wmesh_status_t bms_mirrored_local_coordinates_triangle(wmesh_int_t 			signed_rotation_,
						       wmesh_int_t 			c_storage_,
						       wmesh_int_t 			c_m_,
						       wmesh_int_t 			c_n_,
						       const T * __restrict__ 	c_v_,
						       wmesh_int_t 			c_ld_,
						       wmesh_int_t 			x_storage_,
						       wmesh_int_t 			x_m_,
						       wmesh_int_t 			x_n_,
						       T * __restrict__ 		x_v_,
						       wmesh_int_t 			x_ld_);;

template<typename T>
wmesh_status_t bms_mirrored_local_coordinates_quadrilateral(wmesh_int_t 			signed_rotation_,
							    wmesh_int_t 			c_storage_,
							    wmesh_int_t 			c_m_,
							    wmesh_int_t 			c_n_,
							    const T * __restrict__ 	c_v_,
							    wmesh_int_t 			c_ld_,
							    wmesh_int_t 			x_storage_,
							    wmesh_int_t 			x_m_,
							    wmesh_int_t 			x_n_,
							    T * __restrict__ 		x_v_,
							    wmesh_int_t 			x_ld_);

template<typename T>
wmesh_status_t bms_mirrored_local_coordinates(wmesh_int_t 	       		element_,
					      wmesh_int_t 			signed_rotation_,
					      wmesh_int_t 			c_storage_,
					      wmesh_int_t 			c_m_,
					      wmesh_int_t 			c_n_,
					      const T * __restrict__ 		c_v_,
					      wmesh_int_t 			c_ld_,
					      wmesh_int_t 			x_storage_,
					      wmesh_int_t 			x_m_,
					      wmesh_int_t 			x_n_,
					      T * __restrict__ 			x_v_,
					      wmesh_int_t 			x_ld_);

template<typename T>
wmesh_status_t
bms_template_cubature(wmesh_int_t	element_,
		      wmesh_int_t	family_,
		      wmesh_int_t	n1d_,

		      wmesh_int_t	c_storage_,
		      wmesh_int_t	c_m_,
		      wmesh_int_t	c_n_,
		      T*__restrict__ 	c_v_,
		      wmesh_int_t	c_ld_,
	     
		      wmesh_int_t	w_n_,
		      T*__restrict__ 	w_v_,
		      wmesh_int_t	w_inc_,
	     
		      wmesh_int_t	rwork_n_,
		      T* __restrict__ 	rwork_);

