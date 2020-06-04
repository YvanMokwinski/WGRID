#pragma once

template <wmesh_int_t ELEMENT,typename T>
struct bms_template_shape_vandermonde
{
  template <typename T>
  wmesh_status_t buffer_size(wmesh_int_t 		ndofs_,
			     wmesh_int_t 		num_points_,
			     wmesh_int_p 		iw_n_,
			     wmesh_int_p 		rw_n_);
  
  template <typename T>
  wmesh_status_t eval(wmesh_int_t 		degree_,
		      wmesh_int_t 		ndofs_,
		       
		      wmesh_int_t 		lcoo_storage_,
		      wmesh_int_t 		lcoo_m_,
		      wmesh_int_t 		lcoo_n_,
		      const T * __restrict__	lcoo_,
		      wmesh_int_t 		lcoo_ld_,
		       
		      wmesh_int_t 		v_storage_,
		      T*__restrict__		v_,
		      wmesh_int_t 		v_ld_,
		      
		      wmesh_int_t 		iw_n_,
		      wmesh_int_p		iw_,
		      wmesh_int_t 		rw_n_,
		      T*__restrict__		rw_);
};

template <typename T>
struct bms_template_shape_vandermonde<WMESH_ELEMENT_TRIANGLE,T>
{
  wmesh_status_t eval(wmesh_int_t 		degree_,
		      wmesh_int_t 		ndofs_,
		       
		      wmesh_int_t 		lcoo_storage_,
		      wmesh_int_t 		lcoo_m_,
		      wmesh_int_t 		lcoo_n_,
		      const T * __restrict__	lcoo_,
		      wmesh_int_t 		lcoo_ld_,
		       
		      wmesh_int_t 		v_storage_,
		      T*__restrict__		v_,
		      wmesh_int_t 		v_ld_,
		      
		      wmesh_int_t 		iw_n_,
		      wmesh_int_p		iw_,
		      wmesh_int_t 		rw_n_,
		      T*__restrict__		rw_)
  {
    const T * lr = lcoo_;
    const T * ls = lcoo_ + ( (lcoo_storage_ == WMESH_STORAGE_INTERLEAVE) ? 1 : lcoo_ld_ );
    wmesh_int_t lr_inc = ( (lcoo_storage_ == WMESH_STORAGE_INTERLEAVE) ? lcoo_ld_ : 1 );
    wmesh_int_t ls_inc = ( (lcoo_storage_ == WMESH_STORAGE_INTERLEAVE) ? lcoo_ld_ : 1 );
    wmesh_int_t num_points = ( (lcoo_storage_ == WMESH_STORAGE_INTERLEAVE) ? lcoo_n_ : lcoo_m );
    for (wmesh_int_t i=0;i<num_points;++i)
      {
	v_[i*v_ld_+0] = 1.0;
      }
    if (degree_ > 0)
      {
	
	for (wmesh_int_t i=0;i<num_points;++i)
	  {
	    v_[i*v_ld_+1] = lr[i*lr_inc];
	  }

	for (wmesh_int_t i=0;i<num_points;++i)
	  {
	    v_[i*v_ld_+2] = ls[i*ls_inc];
	  }
	
	if (degree_ > 1)
	  {
	    for (wmesh_int_t d = 2;d<=degree_;++d)
	    {
	      wmesh_int_t pshift = ((degree_-1)*(degree_))/2;
	      wmesh_int_t shift  = ((degree_)*(degree_+1))/2;
	      for (wmesh_int_t k=0;k<d ;++k)
		{
		  for (wmesh_int_t i=0;i<num_points;++i)
		    {
		      v_[i*v_ld_+ (shift+k)] = v_[i*v_ld_+ (pshift+k)] * lr[i*lr_inc];
		    }
		}
	      for (wmesh_int_t k=d;k<=d ;++k)
		{
		  for (wmesh_int_t i=0;i<num_points;++i)
		    {
		      v_[i*v_ld_+ (shift+k)] = v_[i*v_ld_+ (shift-1)] * ls[i*lr_inc];
		    }
		}
	    }	    
	  }	
      }
    return WMESH_STATUS_SUCCESS;
  }
};

template <typename T>
struct bms_template_shape_vandermonde<WMESH_ELEMENT_TRIANGLE,T>
{
  wmesh_status_t eval(wmesh_int_t 		degree_,
		      wmesh_int_t 		ndofs_,
		       
		      wmesh_int_t 		dofcoo_storage_,
		      wmesh_int_t 		dofcoo_m_,
		      wmesh_int_t 		dofcoo_n_,
		      const T * __restrict__	dofcoo_,
		      wmesh_int_t 		dofcoo_ld_,
		       
		      wmesh_int_t 		lcoo_storage_,
		      wmesh_int_t 		lcoo_m_,
		      wmesh_int_t 		lcoo_n_,
		      const T * __restrict__	lcoo_,
		      wmesh_int_t 		lcoo_ld_,
		       
		      wmesh_int_t 		v_storage_,
		      T*__restrict__		v_,
		      wmesh_int_t 		v_ld_,
		      
		      wmesh_int_t 		iw_n_,
		      wmesh_int_p		iw_,
		      wmesh_int_t 		rw_n_,
		      T*__restrict__		rw_)
  {
    WMESH_CHECK(work_n_ >= 5 * lcoo_n_);
    static constexpr double s1 = 1.0;
    static constexpr double s2 = 2.0;
    
    double * x 	= work_ + 0*lcoo_n_;
    double * yi 	= work_ + 1*lcoo_n_;
    double * yj 	= work_ + 2*lcoo_n_;
    
    work_ 	= work_ + 3*lcoo_n_;
    work_n_	-=  3*lcoo_n_;

    wmesh_int_t idx = 0;
    for (wmesh_int_t i=0;i<=d_;++i)
      {
	for (wmesh_int_t k = 0; k < lcoo_n_;++k)
	  {
	    double r 	= lcoo_[lcoo_ld_*k+0];
	    double s 	= lcoo_[lcoo_ld_*k+1];
	    x[k] 		= (s < s1) ? (r * s2) / (s1 - s) - s1 : 0.0;
	  }

	bms_djacobip(0,
		     0,
		     i,		
		     lcoo_n_,
		     x,
		     1,
		     yi,
		     1,
		     work_n_,
		     work_);
      
	for (wmesh_int_t k = 0; k < lcoo_n_;++k)
	  {
	    double s 	= lcoo_[lcoo_ld_*k+1];
	    yi[k] 	*= wmesh_math<double>::xpow(s1 - s, i);
	  }
      
	for (wmesh_int_t k = 0; k < lcoo_n_;++k)
	  {
	    double s 	= lcoo_[lcoo_ld_*k+1];
	    x[k] 		= s2*s - s1;
	  }
      
	for (wmesh_int_t j=0;j<=d_-i;++j)
	  {
	    bms_djacobip(2*i+1,
			 0,
			 j,		
			 lcoo_n_,
			 x,
			 1,
			 yj,
			 1,
			 work_n_,
			 work_);	  
	    //	  std::cout << "copy " << std::endl;
	    if (WMESH_STORAGE_BLOCK == v_storage_)
	      {
		for (wmesh_int_t k = 0; k < lcoo_n_;++k)
		  {
		    v_[v_ld_*idx+k] = yi[k] * yj[k];
		  }
	      }
	    else if (WMESH_STORAGE_INTERLEAVE == v_storage_)
	      {
		for (wmesh_int_t k = 0; k < lcoo_n_;++k)
		  {
		    v_[v_ld_*k+idx] = yi[k] * yj[k];
		  }
	      }
	    else
	      {
		WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
	      }
	    ++idx;
	  }

      
      }
    return WMESH_STATUS_SUCCESS;

  }
};

#if 0
wmesh_status_t bms_vandermonde_triangle(wmesh_int_t 	d_,
					wmesh_int_t 	lcoo_m_,
					wmesh_int_t 	lcoo_n_,
					const double * 	lcoo_,
					wmesh_int_t 	lcoo_ld_,
					wmesh_int_t 	v_storage_,
					double* 	v_,
					wmesh_int_t 	v_ld_,
					wmesh_int_t 	work_n_,
					double* 	work_)
{
  WMESH_CHECK(work_n_ >= 5 * lcoo_n_);
  static constexpr double s1 = 1.0;
  static constexpr double s2 = 2.0;
  
  double * x 	= work_ + 0*lcoo_n_;
  double * yi 	= work_ + 1*lcoo_n_;
  double * yj 	= work_ + 2*lcoo_n_;

  work_ 	= work_ + 3*lcoo_n_;
  work_n_	-=  3*lcoo_n_;

  wmesh_int_t idx = 0;
  for (wmesh_int_t i=0;i<=d_;++i)
    {
      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
	{
	  double r 	= lcoo_[lcoo_ld_*k+0];
	  double s 	= lcoo_[lcoo_ld_*k+1];
	  x[k] 		= (s < s1) ? (r * s2) / (s1 - s) - s1 : 0.0;
	}

      bms_djacobip(0,
		   0,
		   i,		
		   lcoo_n_,
		   x,
		   1,
		   yi,
		   1,
		   work_n_,
		   work_);
      
      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
	{
	  double s 	= lcoo_[lcoo_ld_*k+1];
	  yi[k] 	*= wmesh_math<double>::xpow(s1 - s, i);
	}
      
      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
	{
	  double s 	= lcoo_[lcoo_ld_*k+1];
	  x[k] 		= s2*s - s1;
	}
      
      for (wmesh_int_t j=0;j<=d_-i;++j)
	{
	  bms_djacobip(2*i+1,
		       0,
		       j,		
		       lcoo_n_,
		       x,
		       1,
		       yj,
		       1,
		       work_n_,
		       work_);	  
	  //	  std::cout << "copy " << std::endl;
	  if (WMESH_STORAGE_BLOCK == v_storage_)
	    {
	      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
		{
		  v_[v_ld_*idx+k] = yi[k] * yj[k];
		}
	    }
	  else if (WMESH_STORAGE_INTERLEAVE == v_storage_)
	    {
	      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
		{
		  v_[v_ld_*k+idx] = yi[k] * yj[k];
		}
	    }
	  else
	    {
	      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
	    }
	  ++idx;
	}

      
    }
  return WMESH_STATUS_SUCCESS;
}
#endif





template <wmesh_int_t ELEMENT, typename T>
wmesh_status_t bms_template_shape_vandermonde(wmesh_int_t 		d_,
					      wmesh_int_t 		ndofs_,
					      
					      wmesh_int_t 		lcoo_storage_,
					      wmesh_int_t 		lcoo_m_,
					      wmesh_int_t 		lcoo_n_,
					      const T * __restrict__	lcoo_,
					      wmesh_int_t 		lcoo_ld_,
					      
					      wmesh_int_t 		v_storage_,
					      T*__restrict__		v_,
					      wmesh_int_t 		v_ld_,
					      
					      wmesh_int_t 		rw_n_,
					      T*__restrict__		rw_)
{
  WMESH_CHECK_POINTER(lcoo_);
  wmesh_int_t num_points = (lcoo_storage_ == WMESH_STORAGE_INTERLEAVE) ? lcoo_n_ : lcoo_m_;  
  wmesh_status_t status;
  double *vandermonde_dofs = (double*)malloc(sizeof(double)*ndofs_ * ndofs_);
  double *vandermonde_eval = (double*)malloc(sizeof(double)*ndofs_ * num_points);

  //  WMESH_CHECK(work_n_ >= 5 * lcoo_n_);
  static constexpr double s1 = 1.0;
  static constexpr double s2 = 2.0;
  
  double * x 	= work_ + 0*lcoo_n_;
  double * yi 	= work_ + 1*lcoo_n_;
  double * yj 	= work_ + 2*lcoo_n_;  
  work_ 	= work_ + 3*lcoo_n_;
  work_n_      -= 3*lcoo_n_;
  
  wmesh_int_t idx = 0;
  for (wmesh_int_t i=0;i<=d_;++i)
    {
      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
	{
	  T r = lcoo_[lcoo_ld_*k+0];
	  T s = lcoo_[lcoo_ld_*k+1];
	  x[k] 		= (s < s1) ? (r * s2) / (s1 - s) - s1 : 0.0;
	}
      
      bms_djacobip(0,
		   0,
		   i,		
		   lcoo_n_,
		   x,
		   1,
		   yi,
		   1,
		   work_n_,
		   work_);
      
      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
	{
	  T s 	= lcoo_[lcoo_ld_*k+1];
	  yi[k] 	*= wmesh_math<double>::xpow(s1 - s, i);
	}
      
      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
	{
	  T s 	= lcoo_[lcoo_ld_*k+1];
	  x[k] 		= s2*s - s1;
	}
      
      for (wmesh_int_t j=0;j<=d_-i;++j)
	{
	  bms_djacobip(2*i+1,
		       0,
		       j,		
		       lcoo_n_,
		       x,
		       1,
		       yj,
		       1,
		       work_n_,
		       work_);	  
	  //	  std::cout << "copy " << std::endl;
	  if (WMESH_STORAGE_BLOCK == v_storage_)
	    {
	      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
		{
		  v_[v_ld_*idx+k] = yi[k] * yj[k];
		}
	    }
	  else if (WMESH_STORAGE_INTERLEAVE == v_storage_)
	    {
	      for (wmesh_int_t k = 0; k < lcoo_n_;++k)
		{
		  v_[v_ld_*k+idx] = yi[k] * yj[k];
		}
	    }
	  else
	    {
	      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
	    }
	  ++idx;
	}

      
    }
  return WMESH_STATUS_SUCCESS;


  

  wmesh_int_t rw_n = 0;
  switch(element_)
    {      
    case WMESH_ELEMENT_NODE:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
      }
      
    case WMESH_ELEMENT_EDGE:
      {
	rw_n = 2*ndofs;
	break;
      }
      
    case WMESH_ELEMENT_TRIANGLE:
      {
	rw_n 	= 5 * ndofs;
	break;
      }
      
    case WMESH_ELEMENT_QUADRILATERAL:
      {
	rw_n = 10*ndofs;
	break;
      }
      
    case WMESH_ELEMENT_TETRAHEDRON:
      {
	rw_n = 10*ndofs;
	break;
      }
    case WMESH_ELEMENT_PYRAMID:
      {
	rw_n = 10*ndofs;
	break;
      }
    case WMESH_ELEMENT_WEDGE:
      {
	rw_n = 10*ndofs;
	break;
      }
    case WMESH_ELEMENT_HEXAHEDRON:
      {
	rw_n = 10*ndofs;
	break;
      }
      
    default:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
	break;
      }

    }
  
  T *    rw   	= (double*)malloc(sizeof(double)*rw_n);
  if (!rw)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
    }

  wmesh_int_t rs_m   	= 2;
  wmesh_int_t rs_ld  	= 2;    
  wmesh_int_t rw_n 	= 5 * ndofs;  
  switch(element_)
    {      
    case WMESH_ELEMENT_NODE:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
	break;
      }
      
    case WMESH_ELEMENT_EDGE:
      {
	break;
      }
      
    case WMESH_ELEMENT_TRIANGLE:
      {
	status = bms_vandermonde_triangle(degree_,
					  rs_m,
					  ndofs_,
					  rs,
					  rs_ld,
					  1,
					  vandermonde_dofs,
					  ndofs_,
					  rw_n,
					  rw);
	WMESH_STATUS_CHECK(status);
	break;
      }

    case WMESH_ELEMENT_QUADRILATERAL:
      {
	status = bms_vandermonde_quadrilateral(degree_,
					       rs_m,
					       ndofs_,
					       rs,
					       rs_ld,
					       1,
					       vandermonde_dofs,
					       ndofs_,
					       rw_n,
					       rw);
	WMESH_STATUS_CHECK(status);
	break;
      }
      
    case WMESH_ELEMENT_TETRAHEDRON:
      {
	status = bms_vandermonde_tetrahedron(degree_,
					     rs_m,
					     ndofs_,
					     rs,
					     rs_ld,
					     1,
					     vandermonde_dofs,
					     ndofs_,
					     rw_n,
					     rw);
	break;
      }
      
    case WMESH_ELEMENT_PYRAMID:
      {
	status = bms_vandermonde_pyramid(degree_,
					 rs_m,
					 ndofs_,
					 rs,
					 rs_ld,
					 1,
					 vandermonde_dofs,
					 ndofs_,
					 rw_n,
					 rw);
	break;
      }
    case WMESH_ELEMENT_WEDGE:
      {
	status = bms_vandermonde_wedge(degree_,
				       rs_m,
				       ndofs_,
				       rs,
				       rs_ld,
				       1,
				       vandermonde_dofs,
				       ndofs_,
				       rw_n,
				       rw);
	break;
      }
    case WMESH_ELEMENT_HEXAHEDRON:
      {
	status = bms_vandermonde_hexahedron(degree_,
					    rs_m,
					    ndofs_,
					    rs,
					    rs_ld,
					    1,
					    vandermonde_dofs,
					    ndofs_,
					    rw_n,
					    rw);
	break;
      }
      
    default:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
	break;
      }
    }


  
  template <typename T>
    wmesh_status_t bms_template_shape_vandermonde(wmesh_int_t 		element_,
						  wmesh_int_t 		num_dofs_,
						  wmesh_int_t 		num_points_,
						  wmesh_int_p 		iw_n_,
						  wmesh_int_p 		rw_n_)
  {
    WMESH_CHECK_POINTER(iw_n_);
    WMESH_CHECK_POINTER(rw_n_);
    iw_n_[0] = ndofs;
    rw_n_[0] = ndofs_*ndofs_ + ndofs_*num_points + ndofs_*ndofs_;
    return WMESH_STATUS_SUCCESS;  
  }
 
 template <typename T>
   wmesh_status_t bms_template_shape_vandermonde(wmesh_int_t 		element_,
						 wmesh_int_t 		degree_,
						 wmesh_int_t 		ndofs_,
					      
						 wmesh_int_t 		dofcoo_storage_,
						 wmesh_int_t 		dofcoo_m_,
						 wmesh_int_t 		dofcoo_n_,
						 const T * __restrict__	dofcoo_,
						 wmesh_int_t 		dofcoo_ld_,
					      
						 wmesh_int_t 		lcoo_storage_,
						 wmesh_int_t 		lcoo_m_,
						 wmesh_int_t 		lcoo_n_,
						 const T * __restrict__	lcoo_,
						 wmesh_int_t 		lcoo_ld_,
					      
						 wmesh_int_t 		v_storage_,
						 T*__restrict__		v_,
						 wmesh_int_t 		v_ld_,

						 wmesh_int_t 		iw_n_,
						 wmesh_int_p		iw_,
						 wmesh_int_t 		rw_n_,
						 T*__restrict__		rw_)
 {
   WMESH_CHECK_POINTER(lcoo_);
   wmesh_int_t num_points = (lcoo_storage_ == WMESH_STORAGE_INTERLEAVE) ? lcoo_n_ : lcoo_m_;  
   wmesh_int_t num_dofs = (dofcoo_storage_ == WMESH_STORAGE_INTERLEAVE) ? dofcoo_n_ : dofcoo_m_;
   WMESH_CHECK(num_dofs == ndofs_);
  
   wmesh_status_t status;
   wmesh_int_t required_iw_n;
   wmesh_int_t required_rw_n;
   status =  bms_template_shape_vandermonde(element_,
					    num_dofs,
					    num_points,
					    &required_iw_n,
					    &required_rw_n);
   WMESH_STATUS_CHECK( status );
  
   if (required_rw_n > rw_n_)
     {
       WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
     }
   if (required_iw_n > iw_n_)
     {
       WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
     }

  
   wmesh_int_p perm 	= iw_;
   iw_ = perm + ndofs_;
   iw_n_ -= ndofs_;

   T *vandermonde_dofs 	= rw_;
   T *vandermonde_eval 	= rw_ + ndofs_ * ndofs_;
   T *id 			= rw_ + ndofs_ * ndofs_ + num_points * ndofs_;
   rw_ 				= id + ndofs_ * ndofs_ * num_points * ndofs_ + ndofs_ * ndofs_;
   rw_n_ 			-= ndofs_ * ndofs_  + num_points * ndofs_ + ndofs_ * ndofs_;
   for (wmesh_int_t i=0;i<ndofs*ndofs;++i) id[i] = 0.0;
   for (wmesh_int_t i=0;i<ndofs;++i) id[i*ndofs+i] = 1.0;
  

   bms_template_vandermonde(element_,
			    degree_,
			    ndofs_,
			    dofcoo_m_,
			    dofcoo_n_,
			    dofcoo_,
			    dofcoo_ld_,
			    WMESH_STORAGE_INTERLEAVE,
			    vandermonde_dofs,
			    ndofs,
			    rw_n_,
			    rw_);

   bms_template_vandermonde(element_,
			    degree_,
			    ndofs_,
			    lcoo_m_,
			    lcoo_n_,
			    lcoo_,
			    lcoo_ld_,
			    WMESH_STORAGE_INTERLEAVE,
			    vandermonde_eval,
			    ndofs,
			    rw_n_,
			    rw_);
  

   {
     wmesh_int_t info_lapack;  
     LAPACK_dgesv((wmesh_int_p)&ndofs,
		  (wmesh_int_p)&ev_n,
		  vandermonde_dofs,
		  (wmesh_int_p)&ndofs,
		  perm,
		  vandermonde_eval,
		  (wmesh_int_p)&ndofs,
		  (wmesh_int_p)&info_lapack);
     WMESH_CHECK( 0 == info_lapack );
   }
  
   for (int i=0;i<ndofs;++i)
     {
       for (int j=0;j<ev_n;++j)
	 {
	   std::cout << " " << vandermonde_eval[ndofs*j+i];	  
	 }
       std::cout << std::endl;
     }

 }


}

