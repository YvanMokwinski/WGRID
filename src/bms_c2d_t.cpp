#include <limits>
#include <iostream>
#include <array>
#include "wmesh_t.hpp"

//
// (x,y) = r*(1,0) + s*(0,1)                 // 0,1,2   x=r       y=s        
// (x,y) = (1-r-s)*(1,0) + r*(0,1) + s*(0,0) // 1,2,0   x=1-r-s   y=r       
// (x,y) = (1-r-s)*(0,1) + r*(0,0) + s*(1,0) // 2,0,1   x=x       y=1-r-s
// 
// 
// (x,y) = (1-r-s)*(0,0) + r*(0,1) + s*(1,0) // 0,2,1   x=s       y=r
// (x,y) = (1-r-s)*(0,1) + r*(1,0) + s*(0,0) // 2,1,0   x=r       y=1-r-s
// (x,y) = (1-r-s)*(1,0) + r*(0,0) + s*(0,1) // 1,0,2   x=1-r-s   y=s
//
static inline void  orientation_triangle_permutation(wmesh_int_t 	orientation_,
						     wmesh_int_t 	degree_,
						     wmesh_int_p	perm_)
{
  
  //  std::cout << "degree_ "<<degree_ << std::endl;
  const wmesh_int_t n = degree_-1;
  wmesh_int_t k = 0;
  for (wmesh_int_t i=1;i<=n;++i)
    {
      //      printf("i = %d\n",i);
      for (wmesh_int_t j=1;j<=n-i;++j)
	{
	  wmesh_int_t
	    r = 0,
	    s = 0;
	  switch(orientation_)
	    {
	    case -1:
	      {
		r = i;
		s = j;
		break;
	      }
	    case 1:
	      {
		r = j;
		s = i;
		break;
	      }
	    case -2:
	      {
		r = j;
		s = degree_ - i - j;
		break;
	      }
	    case 2:
	      {
		r = degree_ - i - j;
		s = j;
		break;
	      }
	    case -3:
	      {
		r = degree_ - i - j;
		s = i;
		break;
	      }
	    case 3:
	      {
		r = i;
		s = degree_- i - j;
		break;
	      }
	    }
	    
	  perm_[k++] = ((n-1)*n)/2 - ((n-r)*(n-r+1))/2 + (s-1);
	  //	  std::cout << "perm[" << k-1 << "] = " << perm_[k-1]<<std::endl;
	  //   3
	  //   2 6 
	  //   1 5 8
	  //   0 4 7 9 

	  //   4
	  //   3 8 
	  //   2 7 11
	  //   1 6 10 13
	  //   0 5 9  12 14
	}
    }
};

static inline wmesh_int_t orientation_triangle(wmesh_int_p icnc_,
					       wmesh_int_p jcnc_) noexcept
{
  const auto i0 = icnc_[0];
  const auto i1 = icnc_[1];
  const auto i2 = icnc_[2];
    
  const auto j0 = jcnc_[0];
  const auto j1 = jcnc_[1];
  const auto j2 = jcnc_[2];
    
  if ( (i0==j0) && (i1==j1) && (i2==j2) )
    {
      return 1;
    }    
  else if ( (i0==j1) && (i1==j2) && (i2==j0) )
    {
      return 2;
    }
  else if ( (i0==j2) && (i1==j0) && (i2==j1) )
    {
      return 3;
    }
  else if ( (i0==j0) && (i1==j2) && (i2==j1) )
    {
      return -1;
    }    
  else if ( (i0==j1) && (i1==j0) && (i2==j2) )
    {
      return -2;
    }
  else if ( (i0==j2) && (i1==j1) && (i2==j0) )
    {
      return -3;
    }
    
  return 0;
};
  
  



static  inline void get_c2t(wmesh_int_t 	m_,
			    const_wmesh_int_p 	c2t_,
			    wmesh_int_t 	c2t_ld_,
			    wmesh_int_t 	cell_idx_,
			    wmesh_int_p 	cnc_)
{
  for (wmesh_int_t i=0;i<m_;++i)
    {
      cnc_[i] = c2t_[cell_idx_*c2t_ld_+i];      
    }
};


static inline wmesh_status_t bms_c2d_t_calculate(wmesh_int_t 			c2n_ptr_,
									   wmesh_int_t 			c2n_m_,
								      wmesh_int_t 			c2n_n_,
								      const_wmesh_int_p		c2n_v_,
								      wmesh_int_t 			c2n_ld_,
								  
									   wmesh_int_t 			c2f_t_ptr_,
									   wmesh_int_t 			c2f_t_m_,
									   wmesh_int_t 			c2f_t_n_,
									   const_wmesh_int_p		c2f_t_v_,
									   wmesh_int_t 			c2f_t_ld_,
								  
									   wmesh_int_t 			c2d_t_ptr_,
									   wmesh_int_t 			c2d_t_m_,
									   wmesh_int_t 			c2d_t_n_,
									   wmesh_int_p			c2d_t_v_,
									   wmesh_int_t 			c2d_t_ld_,
									   
									   wmesh_int_t 			s_t2n_ptr_,
									   wmesh_int_t 			s_t2n_m_,
									   wmesh_int_t 			s_t2n_n_,
									   const_wmesh_int_p 		s_t2n_v_,
									   wmesh_int_t 			s_t2n_ld_,

									   wmesh_int_t 			num_triangles_,
									   wmesh_int_t 			degree_,
									   wmesh_int_t 			ndofs_per_triangle_,
									   wmesh_int_t 			dof_idx_)
{

  wmesh_int_t
    t2n[3],
    t2n_oriented[3],
    c2t[12],
    c2n[8];
  wmesh_int_p permloc = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*degree_*degree_);
  for (wmesh_int_t cell_idx = 0;cell_idx < c2n_n_;++cell_idx)
    {

      //
      // Extract nodes.
      //
      get_c2n(c2n_m_,
	      c2n_v_ + c2n_ptr_,
	      c2n_ld_,
	      cell_idx,
	      c2n);

      //
      // Extract triangle indexing.
      //
      get_c2t(c2f_t_m_,
	      c2f_t_v_ + c2f_t_ptr_,
	      c2f_t_ld_,
	      cell_idx,
	      c2t);
      
      //
      // Loop over the triangles.
      //
      for (wmesh_int_t ltriangle_idx = 0;ltriangle_idx<c2f_t_m_;++ltriangle_idx)
	{
	  wmesh_int_t triangle_idx = c2t[ltriangle_idx];	  
	  //std::cout << "triangle_idx " << triangle_idx << std::endl;
	  get_t2n(c2n,
		  ltriangle_idx,
		  t2n,
		  s_t2n_m_,
		  s_t2n_n_,
		  s_t2n_v_ + s_t2n_ptr_,
		  s_t2n_ld_);

	  for (int i=0;i<3;++i) t2n_oriented[i] = t2n[i];
	  for (int i=0;i<3;++i)
	    {
	      if ( (t2n_oriented[0] <t2n_oriented[1]) && (t2n_oriented[0] < t2n_oriented[2]) )
		{
		  break;
		}
	      else
		{
		  wmesh_int_t tmp = t2n_oriented[0];
		  t2n_oriented[0] = t2n_oriented[1];
		  t2n_oriented[1] = t2n_oriented[2];
		  t2n_oriented[2] = tmp;
		}
	    }
	  
	  wmesh_int_t orientation = orientation_triangle(t2n_oriented,
							 t2n);

	  if (t2n_oriented[1] > t2n_oriented[2])
	    {
	      orientation *= -1;
	    }
	  //	  std::cout << "degree "<< degree_ << std::endl;
	  //	  std::cout << "orientation "<< orientation << std::endl;
	  orientation_triangle_permutation	(orientation,
						 degree_,
						 permloc);
	  
	  for (wmesh_int_t ldof=0;ldof<ndofs_per_triangle_;++ldof)
	    {
	      // std::cout << permloc[ldof] << std::endl;
	      c2d_t_v_[c2d_t_ptr_ + cell_idx * c2d_t_ld_ + ltriangle_idx * ndofs_per_triangle_ + ldof]
		= dof_idx_ + triangle_idx * ndofs_per_triangle_ + permloc[ldof];
	    }
	}
    }
  free(permloc);
  return WMESH_STATUS_SUCCESS;

};


extern "C"
{
  wmesh_status_t  bms_c2d_t(wmesh_int_t 		c2n_size_,
						 const_wmesh_int_p 	c2n_ptr_,
						 const_wmesh_int_p 	c2n_m_,
						 const_wmesh_int_p 	c2n_n_,
						 const_wmesh_int_p 	c2n_v_,
						 const_wmesh_int_p 	c2n_ld_,
						 
						 wmesh_int_t 		c2f_t_size_,
						 const_wmesh_int_p 	c2f_t_ptr_,
						 const_wmesh_int_p 	c2f_t_m_,
						 const_wmesh_int_p 	c2f_t_n_,
						 const_wmesh_int_p 	c2f_t_v_,
						 const_wmesh_int_p 	c2f_t_ld_,
						 
						 wmesh_int_t 		c2d_t_size_,
						 const_wmesh_int_p 	c2d_t_ptr_,
						 const_wmesh_int_p 	c2d_t_m_,
						 const_wmesh_int_p 	c2d_t_n_,
						 wmesh_int_p 		c2d_t_v_,
						 const_wmesh_int_p 	c2d_t_ld_,

						 wmesh_int_t 		s_t2n_size_,
						 const_wmesh_int_p 	s_t2n_ptr_,
						 const_wmesh_int_p	s_t2n_m_,
						 const_wmesh_int_p	s_t2n_n_,
						 const_wmesh_int_p 	s_t2n_v_,
						 const_wmesh_int_p	s_t2n_ld_,
						      
						 wmesh_int_t 		num_triangles_,
						 wmesh_int_t 		degree_,
						 wmesh_int_t 		num_dofs_per_triangle_,
						 wmesh_int_t 		dof_idx_origin_)

  {
    if (0 == num_dofs_per_triangle_)
      {
	return WMESH_STATUS_SUCCESS;
      }
    if (0 == num_triangles_)
      {
	return WMESH_STATUS_SUCCESS;
      }
    
    WMESH_CHECK_POINTER(c2n_ptr_);    
    WMESH_CHECK_POINTER(c2n_m_);    
    WMESH_CHECK_POINTER(c2n_n_);
    WMESH_CHECK_POINTER(c2n_v_);
    WMESH_CHECK_POINTER(c2n_ld_);

    WMESH_CHECK_POINTER(c2f_t_ptr_);    
    WMESH_CHECK_POINTER(c2f_t_m_);    
    WMESH_CHECK_POINTER(c2f_t_n_);
    WMESH_CHECK_POINTER(c2f_t_v_);
    WMESH_CHECK_POINTER(c2f_t_ld_);

    WMESH_CHECK_POINTER(c2d_t_ptr_);    
    WMESH_CHECK_POINTER(c2d_t_m_);    
    WMESH_CHECK_POINTER(c2d_t_n_);
    WMESH_CHECK_POINTER(c2d_t_v_);
    WMESH_CHECK_POINTER(c2d_t_ld_);

    WMESH_CHECK_POINTER(s_t2n_ptr_);    
    WMESH_CHECK_POINTER(s_t2n_m_);    
    WMESH_CHECK_POINTER(s_t2n_n_);
    WMESH_CHECK_POINTER(s_t2n_v_);
    WMESH_CHECK_POINTER(s_t2n_ld_);
    
    for (wmesh_int_t cell_type=0;cell_type<c2n_size_;++cell_type)
      {	
	wmesh_status_t status = bms_c2d_t_calculate(c2n_ptr_[cell_type],
									 c2n_m_[cell_type],
									 c2n_n_[cell_type],
									 c2n_v_,
									 c2n_ld_[cell_type],
									 
									 c2f_t_ptr_[cell_type],
									 c2f_t_m_[cell_type],
									 c2f_t_n_[cell_type],
									 c2f_t_v_,
									 c2f_t_ld_[cell_type],
									 
									 c2d_t_ptr_[cell_type],
									 c2d_t_m_[cell_type],
									 c2d_t_n_[cell_type],
									 c2d_t_v_,
									 c2d_t_ld_[cell_type],
									 
									 s_t2n_ptr_[cell_type],
									 s_t2n_m_[cell_type],
									 s_t2n_n_[cell_type],
									 s_t2n_v_,
									 s_t2n_ld_[cell_type],

									 num_triangles_,
									 degree_,
									 num_dofs_per_triangle_,
									 dof_idx_origin_);
	WMESH_STATUS_CHECK(status);
      }
    return WMESH_STATUS_SUCCESS;
  }
};
