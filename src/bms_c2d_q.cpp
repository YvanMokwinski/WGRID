#include <limits>
#include <iostream>
#include <array>
#include "wmesh_t.hpp"


static inline void  orientation_quadrilateral_permutation(wmesh_int_t 	orientation_,
							  wmesh_int_t	degree_,
							  wmesh_int_p	perm_) noexcept
{    
  const wmesh_int_t n = degree_-1;
  wmesh_int_t k = 0;
  for (wmesh_int_t i=1;i<=n;++i)
    {
      for (wmesh_int_t j=1;j<=n;++j)
	{
	  wmesh_int_t r = 0;
	  wmesh_int_t s = 0;
	  switch(orientation_)
	    {
	      /* on a interchange 2 et -2 attends toi a faire de meme avec 3 et 4*/
	    case 1:
	      {
		r = j;
		s = i;
		/* (j,i) 
		   0,3,2,1
		*/
		break;
	      }
	    case -2:
	      {
		/* (j,k-i)
		   1,2,3,0
		*/
		r = j;
		s = degree_-i;		      
		break;
	      }
	    case -3:
	      {
		/* (k-i,k-j) 
		   2,3,0,1
		*/
		r = degree_-i;
		s = degree_-j;
		break;
	      }
	    case -4:
	      {
		r = degree_-j;
		s = i;
		/* (k-j,i) 
		   3,0,1,2*/
		break;
	      }
		
	    case -1:
	      {
		/* (i,j)
		   0,1,2,3
		*/
		r = i;
		s = j;
		break;
	      }
		
	    case 2:
	      {
		r = degree_-i;
		s = j;
		/* (k-i,j)
		   1,0,3,2
		*/
		break;
	      }	      
	    case 3:
	      {
		r = degree_-j;
		s = degree_-i;
		/* (k-j,k-i)
		   2,1,0,3
		*/
		break;
	      }
	    case 4:
	      {
		r = i;
		s = degree_-j;
		/* (i,k-j)
		   3,2,1,0
		*/
		break;
	      }
		
	    }
	    
	  perm_[k++] = n*(r-1)+(s-1);
	} 	
    }     
};
  
static inline wmesh_int_t orientation_quadrilateral(const_wmesh_int_p icnc_,
						    const_wmesh_int_p jcnc_)
{
#if 0
  static const int_t facequad_orientation[]={0,1,2,3,/*1*/
					     1,2,3,0,/*2*/
					     2,3,0,1,/*3*/
					     3,0,1,2,/*4*/
					     
					     0,3,2,1,/*-1*/
					     1,0,3,2,/*-2*/
					     2,1,0,3,/*-3*/
					     3,2,1,0/*-4*/};
#endif

  
  const auto i0 = icnc_[0];
  const auto i1 = icnc_[1];
  const auto i2 = icnc_[2];
  const auto i3 = icnc_[3];
  
  const auto j0 = jcnc_[0];
  const auto j1 = jcnc_[1];
  const auto j2 = jcnc_[2];
  const auto j3 = jcnc_[3];
  
  if ( (i0==j0) && (i1==j1) && (i2==j2) && (i3==j3) )
    {
      return 1;
    }    
  else if ( (i0==j1) && (i1==j2) && (i2==j3) && (i3==j0) )
    {
      return 2;
    }
  else if ( (i0==j2) && (i1==j3) && (i2==j0) && (i3==j1) )
    {
      return 3;
    }
  else if ( (i0==j3) && (i1==j0) && (i2==j1) && (i3==j2) )
    {
      return 4;
    }  
  else if ( (i0==j0) && (i1==j3) && (i2==j2) && (i3==j1) )
    {
      return -1;
    }
  else if ( (i0==j1) && (i1==j0) && (i2==j3) && (i3==j2) )
    {
      return -2;
    }
  else if ( (i0==j2) && (i1==j1) && (i2==j0) && (i3==j3) )
    {
      return -3;
    }
  else if ( (i0==j3) && (i1==j2) && (i2==j1) && (i3==j0) )
    {
      return -4;
    }
  
  return 0;
};



static  inline void get_c2q(wmesh_int_t 	m_,
			    const_wmesh_int_p 	c2q_,
			    wmesh_int_t 	c2q_ld_,
			    wmesh_int_t 	cell_idx_,
			    wmesh_int_p 	cnc_)
{
  for (wmesh_int_t i=0;i<m_;++i)
    {
      cnc_[i] = c2q_[cell_idx_*c2q_ld_+i];      
    }
};


static inline wmesh_status_t bms_c2d_q_calculate(wmesh_int_t 			c2n_ptr_,
									   wmesh_int_t 			c2n_m_,
									   wmesh_int_t 			c2n_n_,
									   const_wmesh_int_p		c2n_v_,
									   wmesh_int_t 			c2n_ld_,
								  
									   wmesh_int_t 			c2f_q_ptr_,
									   wmesh_int_t 			c2f_q_m_,
									   wmesh_int_t 			c2f_q_n_,
									   const_wmesh_int_p		c2f_q_v_,
									   wmesh_int_t 			c2f_q_ld_,
								  
									   wmesh_int_t 			c2d_q_ptr_,
									   wmesh_int_t 			c2d_q_m_,
									   wmesh_int_t 			c2d_q_n_,
									   wmesh_int_p			c2d_q_v_,
									   wmesh_int_t 			c2d_q_ld_,
									   
									   wmesh_int_t 			s_q2n_ptr_,
									   wmesh_int_t 			s_q2n_m_,
									   wmesh_int_t 			s_q2n_n_,
									   const_wmesh_int_p 		s_q2n_v_,
									   wmesh_int_t 			s_q2n_ld_,

									   wmesh_int_t 			num_quadrilaterals_,
									   wmesh_int_t 			degree_,
									   wmesh_int_t 			ndofs_per_quadrilateral_,
									   wmesh_int_t 			dof_idx_)
{

  wmesh_int_t
    q2n[4],
    q2n_oriented[4],
    c2q[12],
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
      // Extract quadrilateral indexing.
      //
      get_c2q(c2f_q_m_,
	      c2f_q_v_ + c2f_q_ptr_,
	      c2f_q_ld_,
	      cell_idx,
	      c2q);
      
      //
      // Loop over the quadrilaterals.
      //
      for (wmesh_int_t lquadrilateral_idx = 0;lquadrilateral_idx<c2f_q_m_;++lquadrilateral_idx)
	{
	  wmesh_int_t quadrilateral_idx = c2q[lquadrilateral_idx];	  

	  get_q2n(c2n,
		  lquadrilateral_idx,
		  q2n,
		  s_q2n_m_,
		  s_q2n_n_,
		  s_q2n_v_ + s_q2n_ptr_,
		  s_q2n_ld_);

	  for (int i=0;i<4;++i) q2n_oriented[i] = q2n[i];
	  for (int i=0;i<4;++i)
	    {
	      if ( (q2n_oriented[0] <q2n_oriented[1]) && (q2n_oriented[0] < q2n_oriented[2]) && (q2n_oriented[0] < q2n_oriented[3]) )
		{
		  break;
		}
	      else
		{
		  wmesh_int_t tmp = q2n_oriented[0];
		  q2n_oriented[0] = q2n_oriented[1];
		  q2n_oriented[1] = q2n_oriented[2];
		  q2n_oriented[2] = q2n_oriented[3];
		  q2n_oriented[3] = tmp;
		}
	    }
#if 0
	  // 2 3 7 6
	  std::cout << "refquad "
		    << q2n_oriented[0]
		    << " " 
		    << q2n_oriented[1]
		    << " " 
		    << q2n_oriented[2] 
		    << " "
		    << q2n_oriented[3] 
		    << " "
		    << std::endl;

	  std::cout << "quad "
		    << q2n[0]
		    << " " 
		    << q2n[1]
		    << " " 
		    << q2n[2] 
		    << " "
		    << q2n[3] 
		    << " "
		    <<std::endl;
#endif
	  wmesh_int_t orientation = orientation_quadrilateral(q2n_oriented,
							      q2n);

	  if (q2n_oriented[1] > q2n_oriented[3])
	    orientation *= -1;
	  
#if 0
	  std::cout << "orientation "
		    << orientation
		    << " " <<std::endl;
#endif	  
	  
	  orientation_quadrilateral_permutation(orientation,
						degree_,
						permloc);
	  //  printf("ddd %d %d\n",quadrilateral_idx,ndofs_per_quadrilateral_);
	  for (wmesh_int_t ldof=0;ldof<ndofs_per_quadrilateral_;++ldof)
	    {
	      c2d_q_v_[c2d_q_ptr_ + cell_idx * c2d_q_ld_ + lquadrilateral_idx * ndofs_per_quadrilateral_ + ldof]
		= dof_idx_ + quadrilateral_idx * ndofs_per_quadrilateral_ + permloc[ldof];
	    }
	}
    }
  free(permloc);
  return WMESH_STATUS_SUCCESS;

};


extern "C"
{
  wmesh_status_t  bms_c2d_q(wmesh_int_t 		c2n_size_,					     
						      const_wmesh_int_p 	c2n_ptr_,
						      const_wmesh_int_p 	c2n_m_,
						      const_wmesh_int_p 	c2n_n_,
						      const_wmesh_int_p 	c2n_v_,
						      const_wmesh_int_p 	c2n_ld_,
						 
						      wmesh_int_t 		c2f_q_size_,					     
						      const_wmesh_int_p 	c2f_q_ptr_,
						      const_wmesh_int_p 	c2f_q_m_,
						      const_wmesh_int_p 	c2f_q_n_,
						      const_wmesh_int_p 	c2f_q_v_,
						      const_wmesh_int_p 	c2f_q_ld_,
						 
						      wmesh_int_t 		c2d_q_size_,					     
						      const_wmesh_int_p 	c2d_q_ptr_,
						      const_wmesh_int_p 	c2d_q_m_,
						      const_wmesh_int_p 	c2d_q_n_,
						      wmesh_int_p 		c2d_q_v_,
						      const_wmesh_int_p 	c2d_q_ld_,
						 
						      wmesh_int_t 		s_q2n_size_,					     
						      const_wmesh_int_p 	s_q2n_ptr_,
						      const_wmesh_int_p		s_q2n_m_,
						      const_wmesh_int_p		s_q2n_n_,
						      const_wmesh_int_p 	s_q2n_v_,
						      const_wmesh_int_p		s_q2n_ld_,
						      
						      wmesh_int_t 		num_quadrilaterals_,
						      wmesh_int_t 		degree_,
						      wmesh_int_t 		num_dofs_per_quadrilateral_,
						      wmesh_int_t 		dof_idx_origin_)

  {
    if (0 == num_dofs_per_quadrilateral_)
      {
	return WMESH_STATUS_SUCCESS;
      }
    if (0 == num_quadrilaterals_)
      {
	return WMESH_STATUS_SUCCESS;
      }
    
    WMESH_CHECK_POINTER(c2n_ptr_);    
    WMESH_CHECK_POINTER(c2n_m_);    
    WMESH_CHECK_POINTER(c2n_n_);
    WMESH_CHECK_POINTER(c2n_v_);
    WMESH_CHECK_POINTER(c2n_ld_);

    WMESH_CHECK_POINTER(c2f_q_ptr_);    
    WMESH_CHECK_POINTER(c2f_q_m_);    
    WMESH_CHECK_POINTER(c2f_q_n_);
    WMESH_CHECK_POINTER(c2f_q_v_);
    WMESH_CHECK_POINTER(c2f_q_ld_);

    WMESH_CHECK_POINTER(c2d_q_ptr_);    
    WMESH_CHECK_POINTER(c2d_q_m_);    
    WMESH_CHECK_POINTER(c2d_q_n_);
    WMESH_CHECK_POINTER(c2d_q_v_);
    WMESH_CHECK_POINTER(c2d_q_ld_);
    
    WMESH_CHECK_POINTER(s_q2n_ptr_);    
    WMESH_CHECK_POINTER(s_q2n_m_);    
    WMESH_CHECK_POINTER(s_q2n_n_);
    WMESH_CHECK_POINTER(s_q2n_v_);
    WMESH_CHECK_POINTER(s_q2n_ld_);
    
    for (wmesh_int_t cell_type=0;cell_type<c2n_size_;++cell_type)
      {	
	wmesh_status_t status = bms_c2d_q_calculate(c2n_ptr_[cell_type],
									      c2n_m_[cell_type],
									      c2n_n_[cell_type],
									      c2n_v_,
									      c2n_ld_[cell_type],
									 
									      c2f_q_ptr_[cell_type],
									      c2f_q_m_[cell_type],
									      c2f_q_n_[cell_type],
									      c2f_q_v_,
									      c2f_q_ld_[cell_type],
									 
									      c2d_q_ptr_[cell_type],
									      c2d_q_m_[cell_type],
									      c2d_q_n_[cell_type],
									      c2d_q_v_,
									      c2d_q_ld_[cell_type],
									 
									      s_q2n_ptr_[cell_type],
									      s_q2n_m_[cell_type],
									      s_q2n_n_[cell_type],
									      s_q2n_v_,
									      s_q2n_ld_[cell_type],

									      num_quadrilaterals_,
									      degree_,
									      num_dofs_per_quadrilateral_,
									      dof_idx_origin_);
	WMESH_STATUS_CHECK(status);
      }
    return WMESH_STATUS_SUCCESS;
  }
};
