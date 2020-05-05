#include <stdlib.h>
#include <iostream>
#include "wmesh.hpp"

static inline wmesh_int_t a_max4(const_wmesh_int_p s_,wmesh_int_p p_)
{
  wmesh_int_t mx = s_[0];
  p_[0] = 0;
  if (mx < s_[1]) {mx = s_[1];p_[0] = 1;}
  if (mx < s_[2]) {mx = s_[2];p_[0] = 2;}
  if (mx < s_[3]) {mx = s_[3];p_[0] = 3;}
  return mx;
}

static inline wmesh_int_t a_min4(const_wmesh_int_p s_,wmesh_int_p p_)
{
  wmesh_int_t mx = s_[0];
  p_[0] = 0;
  if (mx > s_[1]) {mx = s_[1];p_[0] = 1;}
  if (mx > s_[2]) {mx = s_[2];p_[0] = 2;}
  if (mx > s_[3]) {mx = s_[3];p_[0] = 3;}
  return mx;
}

static inline wmesh_int_t a_max4(const_wmesh_int_p s_)
{
  wmesh_int_t mx = s_[0];
  if (mx < s_[1]) mx = s_[1];
  if (mx < s_[2]) mx = s_[2];
  if (mx < s_[3]) mx = s_[3];
  return mx;
}

static inline wmesh_int_t a_min4(const_wmesh_int_p s_)
{
  wmesh_int_t mx = s_[0];
  if (mx > s_[1]) mx = s_[1];
  if (mx > s_[2]) mx = s_[2];
  if (mx > s_[3]) mx = s_[3];
  return mx;
}

static wmesh_status_t extrude_tetra(wmesh_int_t		degree_,				   
				    const_wmesh_int_p 	tr2n_,
				    wmesh_int_t 	tr2n_ld_,
				    wmesh_int_p 	tet2n_,
				    wmesh_int_t 	tet2n_ld_)
{
  WMESH_POINTER_CHECK(tr2n_);
  WMESH_POINTER_CHECK(tet2n_);
  const wmesh_int_t s_numTriangles = (degree_ + 1)*(degree_ + 1);
  wmesh_int_t count	= 0;
  wmesh_int_t dec0 	= 0;
  wmesh_int_t dec1 	= ((degree_ + 1)*(degree_ + 2))/2;  
  wmesh_int_t N 	= s_numTriangles;
  
  for (wmesh_int_t itvl=0;itvl<degree_;++itvl)
    {
      N=(degree_-itvl-1)*(degree_-itvl-1);
      for (wmesh_int_t i=0;i<N;++i)
	{

	  
	  wmesh_int_t mx = a_min4(&tr2n_[i*tr2n_ld_ + 0]);

	  tet2n_[count*tet2n_ld_ + 0] = dec0 + tr2n_[i*tr2n_ld_ + 0];
	  tet2n_[count*tet2n_ld_ + 1] = dec0 + tr2n_[i*tr2n_ld_ + 1];
	  tet2n_[count*tet2n_ld_ + 2] = dec0 + tr2n_[i*tr2n_ld_ + 2];
	  tet2n_[count*tet2n_ld_ + 3] = dec1 + mx;
	  ++count;

	  wmesh_int_t mx2 = a_max4(&tr2n_[i*tr2n_ld_ + 0]);
	  
	  tet2n_[count*tet2n_ld_ + 0] = dec1+tr2n_[i*tr2n_ld_ + 0];
	  tet2n_[count*tet2n_ld_ + 2] = dec1+tr2n_[i*tr2n_ld_ + 1];
	  tet2n_[count*tet2n_ld_ + 1] = dec1+tr2n_[i*tr2n_ld_ + 2];
	  tet2n_[count*tet2n_ld_ + 3] = dec0+mx2;
	  ++count;
	} 
	

      for (wmesh_int_t i=0;i<N;++i)
	{
	  wmesh_int_t j1,j2,j3;
	  wmesh_int_t mx1 = a_min4(&tr2n_[i*tr2n_ld_ + 0],&j1);
	  wmesh_int_t mx2 = a_max4(&tr2n_[i*tr2n_ld_ + 0],&j2);

	  j3 = -1;
	  if ((j1==0)&&(j2==1))
	    {
	      j3 = 2;
	    }
	  else if ((j1==1)&&(j2==2))
	    {
	      j3 = 0;
	    }
	  else if ((j1==0)&&(j2==2))
	    {
	      j3 = 1;
	    }
	  else if ((j1==1)&&(j2==0))
	    {
	      j3 = 2;
	    }
	  else if ((j1==2)&&(j2==1))
	    {
	      j3 = 0;
	    }
	  else if ((j1==2)&&(j2==0))
	    {
	      j3 = 1;
	    }

	  tet2n_[count*tet2n_ld_ + 0] = dec0 + tr2n_[i*tr2n_ld_ + j3];
	  tet2n_[count*tet2n_ld_ + 1] = dec0 + mx2;
	  tet2n_[count*tet2n_ld_ + 2] = dec1 + tr2n_[i*tr2n_ld_ + j3];
	  tet2n_[count*tet2n_ld_ + 3] = dec1 + mx1;
	  if (tet2n_[count*tet2n_ld_ + 0]!=tet2n_[count*tet2n_ld_ + 1]-1)
	    {
	      tet2n_[count*tet2n_ld_ + 0] = dec0+tr2n_[i*tr2n_ld_ + j3];
	      tet2n_[count*tet2n_ld_ + 2] = dec0+mx2;
	      tet2n_[count*tet2n_ld_ + 1] = dec1+tr2n_[i*tr2n_ld_ + j3];
	      tet2n_[count*tet2n_ld_ + 3] = dec1+mx1;		  
	    }
	  ++count;
	} 

      wmesh_int_t n = N;
      N = (degree_-itvl)*(degree_-itvl);		
      for (wmesh_int_t i=n;i<N;++i)
	{
	  wmesh_int_t mx = a_min4(&tr2n_[i*tr2n_ld_ + 0]);
	  tet2n_[count*tet2n_ld_ + 0] = dec0+tr2n_[i*tr2n_ld_ + 0];
	  tet2n_[count*tet2n_ld_ + 1] = dec0+tr2n_[i*tr2n_ld_ + 1];
	  tet2n_[count*tet2n_ld_ + 2] = dec0+tr2n_[i*tr2n_ld_ + 2];
	  tet2n_[count*tet2n_ld_ + 3] = dec1+mx;
	  ++count;
	} 

      for (wmesh_int_t i=n+1;i<N;i+=2)
	{
	  wmesh_int_t j1,j2,j3;
	  wmesh_int_t mx1 = a_min4(&tr2n_[i*tr2n_ld_ + 0],&j1);
	  wmesh_int_t mx2 = a_max4(&tr2n_[i*tr2n_ld_ + 0],&j2);
	  j3=-1;
	  if ((j1==0)&&(j2==1))
	    j3 = 2;
	  else if ((j1==1)&&(j2==2))
	    j3 = 0;
	  else if ((j1==0)&&(j2==2))
	    j3 = 1;
	  else if ((j1==1)&&(j2==0))
	    j3 = 2;
	  else if ((j1==2)&&(j2==1))
	    j3 = 0;
	  else if ((j1==2)&&(j2==0))
	    j3 = 1;
	  tet2n_[count*tet2n_ld_ + 0] = dec0+tr2n_[i*tr2n_ld_ + j3];
	  tet2n_[count*tet2n_ld_ + 1] = dec0+mx2;
	  tet2n_[count*tet2n_ld_ + 2] = dec1+tr2n_[i*tr2n_ld_ + j3];
	  tet2n_[count*tet2n_ld_ + 3] = dec1+mx1;
	  if (tet2n_[count*tet2n_ld_ + 0]!=tet2n_[count*tet2n_ld_ + 1]-1)
	    {		
	      tet2n_[count*tet2n_ld_ + 0] = dec0+tr2n_[i*tr2n_ld_ + j3];
	      tet2n_[count*tet2n_ld_ + 2] = dec0+mx2;
	      tet2n_[count*tet2n_ld_ + 1] = dec1+tr2n_[i*tr2n_ld_ + j3];
	      tet2n_[count*tet2n_ld_ + 3] = dec1+mx1;		
	    }
	  ++count;
	} 
      
      dec0 = dec1;
      dec1 = dec1+((degree_-itvl)*(degree_-itvl+1))/2;
    }
  return WMESH_STATUS_SUCCESS;
};

static inline wmesh_status_t wmesh_treilli_calculate_c2n(wmesh_int_t 		cell_type_,
							 wmesh_int_t 		degree_,
							 const_wmesh_int_p 	c2n_ptr_,
							 const_wmesh_int_p 	c2n_m_,
							 const_wmesh_int_p 	c2n_n_,
							 wmesh_int_p 		c2n_v_,
							 const_wmesh_int_p 	c2n_ld_,
							 wmesh_int_t 		work_n_,
							 wmesh_int_p 		work_)
{
  WMESH_POINTER_CHECK(c2n_ptr_);
  WMESH_POINTER_CHECK(c2n_m_);
  WMESH_POINTER_CHECK(c2n_n_);
  WMESH_POINTER_CHECK(c2n_v_);
  WMESH_POINTER_CHECK(c2n_ld_);
  
  if (cell_type_ == 0 && work_n_ < (degree_+1)*(degree_+1) * 3)
    {
      return WMESH_STATUS_ERROR_WORKSPACE;      
    }
  
  WMESH_POINTER_CHECK(c2n_v_);

  wmesh_status_t status;
  
  //
  // Tetrahedra.
  //
  if (cell_type_==0)
    {
      //
      // Needs workspace.
      //
      WMESH_POINTER_CHECK(work_);      
#define _dec(_i,_j) ((  (_i)+(_j) + 1 )*( (_i)+(_j) ))/2+(_i)
      { wmesh_int_p p 	= work_;      
	for (wmesh_int_t j=0;j<(degree_+1);j++)
	  {	  
	    for (wmesh_int_t i=0;i<j;i++)
	      {
		
		*p++ 	= _dec(i,j-i);
		*p++ 	= _dec(i+1,j-i);
		*p++ 	= _dec(i,j+1-i); 
		
		*p++  	= _dec(i,j-i);
		*p++ 	= _dec(i+1,j-1-i);
		*p++ 	= _dec(i+1,j-i);
	      }
	    
	    *p++	=  _dec(j,0);
	    *p++	=  _dec(j+1,0);
	    *p++	=  _dec(j,1);	      
	  } }      
#undef _dec
      
      status  = extrude_tetra(degree_,			     
			      work_,
			      3,
			      c2n_v_ + c2n_ptr_[cell_type_],
			      c2n_ld_[cell_type_]);
      
      WMESH_STATUS_CHECK(status);
      
      return WMESH_STATUS_SUCCESS;
    }
  
  if (cell_type_==1)
    {
      wmesh_int_p c2n_t		= c2n_v_ + c2n_ptr_[0];
      wmesh_int_t c2n_t_shift 	= c2n_ld_[0];
      wmesh_int_p c2n_p 	= c2n_v_ + c2n_ptr_[cell_type_];
      wmesh_int_t c2n_p_shift 	= c2n_ld_[cell_type_];
      
#define _dec(_i,_j,_k) (degree_+1)*(degree_+1)*(_k) + (degree_+1) * (_i) + (_j)
      
      for (wmesh_int_t k=0;k<degree_;++k)
	{
	  //
	  // Pyramids
	  //
	  {
	    for (wmesh_int_t i=0;i<degree_-k;++i)
	      {
		for (wmesh_int_t j=0;j<degree_-k;++j)
		  {
		    c2n_p[0] = _dec(i,j,k);
		    c2n_p[1] = _dec(i+1,j,k);
		    c2n_p[2] = _dec(i+1,j+1,k);
		    c2n_p[3] = _dec(i,j+1,k);
		    c2n_p[4] = _dec(i,j,k+1);
		    c2n_p += c2n_p_shift;
		  }
	      }
	    for (wmesh_int_t i=0;i<degree_-1-k;++i)
	      {
		for (wmesh_int_t j=0;j<degree_-1-k;++j)
		  {
		    c2n_p[0] = _dec(i,j,k+1);
		    c2n_p[1] = _dec(i,j+1,k+1);
		    c2n_p[2] = _dec(i+1,j+1,k+1);
		    c2n_p[3] = _dec(i+1,j,k+1);
		    c2n_p[4] = _dec(i+1,j+1,k);
		    c2n_p += c2n_p_shift;
		  }
	      }
	  }
	  
	  //
	  // Tetrahedra
	  //
	  {
	    for (wmesh_int_t i=0;i<degree_-1-k;++i)
	      {
		for (wmesh_int_t j=0;j<degree_-k;++j)
		  {
		    c2n_t[0] = _dec(i+1,j,k);
		    c2n_t[1] = _dec(i+1,j+1,k);
		    c2n_t[2] = _dec(i,j,k+1);
		    c2n_t[3] = _dec(i+1,j,k+1);
		    c2n_t += c2n_t_shift;
		  }
	      }

	    for (wmesh_int_t i=0;i<degree_-k;++i)
	      {
		for (wmesh_int_t j=0;j<degree_-1-k;++j)
		  {
		    c2n_t[0] = _dec(i,j+1,k);
		    c2n_t[1] = _dec(i+1,j+1,k);
		    c2n_t[2] = _dec(i,j+1,k+1);
		    c2n_t[3] = _dec(i,j,k+1);
		    c2n_t += c2n_t_shift;
		  }
	      }
	  }
	  
	}
#undef _dec
      return WMESH_STATUS_SUCCESS;
    }
  
  if (cell_type_==2)
    {
#define _dec(_i,_j,_k) (_k)* ( ((degree_+1)*(degree_+2)) /2) + ((  (_i)+(_j) + 1 )*( (_i)+(_j) ))/2 + (_i)
      wmesh_int_p c2n 		= c2n_v_ + c2n_ptr_[cell_type_];
      wmesh_int_t c2n_shift 	= c2n_ld_[cell_type_];
      for (wmesh_int_t k=0;k<degree_;++k)
	{
	  for (wmesh_int_t j=0;j<degree_;++j)
	    {
	      for (wmesh_int_t i=0;i<j;++i)
		{
		  c2n[0] 	= _dec(i,j-i,k);
		  c2n[1] 	= _dec(i+1,j-i,k);
		  c2n[2] 	= _dec(i,j+1-i,k); 
		  c2n[3] 	= _dec(i,j-i,k+1);
		  c2n[4] 	= _dec(i+1,j-i,k+1);
		  c2n[5] 	= _dec(i,j+1-i,k+1); 
		  c2n += c2n_shift;
		  
		  c2n[0] 	= _dec(i,j-i,k);
		  c2n[1] 	= _dec(i+1,j-1-i,k);
		  c2n[2] 	= _dec(i+1,j-i,k);
		  c2n[3] 	= _dec(i,j-i,k+1);
		  c2n[4] 	= _dec(i+1,j-1-i,k+1);
		  c2n[5] 	= _dec(i+1,j-i,k+1);
		  c2n += c2n_shift;
		} 
	      
	      c2n[0]		=  _dec(j,0,k);
	      c2n[1]		=  _dec(j+1,0,k);
	      c2n[2]		=  _dec(j,1,k);	      
	      c2n[3]		=  _dec(j,0,k+1);
	      c2n[4]		=  _dec(j+1,0,k+1);
	      c2n[5]		=  _dec(j,1,k+1);	      
	      c2n += c2n_shift;
	    }
	}
#undef _dec
  
      return WMESH_STATUS_SUCCESS;
    }
  
  if (cell_type_==3)
    {      
#define _dec(_i,_j,_k) (degree_+1)*(degree_+1)*(_k) + (degree_+1) * (_i) + (_j)
      wmesh_int_p c2n 		= c2n_v_ + c2n_ptr_[cell_type_];
      wmesh_int_t c2n_shift 	= c2n_ld_[cell_type_];
      for (wmesh_int_t k=0;k<degree_;++k)
	{
	  for (wmesh_int_t i=0;i<degree_;i++)
	    {
	      for (wmesh_int_t j=0;j<degree_;j++)
		{
		  c2n[0] = _dec(i,j,k);
		  c2n[1] = _dec(i+1,j,k);
		  c2n[2] = _dec(i+1,j+1,k);
		  c2n[3] = _dec(i,j+1,k);
		  
		  c2n[4] = _dec(i,j,k+1);
		  c2n[5] = _dec(i+1,j,k+1);
		  c2n[6] = _dec(i+1,j+1,k+1);
		  c2n[7] = _dec(i,j+1,k+1);
#if 0
		  c2n[0] = _dec(i+1,j+1,k);
		  c2n[1] = _dec(i,j+1,k);
		  c2n[2] = _dec(i,j,k);
		  c2n[3] = _dec(i+1,j,k);
		  c2n[4] = _dec(i+1,j+1,k+1);
		  c2n[5] = _dec(i,j+1,k+1);
		  c2n[6] = _dec(i,j,k+1);
		  c2n[7] = _dec(i+1,j,k+1);
#endif
		  c2n += c2n_shift;
		} 
	    }
	}
#undef _dec
      return WMESH_STATUS_SUCCESS;
    }

  return WMESH_STATUS_INVALID_ARGUMENT;
}


static inline wmesh_status_t wmesh_treilli_calculate_icoo(wmesh_int_t cell_type_,
							  wmesh_int_t degree_,
							  wmesh_int_p icoo_,
							  wmesh_int_t icoo_ld_) 
{
  WMESH_POINTER_CHECK(icoo_);
  wmesh_int_t n1d = degree_+1;
  //
  // Tetra.
  //
  if (cell_type_==0)
    {
      wmesh_int_t idx = 0;
      for (wmesh_int_t i=0;i<n1d;i++)
	{
	  for (wmesh_int_t j=0;j<=i;j++)
	    {
	      icoo_[icoo_ld_*idx+0] = j;
	      icoo_[icoo_ld_*idx+1] = i-j;
	      icoo_[icoo_ld_*idx+2] = 0;
	      ++idx;
	    } 
	}
      
      { wmesh_int_t N     = ((degree_+1)*(degree_+2))/2-(degree_+1);
	for (wmesh_int_t irot=1;irot<n1d;++irot)
	  {
	    for (wmesh_int_t i=0;i<N;++i)
	      {
		icoo_[icoo_ld_*idx+0] = icoo_[i*icoo_ld_+0];
		icoo_[icoo_ld_*idx+1] = icoo_[i*icoo_ld_+1];
		icoo_[icoo_ld_*idx+2] = irot;
		++idx;
	      } 
	    N = N-(n1d-irot);
	  } }
      
      return WMESH_STATUS_SUCCESS;
    }

  //
  // Wedge.
  //
  if (cell_type_==2)
    {
      /*
	6 
	3 7
	1 4 8 
	0 2 5 9
      */
      wmesh_int_t idx = 0;
      for (wmesh_int_t k=0;k<n1d;++k)
	{
	  for (wmesh_int_t i=0;i<n1d;i++)
	    {
	      for (wmesh_int_t j=0;j<=i;j++)
		{
		  icoo_[icoo_ld_*idx+0] = j;
		  icoo_[icoo_ld_*idx+1] = i-j;
		  icoo_[icoo_ld_*idx+2] = k;
		  ++idx;
		} 
	    } 
	}
      return WMESH_STATUS_SUCCESS;
    }

  //
  // Pyramid.
  //
  if (cell_type_==1)
    {
      for (wmesh_int_t k=0;k<n1d;++k)
	{
	  for (wmesh_int_t i=0;i<n1d-k;++i)
	    {		
	      for (wmesh_int_t j=0;j<n1d-k;j++)
		{
		  const wmesh_int_t idx = ( k*(n1d)*(n1d) +i * (n1d) + j);
		  icoo_[icoo_ld_*idx+0] = i;
		  icoo_[icoo_ld_*idx+1] = j;
		  icoo_[icoo_ld_*idx+2] = k;
		}
	    } 
	}
      return WMESH_STATUS_SUCCESS;
    }

  
  if (cell_type_==3)
    {      
      for (wmesh_int_t k=0;k<n1d;++k)
	{
	  for (wmesh_int_t i=0;i<n1d;++i)
	    {		
	      for (wmesh_int_t j=0;j<n1d;j++)
		{
		  const wmesh_int_t idx = k * n1d * n1d + i * n1d + j;
		  icoo_[icoo_ld_*idx+0] = i;
		  icoo_[icoo_ld_*idx+1] = j;
		  icoo_[icoo_ld_*idx+2] = k;
		}
	    } 
	}
      return WMESH_STATUS_SUCCESS;
    }
  
  return WMESH_STATUS_INVALID_ARGUMENT;
}
	
	
static wmesh_status_t wmesh_treilli_interpolate_coo(wmesh_int_t cell_type_,
						    wmesh_int_t degree_,
						    wmesh_int_p icoo_,
						    wmesh_int_t icoo_ld_) 
{
  if (cell_type_==3)
    {
      wmesh_int_t idx = 0;
      
      icoo_[idx*icoo_ld_+0] = 0;	  
      icoo_[idx*icoo_ld_+1] = 0;	
      icoo_[idx*icoo_ld_+2] = 0;	
      ++idx;
      
      icoo_[idx*icoo_ld_+0] = degree_;	  
      icoo_[idx*icoo_ld_+1] = 0;	
      icoo_[idx*icoo_ld_+2] = 0;	
      ++idx;

      icoo_[idx*icoo_ld_+0] = degree_;	  
      icoo_[idx*icoo_ld_+1] = degree_;	
      icoo_[idx*icoo_ld_+2] = 0;	
      ++idx;
      
      icoo_[idx*icoo_ld_+0] = 0;	  
      icoo_[idx*icoo_ld_+1] = degree_;	
      icoo_[idx*icoo_ld_+2] = 0;	
      ++idx;
      

      icoo_[idx*icoo_ld_+0] = 0;	  
      icoo_[idx*icoo_ld_+1] = 0;	
      icoo_[idx*icoo_ld_+2] = degree_;	
      ++idx;

      icoo_[idx*icoo_ld_+0] = degree_;	  
      icoo_[idx*icoo_ld_+1] = 0;	
      icoo_[idx*icoo_ld_+2] = degree_;	
      ++idx;

      icoo_[idx*icoo_ld_+0] = degree_;	  
      icoo_[idx*icoo_ld_+1] = degree_;	
      icoo_[idx*icoo_ld_+2] = degree_;	
      ++idx;

      icoo_[idx*icoo_ld_+0] = 0;	  
      icoo_[idx*icoo_ld_+1] = degree_;	
      icoo_[idx*icoo_ld_+2] = degree_;	
      ++idx;

      static const wmesh_int_t hexaedge_cnc[]  = { 0,1,
						   1,2,
						   2,3,
						   3,0,
						   4,5,
						   5,6,
						   6,7,
						   7,4,
						   0,4,
						   1,5,
						   2,6,
						   3,7 }; 
      /* deja reference dans eVolume */
      static const wmesh_int_t hexaface_cnc[]  = {0,3,2,1,
						  4,5,6,7,				 
						  0,1,5,4,
						  1,2,6,5,				 
						  2,3,7,6,
						  3,0,4,7 };
      
      static const wmesh_int_t 	refhexa_icoo[] =  { 0,0,0,
						    1,0,0,
						    1,1,0,
						    0,1,0,
						    0,0,1,
						    1,0,1,
						    1,1,1,
						    0,1,1};
#if 0	
      static const wmesh_int_t 	refhexa_icoo[] =  { 1,1,0,
						    0,1,0,
						    0,0,0,
						    1,0,0,
						    1,1,1,
						    0,1,1,
						    0,0,1,
						    1,0,1};
#endif
      { wmesh_int_t reedge[8];
	for (wmesh_int_t iedge=0;iedge<12;++iedge)
	  {
	    reedge[0] = refhexa_icoo[3*hexaedge_cnc[iedge * 2 + 0]+0];
	    reedge[1] = refhexa_icoo[3*hexaedge_cnc[iedge * 2 + 0]+1];
	    reedge[2] = refhexa_icoo[3*hexaedge_cnc[iedge * 2 + 0]+2];

	    reedge[3] = refhexa_icoo[3*hexaedge_cnc[iedge * 2 + 1]+0];
	    reedge[4] = refhexa_icoo[3*hexaedge_cnc[iedge * 2 + 1]+1];
	    reedge[5] = refhexa_icoo[3*hexaedge_cnc[iedge * 2 + 1]+2];
		
	    for (wmesh_int_t i=0;i<degree_-1;++i)
	      {
		const wmesh_int_t l1 		= (i+1);
		const wmesh_int_t l0 		= degree_-(i+1);
		icoo_[idx*icoo_ld_+0] 		= reedge[0] * l0 + reedge[3] * l1;
		icoo_[idx*icoo_ld_+1] 		= reedge[1] * l0 + reedge[4] * l1;
		icoo_[idx*icoo_ld_+2] 		= reedge[2] * l0 + reedge[5] * l1;
		++idx;
	      } 
	  }
      }

      { wmesh_int_t reface[32];
	for (wmesh_int_t iface=0;iface<6;++iface)
	  {
	    reface[0] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 0]+0];
	    reface[1] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 0]+1];
	    reface[2] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 0]+2];

	    reface[3] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 1]+0];
	    reface[4] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 1]+1];
	    reface[5] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 1]+2];

	    reface[6] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 2]+0];
	    reface[7] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 2]+1];
	    reface[8] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 2]+2];

	    reface[9]  = refhexa_icoo[3*hexaface_cnc[iface * 4 + 3]+0];
	    reface[10] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 3]+1];
	    reface[11] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 3]+2];

	    for (wmesh_int_t i=0;i<degree_-1;++i)
	      {
		for (wmesh_int_t j=0;j<degree_-1;++j)
		  {
			
		    const wmesh_int_t l0 			= (i+1)*(j+1);
		    const wmesh_int_t l1 			= (degree_-(i+1))*(j+1);
		    const wmesh_int_t l2 			= (degree_-(i+1))*(degree_-(j+1));
		    const wmesh_int_t l3 			= (degree_-(j+1))*(i+1);
			
		    const wmesh_int_t xi 			= reface[0] * l0 + reface[3] * l1 + reface[6] * l2 + reface[9] *l3;
		    const wmesh_int_t xj 			= reface[1] * l0 + reface[4] * l1 + reface[7] * l2 + reface[10]*l3;
		    const wmesh_int_t xk 			= reface[2] * l0 + reface[5] * l1 + reface[8] * l2 + reface[11]*l3;
			
		    icoo_[idx*icoo_ld_+0] 	= xi/degree_;	  
		    icoo_[idx*icoo_ld_+1] 	= xj/degree_;	
		    icoo_[idx*icoo_ld_+2] 	= xk/degree_;
		    ++idx;
		  }
	      } 
	  }
      }
	

      /* wmesh_int_tNSwmesh_int_tDE */
      for (wmesh_int_t k=0;k<degree_-1;++k)
	{	    
	  for (wmesh_int_t i=0;i<degree_-1;++i)
	    {
	      for (wmesh_int_t j=0;j<degree_-1;j++)
		{			  
		  icoo_[idx*icoo_ld_+0] = i+1;
		  icoo_[idx*icoo_ld_+1] = j+1;	
		  icoo_[idx*icoo_ld_+2] = k+1;
		  ++idx;
		}
	    }
	}
      return WMESH_STATUS_SUCCESS;
    }  

  if (cell_type_ == 0)	
    {
      wmesh_int_t idx = 0;
      static constexpr const wmesh_int_t reftetra_icoo[] = {0,0,0,
							    1,0,0,
							    0,1,0,
							    0,0,1};
      
      const wmesh_int_t 	n 		= (degree_>0)? degree_-1 : 0;
      icoo_[idx*icoo_ld_+0] = 0;	  
      icoo_[idx*icoo_ld_+1] = 0;	
      icoo_[idx*icoo_ld_+2] = 0;	
      ++idx;
      icoo_[idx*icoo_ld_+0] = degree_;	  
      icoo_[idx*icoo_ld_+1] = 0;	
      icoo_[idx*icoo_ld_+2] = 0;	
      ++idx;
      icoo_[idx*icoo_ld_+0] = 0;	  
      icoo_[idx*icoo_ld_+1] = degree_;	
      icoo_[idx*icoo_ld_+2] = 0;	
      ++idx;
      icoo_[idx*icoo_ld_+0] = 0;	  
      icoo_[idx*icoo_ld_+1] = 0;	
      icoo_[idx*icoo_ld_+2] = degree_;	
      ++idx;
      
      static const wmesh_int_t tetraedge_cnc[] = { 1,2,
						   2,0,
						   0,1,
						   2,3,
						   0,3,
						   1,3};			
      
      static const wmesh_int_t tetraface_cnc[] = {1,2,3,
						  2,0,3,
						  0,1,3,
						  0,2,1};
      
      { wmesh_int_t reedge[8];
	for (wmesh_int_t iedge=0;iedge<6;++iedge)
	  {
	    reedge[0] = reftetra_icoo[3*tetraedge_cnc[iedge * 2 + 0]+0];
	    reedge[1] = reftetra_icoo[3*tetraedge_cnc[iedge * 2 + 0]+1];
	    reedge[2] = reftetra_icoo[3*tetraedge_cnc[iedge * 2 + 0]+2];
	    reedge[3] = reftetra_icoo[3*tetraedge_cnc[iedge * 2 + 1]+0];
	    reedge[4] = reftetra_icoo[3*tetraedge_cnc[iedge * 2 + 1]+1];
	    reedge[5] = reftetra_icoo[3*tetraedge_cnc[iedge * 2 + 1]+2];
	    for (wmesh_int_t i=0;i<n;++i)
	      {
		const wmesh_int_t l1 = (i+1);
		const wmesh_int_t l0 = degree_-(i+1);

		icoo_[idx*icoo_ld_+0] = reedge[0] * l0 + reedge[3] * l1;
		icoo_[idx*icoo_ld_+1] = reedge[1] * l0 + reedge[4] * l1;
		icoo_[idx*icoo_ld_+2] = reedge[2] * l0 + reedge[5] * l1;			    
		++idx;
	      } 
	  } } 
      
      
      if (degree_>2)
	{
	  { wmesh_int_t reface[32];
	    for (wmesh_int_t iface=0;iface<4;++iface)
	      {
		reface[0] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 0]+0];
		reface[1] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 0]+1];
		reface[2] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 0]+2];

		reface[3] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 1]+0];
		reface[4] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 1]+1];
		reface[5] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 1]+2];

		reface[6] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 2]+0];
		reface[7] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 2]+1];
		reface[8] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 2]+2];
		  
		for (wmesh_int_t i=0;i<n;++i)
		  {
		    for (wmesh_int_t j=0;j<n-(i+1);++j)
		      {
			const wmesh_int_t l0 		= degree_-( (i+1) + (j+1) );
			const wmesh_int_t l1 		= (i+1);
			const wmesh_int_t l2 		= (j+1);

			  
			  
			icoo_[idx*icoo_ld_+0] = reface[0] * l0 + reface[3] * l1 + reface[6] * l2;
			icoo_[idx*icoo_ld_+1] = reface[1] * l0 + reface[4] * l1 + reface[7] * l2;
			icoo_[idx*icoo_ld_+2] = reface[2] * l0 + reface[5] * l1 + reface[8] * l2;

			++idx;
		      } 
		  }
	      		  
	      } } 
	}
      /* INSIDE TETRA */

      for (wmesh_int_t  k=0;k<n-1;++k)
	{
	  for (wmesh_int_t  i=0;i<n-(k+1);++i)
	    {		    
	      for (wmesh_int_t j=0;j<n-(i+1)-(k+1);j++)
		{			  
		  icoo_[idx*icoo_ld_+0] = i+1;
		  icoo_[idx*icoo_ld_+1] = j+1;	
		  icoo_[idx*icoo_ld_+2] = k+1;	
		  ++idx;
		} 
	    } 
	}
      return WMESH_STATUS_SUCCESS;
    }
  
  if (cell_type_==2)
    {
      const wmesh_int_t	n = (degree_>0)? degree_-1 : 0;
      wmesh_int_t idx = 0;
      /* SOMMETS */
      icoo_[icoo_ld_*idx+0] = 0;
      icoo_[icoo_ld_*idx+1] = 0;
      icoo_[icoo_ld_*idx+2] = 0;
      ++idx;
      icoo_[icoo_ld_*idx+0] = degree_;
      icoo_[icoo_ld_*idx+1] = 0;
      icoo_[icoo_ld_*idx+2] = 0;
      ++idx;
      icoo_[icoo_ld_*idx+0] = 0;
      icoo_[icoo_ld_*idx+1] = degree_;
      icoo_[icoo_ld_*idx+2] = 0;
      ++idx;

      icoo_[icoo_ld_*idx+0] = 0;
      icoo_[icoo_ld_*idx+1] = 0;
      icoo_[icoo_ld_*idx+2] = degree_;
      ++idx;
      icoo_[icoo_ld_*idx+0] = degree_;
      icoo_[icoo_ld_*idx+1] = 0;
      icoo_[icoo_ld_*idx+2] = degree_;
      ++idx;
      icoo_[icoo_ld_*idx+0] = 0;
      icoo_[icoo_ld_*idx+1] = degree_;
      icoo_[icoo_ld_*idx+2] = degree_;
      ++idx;

      /* EDGE 0 */
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = (i+1);	  
	  icoo_[icoo_ld_*idx+1] = 0;	
	  icoo_[icoo_ld_*idx+2] = 0;	
	  ++idx;	      
	} 
      /* EDGE 1 */
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = degree_-(i+1);	  
	  icoo_[icoo_ld_*idx+1] = (i+1);	
	  icoo_[icoo_ld_*idx+2] = 0;	
	  ++idx;	      
	} 
      /* EDGE 2 */
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = 0;	  
	  icoo_[icoo_ld_*idx+1] = degree_-(i+1);
	  icoo_[icoo_ld_*idx+2] = 0;	
	  ++idx;	      
	} 


      /* EDGE 3 */
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = (i+1);	  
	  icoo_[icoo_ld_*idx+1] = 0;	
	  icoo_[icoo_ld_*idx+2] = degree_;	
	  ++idx;	      
	} 
      /* EDGE 4 */
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = degree_-(i+1);	  
	  icoo_[icoo_ld_*idx+1] = (i+1);	
	  icoo_[icoo_ld_*idx+2] = degree_;	
	  ++idx;	      
	} 
      /* EDGE 5 */
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = 0;	  
	  icoo_[icoo_ld_*idx+1] = degree_-(i+1);
	  icoo_[icoo_ld_*idx+2] = degree_;	
	  ++idx;	      
	} 

      /* EDGE 6 */
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = 0;	  
	  icoo_[icoo_ld_*idx+1] = 0;	
	  icoo_[icoo_ld_*idx+2] = (i+1);	
	  ++idx;	      
	} 
      /* EDGE 7 */
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = degree_;	  
	  icoo_[icoo_ld_*idx+1] = 0;	
	  icoo_[icoo_ld_*idx+2] = (i+1);	
	  ++idx;	      
	} 
      /* EDGE 8 */
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = 0;	  
	  icoo_[icoo_ld_*idx+1] = degree_;
	  icoo_[icoo_ld_*idx+2] = (i+1);	
	  ++idx;	      
	} 

      static const wmesh_int_t wedgequadface_cnc[] = {0,2,5,3,
						      4,5,2,1,
						      0,3,4,1};
	  
      static const wmesh_int_t 	refwedge_icoo[] =  { 0,0,0,
						     1,0,0,
						     0,1,0,
						     0,0,1,
						     1,0,1,
						     0,1,1};

      for (wmesh_int_t iface=0;iface<2;++iface)
	{
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      for (wmesh_int_t j=0;j<n-(i+1);++j)
		{
		  icoo_[icoo_ld_*idx+0] = (i+1);	  
		  icoo_[icoo_ld_*idx+1] = (j+1);
		  icoo_[icoo_ld_*idx+2] = iface*degree_;
		  ++idx;	      
		} 
	    }       
	} 

      { wmesh_int_t reface[32];
	for (wmesh_int_t iface=0;iface<3;++iface)
	  {
	    reface[0] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 0]+0];
	    reface[1] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 0]+1];
	    reface[2] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 0]+2];
		  
	    reface[3] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 1]+0];
	    reface[4] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 1]+1];
	    reface[5] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 1]+2];
		  
	    reface[6] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 2]+0];
	    reface[7] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 2]+1];
	    reface[8] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 2]+2];
		
	    reface[9]  = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 3]+0];
	    reface[10] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 3]+1];
	    reface[11] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 3]+2];


	    for (wmesh_int_t i=0;i<n;++i)
	      {
		for (wmesh_int_t j=0;j<n;++j)
		  {
		    const wmesh_int_t l0 			= (i+1)*(j+1);
		    const wmesh_int_t l1 			= (degree_-(i+1))*(j+1);
		    const wmesh_int_t l2 			= (degree_-(i+1))*(degree_-(j+1));
		    const wmesh_int_t l3 			= (degree_-(j+1))*(i+1);
			    
		    const wmesh_int_t xi 			= reface[0] * l0 + reface[3] * l1 + reface[6] * l2 + reface[9] *l3;
		    const wmesh_int_t xj 			= reface[1] * l0 + reface[4] * l1 + reface[7] * l2 + reface[10]*l3;
		    const wmesh_int_t xk 			= reface[2] * l0 + reface[5] * l1 + reface[8] * l2 + reface[11]*l3;
			    
		    icoo_[icoo_ld_*idx+0] 	= xi/degree_;	  
		    icoo_[icoo_ld_*idx+1] 	= xj/degree_;	
		    icoo_[icoo_ld_*idx+2] 	= xk/degree_;	
		    ++idx;
		  } 
	      } 

	  } }


      /* interior */
      for (wmesh_int_t k=0;k<n;++k)
	{

	  for (wmesh_int_t i=0;i<n;++i)
	    {

	      for (wmesh_int_t j=0;j<n-(i+1);++j)
		{
		  icoo_[icoo_ld_*idx+0] = (i+1);	  
		  icoo_[icoo_ld_*idx+1] = (j+1);
		  icoo_[icoo_ld_*idx+2] = (k+1);
		  ++idx;	      
		} 
	    } 
	} 
      return WMESH_STATUS_SUCCESS;

    }



  
  if (cell_type_==1)
    {
      const wmesh_int_t	n = (degree_>0)? degree_-1 : 0;
      wmesh_int_t idx = 0;
      /* SOMMETS */
      icoo_[icoo_ld_*idx+0] = 0;
      icoo_[icoo_ld_*idx+1] = 0;
      icoo_[icoo_ld_*idx+2] = 0;
      ++idx;
      icoo_[icoo_ld_*idx+0] = degree_;
      icoo_[icoo_ld_*idx+1] = 0;
      icoo_[icoo_ld_*idx+2] = 0;
      ++idx;
      icoo_[icoo_ld_*idx+0] = degree_;
      icoo_[icoo_ld_*idx+1] = degree_;
      icoo_[icoo_ld_*idx+2] = 0;
      ++idx;
      icoo_[icoo_ld_*idx+0] = 0;
      icoo_[icoo_ld_*idx+1] = degree_;
      icoo_[icoo_ld_*idx+2] = 0;
      ++idx;
      icoo_[icoo_ld_*idx+0] = 0;
      icoo_[icoo_ld_*idx+1] = 0;
      icoo_[icoo_ld_*idx+2] = degree_;
      ++idx;
      
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = (i+1);	  
	  icoo_[icoo_ld_*idx+1] = 0;	
	  icoo_[icoo_ld_*idx+2] = 0;	
	  ++idx;	      
	}
      
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = degree_; 
	  icoo_[icoo_ld_*idx+1] = (i+1);	
	  icoo_[icoo_ld_*idx+2] = 0;	
	  ++idx;	      
	} 
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = degree_-(i+1); 
	  icoo_[icoo_ld_*idx+1] = degree_;	
	  icoo_[icoo_ld_*idx+2] = 0;	
	  ++idx;	      
	}
      
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = 0; 
	  icoo_[icoo_ld_*idx+1] = degree_-(i+1);	
	  icoo_[icoo_ld_*idx+2] = 0;	
	  ++idx;	      
	}


      
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = 0; 
	  icoo_[icoo_ld_*idx+1] = 0;	
	  icoo_[icoo_ld_*idx+2] = (i+1);	
	  ++idx;	      
	}
      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = degree_-(i+1); 
	  icoo_[icoo_ld_*idx+1] = 0;	
	  icoo_[icoo_ld_*idx+2] = (i+1);	
	  ++idx;	      
	}

      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = degree_-(i+1);
	  icoo_[icoo_ld_*idx+1] = degree_-(i+1);	
	  icoo_[icoo_ld_*idx+2] = (i+1);	
	  ++idx;	      
	}

      for (wmesh_int_t i=0;i<n;++i)
	{
	  icoo_[icoo_ld_*idx+0] = 0;
	  icoo_[icoo_ld_*idx+1] = degree_-(i+1);	
	  icoo_[icoo_ld_*idx+2] = (i+1);	
	  ++idx;	      
	}


      

      
      static const wmesh_int_t s_q2n[] = {0,3,2,1};
      static const wmesh_int_t s_t2n[] = {0,1,4,
					  1,2,4, 
					  2,3,4, 
					  3,0,4};
      
      static const wmesh_int_t 	refpyr_icoo[] =  { 0,0,0,
						   1,0,0,
						   1,1,0,
						   0,1,0,
						   0,0,1};
      //      std::cout << "idx " << idx << std::endl;           
      { wmesh_int_t reface[32];
	for (wmesh_int_t iface=0;iface<4;++iface)
	  {
	    
	    reface[0] = refpyr_icoo[3*s_t2n[iface * 3 + 0]+0];
	    reface[1] = refpyr_icoo[3*s_t2n[iface * 3 + 0]+1];
	    reface[2] = refpyr_icoo[3*s_t2n[iface * 3 + 0]+2];
		  
	    reface[3] = refpyr_icoo[3*s_t2n[iface * 3 + 1]+0];
	    reface[4] = refpyr_icoo[3*s_t2n[iface * 3 + 1]+1];
	    reface[5] = refpyr_icoo[3*s_t2n[iface * 3 + 1]+2];
		  
	    reface[6] = refpyr_icoo[3*s_t2n[iface * 3 + 2]+0];
	    reface[7] = refpyr_icoo[3*s_t2n[iface * 3 + 2]+1];
	    reface[8] = refpyr_icoo[3*s_t2n[iface * 3 + 2]+2];

	    for (wmesh_int_t i=0;i<n;++i)
	      {
		for (wmesh_int_t j=0;j<n-i-1;++j)
		  {
		    const wmesh_int_t l0 = degree_ - (i+1) - (j+1);
		    const wmesh_int_t l1 = (i+1);
		    const wmesh_int_t l2 = (j+1);
			    
		    const wmesh_int_t xi = reface[0] * l0 + reface[3] * l1 + reface[6] * l2;
		    const wmesh_int_t xj = reface[1] * l0 + reface[4] * l1 + reface[7] * l2;
		    const wmesh_int_t xk = reface[2] * l0 + reface[5] * l1 + reface[8] * l2;
			    
		    icoo_[icoo_ld_*idx+0] = xi;	  
		    icoo_[icoo_ld_*idx+1] = xj;	
		    icoo_[icoo_ld_*idx+2] = xk;
#if 0
		    std::cout << "----" << std::endl;
		    std::cout << "s " << xi << std::endl;
		    std::cout << "s " << xj << std::endl;
		    std::cout << "s " << xk << std::endl;
#endif
		    ++idx;
		  } 
	      } 
	  } }

      { wmesh_int_t reface[32];
	for (wmesh_int_t iface=0;iface<1;++iface)
	  {
	    
	    reface[0] = refpyr_icoo[3*s_q2n[iface * 4 + 0]+0];
	    reface[1] = refpyr_icoo[3*s_q2n[iface * 4 + 0]+1];
	    reface[2] = refpyr_icoo[3*s_q2n[iface * 4 + 0]+2];
		  
	    reface[3] = refpyr_icoo[3*s_q2n[iface * 4 + 1]+0];
	    reface[4] = refpyr_icoo[3*s_q2n[iface * 4 + 1]+1];
	    reface[5] = refpyr_icoo[3*s_q2n[iface * 4 + 1]+2];
		  
	    reface[6] = refpyr_icoo[3*s_q2n[iface * 4 + 2]+0];
	    reface[7] = refpyr_icoo[3*s_q2n[iface * 4 + 2]+1];
	    reface[8] = refpyr_icoo[3*s_q2n[iface * 4 + 2]+2];
		
	    reface[9]  = refpyr_icoo[3*s_q2n[iface * 4 + 3]+0];
	    reface[10] = refpyr_icoo[3*s_q2n[iface * 4 + 3]+1];
	    reface[11] = refpyr_icoo[3*s_q2n[iface * 4 + 3]+2];

	    for (wmesh_int_t i=0;i<n;++i)
	      {
		for (wmesh_int_t j=0;j<n;++j)
		  {
		    const wmesh_int_t l0 = (i+1)*(j+1);
		    const wmesh_int_t l1 = (degree_-(i+1))*(j+1);
		    const wmesh_int_t l2 = (degree_-(i+1))*(degree_-(j+1));
		    const wmesh_int_t l3 = (degree_-(j+1))*(i+1);
#if 0
		    std::cout << l0 << std::endl;
		    std::cout << l1 << std::endl;
		    std::cout << l2 << std::endl;
		    std::cout << l3 << std::endl;
#endif
		    const wmesh_int_t xi = reface[0] * l0 + reface[3] * l1 + reface[6] * l2 + reface[9] *l3;
		    const wmesh_int_t xj = reface[1] * l0 + reface[4] * l1 + reface[7] * l2 + reface[10]*l3;
		    const wmesh_int_t xk = reface[2] * l0 + reface[5] * l1 + reface[8] * l2 + reface[11]*l3;
		    icoo_[icoo_ld_*idx+0] = xi / degree_;	  
		    icoo_[icoo_ld_*idx+1] = xj / degree_;	
		    icoo_[icoo_ld_*idx+2] = xk / degree_;

#if 0
		    std::cout << "idx " << idx << std::endl;
		    std::cout << icoo_[icoo_ld_*idx+0] << std::endl;
		    std::cout << icoo_[icoo_ld_*idx+1] << std::endl;
		    std::cout << icoo_[icoo_ld_*idx+2] << std::endl;
#endif

		    
		    ++idx;
		  } 
	      } 

	  } }
      
      /* 
	 interior 
	 degree 2 0
	 degree 3 1 
	 degree 4 1+4
	 degree 5 1+4+9
      */
      if (degree_>2)
	{
	  for (wmesh_int_t k=1;k<degree_-1;++k)
	    {
	      for (wmesh_int_t i=1;i<=degree_-1-k;++i)
		{
		  for (wmesh_int_t j=1;j<=degree_-1-k;++j)
		    {
		      icoo_[icoo_ld_*idx+0] = i;	  
		      icoo_[icoo_ld_*idx+1] = j;
		      icoo_[icoo_ld_*idx+2] = k;
		      // printf("rita %d %d %d\n",i,j,k);
		      ++idx;	      
		    } 
		} 
	    }
	}
      
      return WMESH_STATUS_SUCCESS;
    }

  return WMESH_STATUS_INVALID_ARGUMENT;

};	  


extern "C"
{
  
  wmesh_status_t wmesh_treilli_buffer_size(wmesh_int_t 	cell_type_,
					   wmesh_int_t 	degree_,
					   wmesh_int_p	work_n_,
					   wmesh_int_p	num_entities_)
  {
    WMESH_POINTER_CHECK(work_n_);
   
    work_n_[0] = 0;
    for (int i=0;i<WMESH_ELEMENT_NODE;++i)
      {
	num_entities_[i] = 0; 
      }
    
    if (cell_type_==0)
      {
	num_entities_[WMESH_ELEMENT_NODE] =  ( (degree_+1)*(degree_+2)*(degree_+3))/6;	
	num_entities_[WMESH_ELEMENT_TETRAHEDRON] = degree_*degree_*degree_;
      }
    else if (cell_type_==1)
      {
	num_entities_[WMESH_ELEMENT_NODE]  = ((degree_+1)*(degree_+2)*(2*degree_+3))/6;
	num_entities_[WMESH_ELEMENT_PYRAMID] = ((degree_-1)*(degree_)*(2*degree_-1))/3 + degree_ * degree_;
	num_entities_[WMESH_ELEMENT_TETRAHEDRON] = (2*(degree_-1)*(degree_)*(degree_+1))/3;
      }
    else if (cell_type_==2)
      {
	num_entities_[WMESH_ELEMENT_NODE] = ( ( (degree_+1) *(degree_+2) ) /2 ) * (degree_+1);
	num_entities_[WMESH_ELEMENT_WEDGE] = degree_*degree_*degree_;
      }
    else if (cell_type_==3)
      {
	num_entities_[WMESH_ELEMENT_NODE] = ( (degree_+1)*(degree_+1)*(degree_+1));
	num_entities_[WMESH_ELEMENT_HEXAHEDRON] = degree_*degree_*degree_;
      }
    else
      {
	return WMESH_STATUS_INVALID_ARGUMENT;
      }
    
    work_n_[0] = 4*(degree_+1)*(degree_+1)*(degree_+1) + 3 * (degree_+1)*(degree_+1) + num_entities_[WMESH_ELEMENT_NODE] * 3;
    return WMESH_STATUS_SUCCESS;
  };
  

  wmesh_status_t wmesh_treilli_calculate(    wmesh_int_t 		cell_type_,
					     wmesh_int_t 		degree_,
					     wmesh_int_t 		num_nodes_,
						 
					     const_wmesh_int_p 		c2n_ptr_,
					     const_wmesh_int_p 		c2n_m_,
					     const_wmesh_int_p 		c2n_n_,
					     wmesh_int_p 		c2n_v_,
					     const_wmesh_int_p 		c2n_ld_,
						 
					     wmesh_int_t 		icoo_m_,
					     wmesh_int_t 		icoo_n_,
					     wmesh_int_p 		icoo_v_,
					     wmesh_int_t 		icoo_ld_,
						 
					     wmesh_int_t		work_n_,
					     wmesh_int_p		work_)
  {

    wmesh_status_t status;    
#if 0
    {
      wmesh_int_t expected_work_n;
      status = wmesh_treilli_buffer_size(cell_type_,
					 degree_,
					 &expected_work_n);    
      if (work_n_ < expected_work_n)
	{
	  return WMESH_STATUS_ERROR_WORKSPACE;
	}
    }
#endif    
    wmesh_int_t perm_n 	= (degree_+1)*(degree_+1)*(degree_+1);
    wmesh_int_t ilagr_n = (degree_+1)*(degree_+1)*(degree_+1) * 3; // num_nodes*3;
    wmesh_int_t work_n 	= 3*(degree_+1)*(degree_+1);
    wmesh_int_p ilagr 	= work_;
    wmesh_int_p perm 	= work_ + (ilagr_n);
    wmesh_int_p work 	= work_ + (ilagr_n + perm_n);

    {
      status = wmesh_treilli_interpolate_coo(cell_type_,
					     degree_,
					     icoo_v_,
					     icoo_ld_);
      //
      // Compute the permutation.
      //
      for (wmesh_int_t i=0;i<num_nodes_;++i)
	{
	  perm[(degree_+1)*(degree_+1)*icoo_v_[icoo_ld_*i+2] + icoo_v_[icoo_ld_*i+0]*(degree_+1) + icoo_v_[icoo_ld_*i+1]] = 1+i;
	} 

      WMESH_STATUS_CHECK(status);
    }

    //
    // Compute the connectivity.
    //
    {
      status = wmesh_treilli_calculate_c2n(cell_type_,
					   degree_,
					   c2n_ptr_,
					   c2n_m_,
					   c2n_n_,
					   c2n_v_,
					   c2n_ld_,
					   work_n,
					   work);
      WMESH_STATUS_CHECK(status);
    }

    {
      //
      // Compute the integers.
      //
      status = wmesh_treilli_calculate_icoo(cell_type_,
					    degree_,
					    ilagr,
					    3);
      WMESH_STATUS_CHECK(status);
    }
    //
    // Adjust the connectivity.
    //

    {
      for (wmesh_int_t subcell_type=0;subcell_type<4;++subcell_type)
	{
	  wmesh_int_t nsubcells = c2n_n_[subcell_type];
	  if (nsubcells > 0)
	    {
	      for (wmesh_int_t idx = 0;idx < nsubcells;++idx)
		{
		  for (wmesh_int_t iv=0;iv<c2n_m_[subcell_type];++iv)
		    {	  
		      const wmesh_int_t l	= c2n_v_[ c2n_ptr_[subcell_type] + c2n_ld_[subcell_type] * idx + iv];
			
		      const wmesh_int_t i	= ilagr[3*l+0];
		      const wmesh_int_t j	= ilagr[3*l+1];
		      const wmesh_int_t k	= ilagr[3*l+2];
			
		      c2n_v_[c2n_ptr_[subcell_type] + c2n_ld_[subcell_type] * idx + iv] = perm[(degree_+1)*(degree_+1)*k+(degree_+1)*i+j];
		    }
		}
	    }
	}	
    }
    
    return WMESH_STATUS_SUCCESS;
  };



  

  wmesh_status_t wmesh_treilli(wmesh_t ** 	mesh__,
			       wmesh_int_t 	cell_type_,
			       wmesh_int_t 	degree_,
			       wmesh_int_t	work_n_,
			       wmesh_int_p	work_)
  {
    wmesh_status_t status;
    wmesh_int_t	num_entities[WMESH_ELEMENT_ALL];
    for (int i=0;i<WMESH_ELEMENT_ALL;++i)
      num_entities[i]=0;
    {
      wmesh_int_t expected_work_n;
      status = wmesh_treilli_buffer_size(cell_type_,
					 degree_,
					 &expected_work_n,
					 num_entities);
      WMESH_STATUS_CHECK(status);
      if (work_n_ < expected_work_n)
	{
	  return WMESH_STATUS_ERROR_WORKSPACE;
	}
    }
    
    wmesh_int_t
      num_nodes,
      num_cells[4];
    
    num_nodes = num_entities[WMESH_ELEMENT_NODE];
    num_cells[0] = num_entities[WMESH_ELEMENT_TETRAHEDRON];
    num_cells[1] = num_entities[WMESH_ELEMENT_PYRAMID];
    num_cells[2] = num_entities[WMESH_ELEMENT_WEDGE];
    num_cells[3] = num_entities[WMESH_ELEMENT_HEXAHEDRON];

    wmesh_int_t icoo_n 	= num_nodes*3;
    wmesh_int_t work_n 	= work_n_ - icoo_n;

    WMESH_POINTER_CHECK(mesh__);

    wmesh_int_t icoo_ld	= 3;
    wmesh_int_p icoo 	= work_;
    wmesh_int_p work 	= work_ + icoo_n;

    wmesh_int_t c2n_ptr[5]{};
    wmesh_int_t c2n_m[4];
    wmesh_int_t c2n_n[4];
    wmesh_int_t c2n_ld[4];
    wmesh_int_p c2n_v;

    c2n_m[0] = 4;
    c2n_m[1] = 5;
    c2n_m[2] = 6;
    c2n_m[3] = 8;
    
    c2n_n[0] = num_cells[0];
    c2n_n[1] = num_cells[1];
    c2n_n[2] = num_cells[2];
    c2n_n[3] = num_cells[3];

    c2n_ld[0] = c2n_m[0];
    c2n_ld[1] = c2n_m[1];
    c2n_ld[2] = c2n_m[2];
    c2n_ld[3] = c2n_m[3];
    
    c2n_ptr[0] = 0;
    c2n_ptr[1] = c2n_ptr[0]+num_cells[0] * c2n_ld[0];
    c2n_ptr[2] = c2n_ptr[1]+num_cells[1] * c2n_ld[1];
    c2n_ptr[3] = c2n_ptr[2]+num_cells[2] * c2n_ld[2];
    c2n_ptr[4] = c2n_ptr[3]+num_cells[3] * c2n_ld[3];
    
    c2n_v = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*c2n_ptr[4]);


    wmesh_treilli_calculate(cell_type_,
			    degree_,
			    num_nodes,
			    
			    c2n_ptr,
			    c2n_m,
			    c2n_n,
			    c2n_v,
			    c2n_ld,

			    3,
			    num_nodes,
			    icoo,
			    icoo_ld,
			    
			    work_n,
			    work);
    
    double * coo = (double*)malloc(sizeof(double)*num_nodes*3);
    double idegree = ((double)1.0)/((double)degree_);
    for (wmesh_int_t i=0;i<num_nodes*3;++i)
      {
	coo[i] = icoo[i]*idegree;
      }
    
    status = wmesh_def(mesh__,
		       num_nodes,
		       4,
		       c2n_ptr,
		       c2n_m,
		       c2n_n,
		       c2n_v,
		       c2n_ld,
		       coo,
		       3);
    wmesh_write_medit(mesh__[0],"debug.mesh");    
    return WMESH_STATUS_SUCCESS;
  };


#if 0
    wmesh_status_t wmesh_treilli(wmesh_t ** 	mesh__,
			       wmesh_int_t 	cell_type_,
			       wmesh_int_t 	degree_,
			       wmesh_int_t	work_n_,
			       wmesh_int_p	work_)
  {
    wmesh_status_t status;
    wmesh_int_t	num_entities[WMESH_ELEMENT_ALL];
    
    {
      wmesh_int_t expected_work_n;
      status = wmesh_treilli_buffer_size(cell_type_,
					 degree_,
					 &expected_work_n,
					 num_entities);    
      if (work_n_ < expected_work_n)
	{
	  return WMESH_STATUS_ERROR_WORKSPACE;
	}
    }
    
    wmesh_int_t
      num_nodes,
      num_cells[4];
    
    num_nodes = num_entities[WMESH_ELEMENT_NODE];
    num_cells[0] = num_entities[WMESH_ELEMENT_TETRAHEDRON];
    num_cells[1] = num_entities[WMESH_ELEMENT_PYRAMID];
    num_cells[2] = num_entities[WMESH_ELEMENT_WEDGE];
    num_cells[3] = num_entities[WMESH_ELEMENT_HEXAHEDRON];

#if 0
    if (cell_type_==0)
      {
	num_nodes = ( (degree_+1)*(degree_+2)*(degree_+3))/6;
	num_cells[cell_type_] = degree_*degree_*degree_;
      }
    else if (cell_type_==1)
      {
	num_nodes = ((degree_+1)*(degree_+2)*(2*degree_+3))/6;

	// num pyr
	num_cells[cell_type_] = ((degree_-1)*(degree_)*(2*degree_-1))/3 + degree_ * degree_;

	// num tet
	num_cells[0] = (2*(degree_-1)*(degree_)*(degree_+1))/3;
	
	// K  NPYR NTET
	// 1   1                        0
	// 2   1 + 4+1                  4
	// 3   1 + 4+1 + 9+4            4 + 12 
	// 4   1 + 1+4 + 4+9 + 9+16     4 + 12 + sum (2*(k-1)*k)
	// 1 6 19 44
	// 4
	// npyramid 2*(((k-1)(k)(2*k-1)/6) + k*k
	// ntet = sum (2*(k-1)*k) = 2*(sum(k*k) - sum(k)) =(2*(k-1)*k*(k+1))/3
      }    
    else if (cell_type_==2)
      {
	num_nodes = ( ( (degree_+1) *(degree_+2) ) / 2 ) * (degree_+1);
	num_cells[cell_type_] = degree_*degree_*degree_;
      }
    else  // if (cell_type_==3)
      {
	num_nodes = ( (degree_+1)*(degree_+1)*(degree_+1));
	num_cells[cell_type_] = degree_*degree_*degree_;
      }
#endif
    
    wmesh_int_t perm_n 	= (degree_+1)*(degree_+1)*(degree_+1);
    wmesh_int_t ilagr_n = (degree_+1)*(degree_+1)*(degree_+1) * 3;
    wmesh_int_t icoo_n 	= num_nodes*3;
    wmesh_int_t work_n 	= 3*(degree_+1)*(degree_+1);

    WMESH_POINTER_CHECK(mesh__);

    wmesh_int_t icoo_ld	= 3;
    wmesh_int_p ilagr 	= work_;
    wmesh_int_p icoo 	= work_ + (ilagr_n);
    wmesh_int_p perm 	= work_ + (ilagr_n + icoo_n);
    wmesh_int_p work 	= work_ + (ilagr_n + icoo_n + perm_n);

    wmesh_int_t c2n_ptr[5]{};
    wmesh_int_t c2n_m[4];
    wmesh_int_t c2n_n[4];
    wmesh_int_t c2n_ld[4];
    wmesh_int_p c2n_v;

    c2n_m[0] = 4;
    c2n_m[1] = 5;
    c2n_m[2] = 6;
    c2n_m[3] = 8;
    
    c2n_n[0] = num_cells[0];
    c2n_n[1] = num_cells[1];
    c2n_n[2] = num_cells[2];
    c2n_n[3] = num_cells[3];

    c2n_ld[0] = c2n_m[0];
    c2n_ld[1] = c2n_m[1];
    c2n_ld[2] = c2n_m[2];
    c2n_ld[3] = c2n_m[3];
    
    c2n_ptr[0] = 0;
    c2n_ptr[1] = c2n_ptr[0]+num_cells[0] * c2n_ld[0];
    c2n_ptr[2] = c2n_ptr[1]+num_cells[1] * c2n_ld[1];
    c2n_ptr[3] = c2n_ptr[2]+num_cells[2] * c2n_ld[2];
    c2n_ptr[4] = c2n_ptr[3]+num_cells[3] * c2n_ld[3];
    
    c2n_v = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*c2n_ptr[4]);

    //
    // Compute the scaled integer coordinates.
    //

    {
      status = wmesh_treilli_interpolate_coo(cell_type_,
					     degree_,
					     icoo,
					     icoo_ld);
      //
      // Compute the permutation.
      //
      for (wmesh_int_t i=0;i<num_nodes;++i)
	{
	  perm[(degree_+1)*(degree_+1)*icoo[3*i+2] + icoo[3*i+0]*(degree_+1) + icoo[3*i+1]] = 1+i;
	} 

      WMESH_STATUS_CHECK(status);
    }

    //
    // Compute the connectivity.
    //
    {
      status = wmesh_treilli_calculate_c2n(cell_type_,
					   degree_,
					   c2n_ptr,
					   c2n_m,
					   c2n_n,
					   c2n_v,
					   c2n_ld,
					   work_n,
					   work);
      WMESH_STATUS_CHECK(status);
    }

    {
      //
      // Compute the integers.
      //
      status = wmesh_treilli_calculate_icoo(cell_type_,
					    degree_,
					    ilagr,
					    3);
      WMESH_STATUS_CHECK(status);
    }
    //
    // Adjust the connectivity.
    //

    {
      for (wmesh_int_t subcell_type=0;subcell_type<4;++subcell_type)
	{
	  wmesh_int_t nsubcells = c2n_n[subcell_type];
	  if (nsubcells > 0)
	    {
	      for (wmesh_int_t idx = 0;idx < nsubcells;++idx)
		{
		  for (wmesh_int_t iv=0;iv<c2n_m[subcell_type];++iv)
		    {	  
		      const wmesh_int_t l	= c2n_v[ c2n_ptr[subcell_type] + c2n_ld[subcell_type] * idx + iv];
			
		      const wmesh_int_t i	= ilagr[3*l+0];
		      const wmesh_int_t j	= ilagr[3*l+1];
		      const wmesh_int_t k	= ilagr[3*l+2];
			
		      c2n_v[c2n_ptr[subcell_type] + c2n_ld[subcell_type] * idx + iv] = perm[(degree_+1)*(degree_+1)*k+(degree_+1)*i+j];
		    }
		}
	    }
	}	
    }
    
    double * coo = (double*)malloc(sizeof(double)*num_nodes*3);
    double idegree = ((double)1.0)/((double)degree_);
    for (wmesh_int_t i=0;i<num_nodes;++i)
      {
	coo[3*i+0] = icoo[3*i+0]*idegree;
	coo[3*i+1] = icoo[3*i+1]*idegree;
	coo[3*i+2] = icoo[3*i+2]*idegree;
      }

#if 0
    std::cout << " " << c2n_n[0] << std::endl;
    std::cout << " " << c2n_n[1] << std::endl;
    std::cout << " " << c2n_n[2] << std::endl;
    std::cout << " " << c2n_n[3] << std::endl;
#endif
    status = wmesh_def(mesh__,
		       num_nodes,
		       4,
		       c2n_ptr,
		       c2n_m,
		       c2n_n,
		       c2n_v,
		       c2n_ld,
		       coo,
		       3);
    
    return WMESH_STATUS_SUCCESS;
  };
#endif
  
};
