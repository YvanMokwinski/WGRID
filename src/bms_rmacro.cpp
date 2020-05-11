#include <stdlib.h>
#include <iostream>
#include "wmesh.hpp"
#include "wmesh_utils.hpp"

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
  WMESH_CHECK_POINTER(tr2n_);
  WMESH_CHECK_POINTER(tet2n_);
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


static inline wmesh_status_t bms_rmacro_c2n(wmesh_int_t 	element_,
					    wmesh_int_t 	degree_,
					    wmesh_int_t 	c2n_size_,
					    const_wmesh_int_p 	c2n_ptr_,
					    const_wmesh_int_p 	c2n_m_,
					    const_wmesh_int_p 	c2n_n_,
					    wmesh_int_p 	c2n_v_,
					    const_wmesh_int_p 	c2n_ld_,
					    wmesh_int_t 	work_n_,
					    wmesh_int_p 	work_)
{
  WMESH_CHECK( element_ > 0 );
  WMESH_CHECK_POINTER(c2n_ptr_);
  WMESH_CHECK_POINTER(c2n_m_);
  WMESH_CHECK_POINTER(c2n_n_);
  WMESH_CHECK_POINTER(c2n_v_);
  WMESH_CHECK_POINTER(c2n_ld_);

  wmesh_int_t s_entity_local_index_of_topodim[WMESH_ELEMENT_ALL] = {0,0,0,1,0,1,2,3};  
  wmesh_int_t cell_type_ = s_entity_local_index_of_topodim[element_];
  
  //
  // 0 0
  // 1 0
  // 2 0
  // 3 1
  // 4 0
  // 5 1
  // 6 2
  // 7 3
  //
  const wmesh_int_t n1d = degree_+1;
  if (cell_type_ == 0 && work_n_ < n1d * n1d * 3)
    {
      std::cerr << "cell_type_ = " << cell_type_ << ", work_n_ = " << work_n_ << ", required_work_n_ " << n1d * n1d * 3 << std::endl;
      return WMESH_STATUS_ERROR_WORKSPACE;      
    }
  
  WMESH_CHECK_POINTER(c2n_v_);

  wmesh_status_t status;
  switch(element_)
    {
      
    case WMESH_ELEMENT_EDGE:
      {
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_TRIANGLE:
      {
#define _dec(_i,_j) (( (_i)+(_j) + 1 )*( (_i)+(_j) ))/2+(_i)
	wmesh_int_p c2n		= c2n_v_ + c2n_ptr_[cell_type_];
	wmesh_int_t c2n_shift 	= c2n_ld_[cell_type_];
	  for (wmesh_int_t j=0;j<degree_;j++)
	    {	    
	      for (wmesh_int_t i=0;i<j;i++)
		{		  
		  c2n[0] = _dec(i,j-i);
		  c2n[1] = _dec(i+1,j-i);
		  c2n[2] = _dec(i,j+1-i); 
		  c2n += c2n_shift;
		  
		  c2n[0] = _dec(i,j-i);
		  c2n[1] = _dec(i+1,j-1-i);
		  c2n[2] = _dec(i+1,j-i);
		  c2n += c2n_shift;
		}
	      
	      c2n[0] = _dec(j,0);
	      c2n[1] = _dec(j+1,0);
	      c2n[2] = _dec(j,1);
	      c2n += c2n_shift;
	    }
#undef _dec  	
	  return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_QUADRILATERAL:
      {
#define _dec(_i,_j) (degree_+1) * (_i) + (_j)
	wmesh_int_p c2n		= c2n_v_ + c2n_ptr_[cell_type_];
	wmesh_int_t c2n_shift 	= c2n_ld_[cell_type_];
	for (wmesh_int_t i=0;i<degree_;i++)
	  {
	    for (wmesh_int_t j=0;j<degree_;j++)
	      {
		c2n[0] = _dec(i,j);
		c2n[1] = _dec(i+1,j);
		c2n[2] = _dec(i+1,j+1);
		c2n[3] = _dec(i,j+1);
		c2n += c2n_shift;
	      } 
	  }
#undef _dec
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_TETRAHEDRON:
      {
	//
	// Needs workspace.
	//
	WMESH_CHECK_POINTER(work_);      
#define _dec(_i,_j) ((  (_i)+(_j) + 1 )*( (_i)+(_j) ))/2+(_i)
	wmesh_int_p p	= work_;      
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
	    }    
#undef _dec
	  
	  status  = extrude_tetra(degree_,			     
				  work_,
				  3,
				  c2n_v_ + c2n_ptr_[cell_type_],
				  c2n_ld_[cell_type_]);
	  
	  WMESH_STATUS_CHECK(status);
	  
	  return WMESH_STATUS_SUCCESS;
      }
    case WMESH_ELEMENT_PYRAMID:
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
    case WMESH_ELEMENT_WEDGE:
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
    case WMESH_ELEMENT_HEXAHEDRON:
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
		    c2n += c2n_shift;
		  } 
	      }
	  }
#undef _dec
	return WMESH_STATUS_SUCCESS;
      }
    case WMESH_ELEMENT_NODE:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	break;
      }
    default:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	break;
      }
    }
  
  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
}

static inline wmesh_status_t bms_rmacro_internal_icoo(wmesh_int_t element_,
						      wmesh_int_t degree_,
						      wmesh_int_p icoo_,
						      wmesh_int_t icoo_ld_) 
{
  WMESH_CHECK_POINTER(icoo_);
  wmesh_int_t n1d = degree_+1;
  switch(element_)
    {      
    case WMESH_ELEMENT_EDGE:
      {
	wmesh_int_t idx = 0;
	for (wmesh_int_t i=0;i<n1d;i++)
	  {
	    icoo_[icoo_ld_*idx+0] = i;
	    ++idx;
	  }
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_TRIANGLE:
      {
	wmesh_int_p icoo = icoo_;
	for (wmesh_int_t i=0;i<n1d;i++)
	  {
	    for (wmesh_int_t j=0;j<=i;j++)
	      {
		icoo[0] = j;
		icoo[1] = i-j;
		icoo += icoo_ld_;
	      } 
	  }
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_QUADRILATERAL:
      {
	wmesh_int_p icoo = icoo_;
	for (wmesh_int_t i=0;i<n1d;++i)
	  {		
	    for (wmesh_int_t j=0;j<n1d;j++)
	      {
		icoo[0] = i;
		icoo[1] = j;
		icoo += icoo_ld_;
	      }
	  }
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_TETRAHEDRON:
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
    case WMESH_ELEMENT_PYRAMID:
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
    case WMESH_ELEMENT_WEDGE:
      {
	//
	//  6 
	//  3 7
	//  1 4 8 
	//  0 2 5 9
	//
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
    case WMESH_ELEMENT_HEXAHEDRON:
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
    case WMESH_ELEMENT_NODE:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	break;
      }
    default:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	break;
      }
    }
  
  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
}


wmesh_status_t bms_fe_triangle_ordering(wmesh_int_p 	icoo_,
					wmesh_int_t 	icoo_ld_,
					wmesh_int_t 	shift_,
					wmesh_int_t 	degree_)
{
  WMESH_CHECK_POINTER(icoo_);
  
  if (degree_ == 0)
    {
      wmesh_int_p icoo = icoo_;  
      icoo[0] = shift_;
      icoo[1] = shift_;
      icoo += icoo_ld_;
      return WMESH_STATUS_SUCCESS;
    }
  
  const wmesh_int_t n = (degree_>0) ? (degree_)-1 : 0;
  wmesh_int_p icoo = icoo_;
  
  icoo[0] = shift_ + 0;
  icoo[1] = shift_ + 0;
  icoo += icoo_ld_;
  
  icoo[0] = shift_ + degree_;
  icoo[1] = shift_ + 0;
  icoo += icoo_ld_;
  
  icoo[0] = shift_ + 0;
  icoo[1] = shift_ + degree_;
  icoo += icoo_ld_;
  
  for (wmesh_int_t i=0;i<n;++i)
    {
      icoo[0] = shift_+ i + 1;
      icoo[1] = shift_+ 0;
      icoo += icoo_ld_;
    } 
  
  for (wmesh_int_t i=0;i<n;++i)
    {
      icoo[0] = shift_+ degree_ - (i + 1);
      icoo[1] = shift_+ i + 1;
      icoo += icoo_ld_;
    }    
  
  for (wmesh_int_t i=0;i<n;++i)
    {
      icoo[0] = shift_ + 0;
      icoo[1] = shift_ + degree_ - (i + 1);
      icoo += icoo_ld_;
    }

  if (degree_>=3)
    {
  wmesh_status_t status =  bms_fe_triangle_ordering(icoo,
    icoo_ld_,
    shift_ + 1,
    degree_-3);
  WMESH_STATUS_CHECK(status);
    }
  
  return WMESH_STATUS_SUCCESS;
}

wmesh_status_t bms_fe_quadrilateral_ordering(wmesh_int_p 	icoo_,
					     wmesh_int_t 	icoo_ld_,
					     wmesh_int_t 	shift_,
					     wmesh_int_t 	degree_)
{
  WMESH_CHECK_POINTER(icoo_);

  wmesh_int_p 		icoo 	= icoo_;
  
  if (degree_ == 0)
    {
      icoo[0] = shift_;
      icoo[1] = shift_;
      icoo += icoo_ld_;
      return WMESH_STATUS_SUCCESS;
    }
  
  
  icoo[0] = shift_ + 0;
  icoo[1] = shift_ + 0;
  icoo += icoo_ld_;
  
  icoo[0] = shift_ + degree_;
  icoo[1] = shift_ + 0;
  icoo += icoo_ld_;
  
  icoo[0] = shift_ + degree_;
  icoo[1] = shift_ + degree_;
  icoo += icoo_ld_;

  icoo[0] = shift_ + 0;
  icoo[1] = shift_ + degree_;
  icoo += icoo_ld_;

  const wmesh_int_t n = (degree_>0) ? degree_ - 1 : 0;
  for (wmesh_int_t i=0;i<n;++i)
    {
      icoo[0] = shift_+ i + 1;
      icoo[1] = shift_+ 0;
      icoo += icoo_ld_;
    } 
  
  for (wmesh_int_t i=0;i<n;++i)
    {
      icoo[0] = shift_+ degree_;
      icoo[1] = shift_+ i + 1;
      icoo += icoo_ld_;
    }    

  for (wmesh_int_t i=0;i<n;++i)
    {
      icoo[0] = shift_+ n - i;
      icoo[1] = shift_+ degree_;
      icoo += icoo_ld_;
    }    

  for (wmesh_int_t i=0;i<n;++i)
    {
      icoo[0] = shift_ + 0;
      icoo[1] = shift_ + n - i;
      icoo += icoo_ld_;
    }

  if (degree_>=2)
    {
  wmesh_status_t status =  bms_fe_quadrilateral_ordering(icoo,
    icoo_ld_,
    shift_ + 1,
    degree_-2);
  WMESH_STATUS_CHECK(status);
    }
  
  return WMESH_STATUS_SUCCESS;
}
	
static wmesh_status_t bms_rmacro_fe_icoo(wmesh_int_t element_,
				     wmesh_int_t degree_,
				     wmesh_int_p icoo_,
				     wmesh_int_t icoo_ld_) 
{
  WMESH_CHECK_POINTER(icoo_);
  wmesh_int_t n1d = degree_+1;
  switch(element_)
    {
      
    case WMESH_ELEMENT_EDGE:
      {
	wmesh_int_p icoo = icoo_;
	icoo[0] = 0; icoo += icoo_ld_;
	icoo[0] = degree_; icoo += icoo_ld_;
	for (wmesh_int_t i=1;i<degree_;++i)
	  {
	    icoo[0] = i; icoo += icoo_ld_;
	  }
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_TRIANGLE:
      {
	bms_fe_triangle_ordering(icoo_,icoo_ld_,0,degree_);
#if 0
	wmesh_int_p icoo = icoo_;
	for (wmesh_int_t d = 0;d<degree_;++d)
	  {  
	    const wmesh_int_t n = (degree_-d>0)? (degree_-d)-1 : 0;
	    
	    icoo[0] = d;
	    icoo[1] = d;
	    icoo += icoo_ld_;
	    
	    icoo[0] = degree_ - d;
	    icoo[1] = d;
	    icoo += icoo_ld_;
	    
	    icoo[0] = d;
	    icoo[1] = degree_ - d;
	    icoo += icoo_ld_;
	    
	    for (wmesh_int_t i=d;i<n;++i)
	      {
		icoo[0] = i+1;
		icoo[1] = d;
		icoo += icoo_ld_;
	      } 
	    
	    for (wmesh_int_t i=d;i<n;++i)
	      {
		icoo[0] = (degree_ - d) - (i+1-d);
		icoo[1] = i + 1;
		icoo += icoo_ld_;
	      }    
	    
	    for (wmesh_int_t i=d;i<n;++i)
	      {
		icoo[0] = d + 0;
		icoo[1] = (degree_-d) - (i+1-d);
		icoo += icoo_ld_;
	      }

	  }
#endif	  

#if 0
	const wmesh_int_t n = (degree_>0)? degree_-1 : 0;

	wmesh_int_p icoo = icoo_;
	icoo[0] = 0;
	icoo[1] = 0;
	icoo += icoo_ld_;

	icoo[0] = degree_;
	icoo[1] = 0;
	icoo += icoo_ld_;

	icoo[0] = 0;
	icoo[1] = degree_;
	icoo += icoo_ld_;

	
	for (wmesh_int_t i=0;i<n;++i)
	  {
	    icoo[0] = i+1;
	    icoo[1] = 0;
	    icoo += icoo_ld_;
	  } 

	for (wmesh_int_t i=0;i<n;++i)
	  {
	    icoo[0] = degree_ - (i+1);
	    icoo[1] = i+1;
	    icoo += icoo_ld_;
	  } 	

	for (wmesh_int_t i=0;i<n;++i)
	  {
	    icoo[0] = 0;
	    icoo[1] = degree_ - (i+1);
	    icoo += icoo_ld_;
	  } 
	
	for (wmesh_int_t i=0;i<n1d;i++)
	  {
	    for (wmesh_int_t j=0;j<=i;j++)
	      {
		icoo[0] = j;
		icoo[1] = i - j;
		icoo += icoo_ld_;
	      } 
	  }
#endif
	
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_QUADRILATERAL:
      {
	bms_fe_quadrilateral_ordering(icoo_,icoo_ld_,0,degree_);
#if 0
	wmesh_int_p icoo = icoo_;
	for (wmesh_int_t d = 0;d<degree_;++d)
	  {  
	    const wmesh_int_t n = (degree_-d>0)? (degree_-d)-1 : 0;
	    
	    icoo[0] = d;
	    icoo[1] = d;
	    icoo += icoo_ld_;
	    
	    icoo[0] = degree_ - d;
	    icoo[1] = d;
	    icoo += icoo_ld_;
	    
	    icoo[0] = degree_ - d;
	    icoo[1] = degree_ - d;
	    icoo += icoo_ld_;
	    
	    icoo[0] = d;
	    icoo[1] = degree_ - d;
	    icoo += icoo_ld_;
	    
	    for (wmesh_int_t i=d;i<n;++i)
	      {
		icoo[0] = i+1;
		icoo[1] = d;
		icoo += icoo_ld_;
	      } 
	    
	    for (wmesh_int_t i=d;i<n;++i)
	      {
		icoo[0] = degree_ - d;
		icoo[1] = i + 1;
		icoo += icoo_ld_;
	      }    
	    
	    for (wmesh_int_t i=d;i<n;++i)
	      {
     	    icoo[0] = (degree_ - d) - (i+1-d);
	    icoo[1] = (degree_ - d);
	    icoo += icoo_ld_;
	  }
	    
	    for (wmesh_int_t i=d;i<n;++i)
	      {
		icoo[0] = d + 0;
		icoo[1] = (degree_-d) - (i+1-d);
		icoo += icoo_ld_;
	      }

	  }
	  
	  if (degree_ %2 ==0)
	    {
	    icoo[0] = degree_/2;
	    icoo[1] = degree_/2;
	    std::cout << icoo[0]<< " ddddd " << icoo[1]<<std::endl;
	    icoo += icoo_ld_;	    
	  }
#endif	  
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_TETRAHEDRON:
      {
	wmesh_int_t idx = 0;
	static constexpr const wmesh_int_t reftetra_icoo[] = {0,0,0,
							      1,0,0,
							      0,1,0,
							      0,0,1};
      
	const wmesh_int_t n = (degree_>0)? degree_-1 : 0;
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
    case WMESH_ELEMENT_PYRAMID:
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

    case WMESH_ELEMENT_WEDGE:
      {
	const wmesh_int_t	n = (degree_>0)? degree_-1 : 0;
	wmesh_int_t idx = 0;

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

    case WMESH_ELEMENT_HEXAHEDRON:
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
    case WMESH_ELEMENT_NODE:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	break;
      }
    default:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	break;
      }
    }
  
  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);

  return WMESH_STATUS_INVALID_ARGUMENT;

};	  


extern "C"
{
  
  wmesh_status_t bms_rmacro_buffer_size(wmesh_int_t 	element_,
					wmesh_int_t 	degree_,
					wmesh_int_p	work_n_,
					wmesh_int_p	num_entities_)
  {
    WMESH_CHECK_POINTER(work_n_);
    WMESH_CHECK_POINTER(num_entities_);        
    work_n_[0] = 0;
    for (int i=0;i<WMESH_ELEMENT_ALL;++i)
      {
	num_entities_[i] = 0; 
      }

    switch(element_)
      {
      case WMESH_ELEMENT_EDGE:
	{
	  num_entities_[WMESH_ELEMENT_NODE] = degree_+1;
	  num_entities_[WMESH_ELEMENT_TRIANGLE] = degree_;
	  break;
	}
      case WMESH_ELEMENT_TRIANGLE:
	{
	  num_entities_[WMESH_ELEMENT_NODE] = ( (degree_+1)*(degree_+2) ) / 2;
	  num_entities_[WMESH_ELEMENT_TRIANGLE] = degree_*degree_;
	  break;
	}
      case WMESH_ELEMENT_QUADRILATERAL:
	{
	  num_entities_[WMESH_ELEMENT_NODE] = ( (degree_+1)*(degree_+1));
	  num_entities_[WMESH_ELEMENT_QUADRILATERAL] = degree_*degree_;
	  break;
	}
      case WMESH_ELEMENT_TETRAHEDRON:
	{
	  num_entities_[WMESH_ELEMENT_NODE] =  ( (degree_+1)*(degree_+2)*(degree_+3))/6;	
	  num_entities_[WMESH_ELEMENT_TETRAHEDRON] = degree_*degree_*degree_;
	  break;
	}
      case WMESH_ELEMENT_PYRAMID:
	{
	  num_entities_[WMESH_ELEMENT_NODE]  = ((degree_+1)*(degree_+2)*(2*degree_+3))/6;
	  num_entities_[WMESH_ELEMENT_PYRAMID] = ((degree_-1)*(degree_)*(2*degree_-1))/3 + degree_ * degree_;
	  num_entities_[WMESH_ELEMENT_TETRAHEDRON] = (2*(degree_-1)*(degree_)*(degree_+1))/3;
	  break;
	}
      case WMESH_ELEMENT_WEDGE:
	{
	  num_entities_[WMESH_ELEMENT_NODE] = ( ( (degree_+1) *(degree_+2) ) /2 ) * (degree_+1);
	  num_entities_[WMESH_ELEMENT_WEDGE] = degree_*degree_*degree_;
	  break;
	}
      case WMESH_ELEMENT_HEXAHEDRON:
	{
	  num_entities_[WMESH_ELEMENT_NODE] = ( (degree_+1)*(degree_+1)*(degree_+1));
	  num_entities_[WMESH_ELEMENT_HEXAHEDRON] = degree_*degree_*degree_;
	  break;
	}
      case WMESH_ELEMENT_NODE:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	  break;
	}
      default:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	  break;
	}
      }
    
    work_n_[0] = 4 * (degree_+1)*(degree_+1)*(degree_+1) + 3 * (degree_+1)*(degree_+1) + num_entities_[WMESH_ELEMENT_NODE] * 3;
    return WMESH_STATUS_SUCCESS;
  };
  
  
  wmesh_status_t bms_rmacro(    wmesh_int_t 		element_,
				wmesh_int_t 		degree_,
				wmesh_int_t 		num_nodes_,

				wmesh_int_t 		c2n_size_,
				const_wmesh_int_p	c2n_ptr_,
				const_wmesh_int_p	c2n_m_,
				const_wmesh_int_p	c2n_n_,
				wmesh_int_p 		c2n_v_,
				const_wmesh_int_p	c2n_ld_,
				
				wmesh_int_t 		icoo_m_,
				wmesh_int_t 		icoo_n_,
				wmesh_int_p 		icoo_v_,
				wmesh_int_t 		icoo_ld_,
				
				wmesh_int_t		work_n_,
				wmesh_int_p		work_)
  {

    
    wmesh_status_t status;

    //
    // Topoloigical dimension.
    //
    wmesh_int_t topodim;
    status = wmesh_element2topodim(element_,
				   &topodim);
    WMESH_STATUS_CHECK(status);
    
    wmesh_int_t
      perm_n,
      ilagr_n,
      work_n;
    
    perm_n = (degree_+1)*(degree_+1)*(degree_+1);
    ilagr_n= (degree_+1)*(degree_+1)*(degree_+1)*3;
    work_n = 3*(degree_+1)*(degree_+1);
#if 0
    switch(topodim)
      {
      case 3:
	{
	  break;
	}
	
      case 2:
	{
	  perm_n = (degree_+1)*(degree_+1);
	  ilagr_n= (degree_+1)*(degree_+1)*2;
	  work_n = 2*(degree_+1);
	  break;	  
	}
	
      case 1:
	{
	  perm_n  = (degree_+1);
	  ilagr_n = (degree_+1);
	  work_n  = (degree_+1);
	  break;
	}
	

      default:
      case 0:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
	}

      }
#endif
    wmesh_int_p ilagr 	= work_;
    wmesh_int_p perm 	= work_ + (ilagr_n);
    wmesh_int_p work 	= work_ + (ilagr_n + perm_n);
    
    status 		= bms_rmacro_fe_icoo(element_,
					 degree_,
					 icoo_v_,
					 icoo_ld_);    
    WMESH_STATUS_CHECK(status);

    const wmesh_int_t N = (degree_+1);
    const wmesh_int_t NxN = N*N;
    wmesh_int_t ijk[3];
    //
    // Compute the permutation.
    //
    for (wmesh_int_t idx=0;idx<num_nodes_;++idx)
      {
	for (wmesh_int_t idim=0;idim<topodim;++idim)
	  {
	    ijk[idim] = icoo_v_[icoo_ld_*idx + idim];
	  }
	if (topodim==3)
	  perm[NxN * ijk[2] + N * ijk[0] + ijk[1]] = 1+idx;
	else if (topodim==2)
	  perm[N * ijk[0] + ijk[1]] = 1+idx;
	else if (topodim==1)
	  perm[ijk[0]] = 1+idx;
      }
    
    //
    // Compute the connectivity.
    //
    status = bms_rmacro_c2n(element_,
			    degree_,
			    c2n_size_,
			    c2n_ptr_,
			    c2n_m_,
			    c2n_n_,
			    c2n_v_,
			    c2n_ld_,
			    work_n,
			    work);    
    WMESH_STATUS_CHECK(status);
    
    //
    // Compute the integers.
    //
    status = bms_rmacro_internal_icoo(element_,
				      degree_,
				      ilagr,
				      topodim);
    WMESH_STATUS_CHECK(status);

    if (topodim==3)
      {
	//
	// Adjust the connectivity.
	//
	const wmesh_int_t k_ld = (degree_+1)*(degree_+1);
	const wmesh_int_t i_ld = (degree_+1);
	for (wmesh_int_t subcell_type = 0;subcell_type < c2n_size_;++subcell_type)
	  {
	    const wmesh_int_t num_subcells = c2n_n_[subcell_type];
	    if (num_subcells > 0)
	      {
		wmesh_int_p c2n_v 				= c2n_v_ + c2n_ptr_[subcell_type];
		const wmesh_int_t c2n_ld 			= c2n_ld_[subcell_type];
		const wmesh_int_t num_nodes_in_subcell 	= c2n_m_[subcell_type];
		for (wmesh_int_t j = 0;j < num_subcells;++j)
		  {
		    for (wmesh_int_t i=0;i<num_nodes_in_subcell;++i)
		      {	  
			const wmesh_int_t l		= c2n_v[c2n_ld * j + i];
			const wmesh_int_t ip	= ilagr[3*l+0];
			const wmesh_int_t jp	= ilagr[3*l+1];
			const wmesh_int_t kp	= ilagr[3*l+2];
			c2n_v[c2n_ld * j + i]  	= perm[k_ld*kp+i_ld*ip+jp];
		      }
		  }
	      }
	  }
      }
    else if (topodim==2)
      {
	//
	// Adjust the connectivity.
	//	
	const wmesh_int_t i_ld = (degree_+1);
	for (wmesh_int_t subcell_type = 0;subcell_type < c2n_size_;++subcell_type)
	  {
	    const wmesh_int_t num_subcells = c2n_n_[subcell_type];
	    if (num_subcells > 0)
	      {
		wmesh_int_p c2n_v 				= c2n_v_ + c2n_ptr_[subcell_type];
		const wmesh_int_t c2n_ld 			= c2n_ld_[subcell_type];
		const wmesh_int_t num_nodes_in_subcell 	= c2n_m_[subcell_type];
		for (wmesh_int_t j = 0;j < num_subcells;++j)
		  {
		    for (wmesh_int_t i=0;i<num_nodes_in_subcell;++i)
		      {	  
			const wmesh_int_t l	= c2n_v[c2n_ld * j + i];
			const wmesh_int_t ip	= ilagr[topodim*l+0];
			const wmesh_int_t jp	= ilagr[topodim*l+1];
			c2n_v[c2n_ld * j + i]  	= perm[i_ld*ip+jp];
		      }
		  }
	      }
	  }
      }
    else
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
      }
    return WMESH_STATUS_SUCCESS;
  };



  

  
};
