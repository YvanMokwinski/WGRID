#include "bms.h"
#include "wmesh-status.h"
#include "wmesh-enums.h"
#include "bms_templates.hpp"
#include <string.h>

template<typename T>
wmesh_status_t bms_ordering_vertices(wmesh_int_t 	element_,
				     T*__restrict__ 	c_)
{
  switch(element_)
    {      
    case WMESH_ELEMENT_NODE:
      {
	static constexpr T ref[1] = {0};
	memcpy(c_,ref,sizeof(T)*1);
	return WMESH_STATUS_SUCCESS;
      }
    case WMESH_ELEMENT_EDGE:
      {
	static constexpr T ref[2] = {0,1};
	memcpy(c_,ref,sizeof(T)*2);
	return WMESH_STATUS_SUCCESS;
      }
    case WMESH_ELEMENT_TRIANGLE:
      {
	static constexpr T ref[6] = {0,0,1,0,0,1};
	memcpy(c_,ref,sizeof(T)*6);
	return WMESH_STATUS_SUCCESS;
      }
    case WMESH_ELEMENT_QUADRILATERAL:
      {
	static constexpr T ref[8] = {0,0,1,0,1,1,0,1};
	memcpy(c_,ref,sizeof(T)*8);
	return WMESH_STATUS_SUCCESS;
      }
    case WMESH_ELEMENT_TETRAHEDRON:
      {
	static constexpr T ref[12] = {0,0,0,
				      1,0,0,
				      0,1,0,
				      0,0,1};
	memcpy(c_,ref,sizeof(T)*12);
	return WMESH_STATUS_SUCCESS;
      }
    case WMESH_ELEMENT_PYRAMID:
      {
	static constexpr T ref[15] = {0,0,0,1,0,0,1,1,0,0,1,0,0,0,1};
	memcpy(c_,ref,sizeof(T)*15);
	return WMESH_STATUS_SUCCESS;
      }
    case WMESH_ELEMENT_WEDGE:
      {
	static constexpr T ref[18] = {0,0,0,
				      1,0,0,
				      0,1,0,
				      0,0,1,
				      1,0,1,
				      0,1,1};	  
	memcpy(c_,ref,sizeof(T)*18);
	return WMESH_STATUS_SUCCESS;
      }
    case WMESH_ELEMENT_HEXAHEDRON:
      {
	static constexpr     T ref[24] = {0,0,0,
					  1,0,0,
					  1,1,0,
					  0,1,0,
					  0,0,1,
					  1,0,1,
					  1,1,1,
					  0,1,1};
	
	memcpy(c_,ref,sizeof(T)*24);
	return WMESH_STATUS_SUCCESS;
      }


    }
  return WMESH_STATUS_INVALID_ENUM;
    
}


template
wmesh_status_t bms_ordering_vertices<wmesh_int_t>(wmesh_int_t 	element_,
						  wmesh_int_t*__restrict__ 	c_);
template
wmesh_status_t bms_ordering_vertices<float>	(wmesh_int_t 	element_,
						 float*__restrict__ 	c_);
template
wmesh_status_t bms_ordering_vertices<double>	(wmesh_int_t 	element_,
					     double*__restrict__ 	c_);


extern "C"
{
  
  static wmesh_status_t bms_ordering_topoid_calculate(wmesh_int_t 	num_dofs_0_,
						      wmesh_int_t 	num_dofs_1_,
						      wmesh_int_t 	num_dofs_2t_,
						      wmesh_int_t 	num_dofs_2q_,
						      wmesh_int_t 	num_dofs_3_,
						      wmesh_int_t 	topoid_n_,
						      wmesh_int_p 	topoid_v_,
						      wmesh_int_t 	topoid_inc_)
    
  {
    
    for (wmesh_int_t i=0;i<num_dofs_0_;++i)
      {
	topoid_v_[0] = 0;
	topoid_v_ += topoid_inc_;
      }

    for (wmesh_int_t i=0;i<num_dofs_1_;++i)
      {
	topoid_v_[0] = 1;
	topoid_v_ += topoid_inc_;
      }
    
    for (wmesh_int_t i=0;i<num_dofs_2t_;++i)
      {
	topoid_v_[0] = 2;
	topoid_v_ += topoid_inc_;
      }
    
    for (wmesh_int_t i=0;i<num_dofs_2q_;++i)
      {
	topoid_v_[0] = 2;
	topoid_v_ += topoid_inc_;
      }
    
    for (wmesh_int_t i=0;i<num_dofs_3_;++i)
      {
	topoid_v_[0] = 3;
	topoid_v_ += topoid_inc_;
      }
    
    return WMESH_STATUS_SUCCESS;
  }




  
  wmesh_status_t bms_ordering_triangle(wmesh_int_t 	degree_,
				       wmesh_int_t 	c_storage_,
				       wmesh_int_t 	c_m_,
				       wmesh_int_t 	c_n_,
				       wmesh_int_p 	c_,
				       wmesh_int_t 	c_ld_,
				       wmesh_int_t 	shift_)
  {
    wmesh_int_t ref_icoo[6];
    {
      wmesh_status_t status = bms_ordering_vertices(WMESH_ELEMENT_TRIANGLE,
						    ref_icoo);
      WMESH_STATUS_CHECK(status);
    }
      
#if 0
    static constexpr   wmesh_int_t ref_icoo[] 	= {0,0,
						   1,0,
						   0,1};
#endif
    
    static constexpr   wmesh_int_t ref_icoo_n 	= 3;
    static constexpr   wmesh_int_t ref_icoo_ld 	= 2;
    
    WMESH_CHECK_POINTER(c_);    
    
    wmesh_int_p ic = c_;  
    wmesh_int_p jc = c_ + ( (WMESH_STORAGE_INTERLEAVE == c_storage_) ? 1 : c_ld_ );  
    if (degree_ == 0)
      {
	ic[0] = shift_;
	jc[0] = shift_;
	return WMESH_STATUS_SUCCESS;
      }
    
    WMESH_CHECK_POSITIVE(degree_);    
    
    const wmesh_int_t n 	= degree_ - 1;
    const wmesh_int_t ic_ld 	= (WMESH_STORAGE_INTERLEAVE == c_storage_) ? c_ld_ : 1;
    const wmesh_int_t jc_ld 	= ic_ld;
    
    for (wmesh_int_t i=0;i<ref_icoo_n;++i)
      {
	*ic = shift_ + ref_icoo[ref_icoo_ld*i+0] * degree_;
	*jc = shift_ + ref_icoo[ref_icoo_ld*i+1] * degree_;
	ic += ic_ld;
	jc += jc_ld;
      }
  
    for (wmesh_int_t i=0;i<n;++i)
      {
	ic[0] = shift_+ i + 1;
	jc[0] = shift_+ 0;
	ic += ic_ld;
	jc += jc_ld;      
      } 
  
    for (wmesh_int_t i=0;i<n;++i)
      {
	ic[0] = shift_+ degree_ - (i + 1);
	jc[0] = shift_+ i + 1;
	ic += ic_ld;
	jc += jc_ld;      
      }    
  
    for (wmesh_int_t i=0;i<n;++i)
      {
	ic[0] = shift_ + 0;
	jc[0] = shift_ + degree_ - (i + 1);
	ic += ic_ld;
	jc += jc_ld;      
      }
  
    if (degree_ >= 3)
      {
	wmesh_status_t status =  bms_ordering_triangle(degree_ - 3,
						       c_storage_,
						       c_m_,
						       c_n_,						       
						       ic,
						       c_ld_,
						       shift_ + 1);
	WMESH_STATUS_CHECK(status);
      }
  
    return WMESH_STATUS_SUCCESS;
  }

  

  wmesh_status_t bms_ordering_quadrilateral(wmesh_int_t 	degree_,
					    wmesh_int_t 	c_storage_,
					    wmesh_int_t		c_m_,
					    wmesh_int_t		c_n_,
					    wmesh_int_p 	c_,
					    wmesh_int_t 	c_ld_,
					    wmesh_int_t 	shift_)
  {
#if 0
    static constexpr wmesh_int_t ref_icoo[] = {0,0,
					       1,0,
					       1,1,
					       0,1};
#endif
    wmesh_int_t ref_icoo[8];
    {
      wmesh_status_t status = bms_ordering_vertices(WMESH_ELEMENT_QUADRILATERAL,
						    ref_icoo);
      WMESH_STATUS_CHECK(status);
    }

    static constexpr wmesh_int_t ref_icoo_n = 4;
    static constexpr wmesh_int_t ref_icoo_ld = 2;
    wmesh_status_t status;
    
    WMESH_CHECK_POINTER(c_);
    
    wmesh_int_p ic 		= c_;  
    wmesh_int_p jc 	    	= c_ + ( (WMESH_STORAGE_INTERLEAVE == c_storage_) ? 1 : c_ld_ );  

    const wmesh_int_t ic_ld 	= (WMESH_STORAGE_INTERLEAVE == c_storage_) ? c_ld_ : 1;
    const wmesh_int_t jc_ld 	= ic_ld;
 
    if (degree_ == 0)
      {
	*ic = shift_;
	*jc = shift_;
	return WMESH_STATUS_SUCCESS;
      }
    
    WMESH_CHECK_POSITIVE(degree_);
    const wmesh_int_t n	= degree_ - 1;
    
    for (wmesh_int_t i=0;i<ref_icoo_n;++i)
      {
	*ic = shift_ + ref_icoo[ref_icoo_ld*i+0] * degree_;
	*jc = shift_ + ref_icoo[ref_icoo_ld*i+1] * degree_;
	ic += ic_ld;
	jc += jc_ld;
      }

    for (wmesh_int_t i=0;i<n;++i)
      {
	*ic = shift_+ i + 1;
	*jc = shift_+ 0;
	ic += ic_ld;
	jc += jc_ld;      
      } 
  
    for (wmesh_int_t i=0;i<n;++i)
      {
	*ic = shift_+ degree_;
	*jc = shift_+ i + 1;
	ic += ic_ld;
	jc += jc_ld;      
      }    

    for (wmesh_int_t i=0;i<n;++i)
      {
	*ic = shift_+ n - i;
	*jc = shift_+ degree_;
	ic += ic_ld;
	jc += jc_ld;      
      }    

    for (wmesh_int_t i=0;i<n;++i)
      {
	*ic = shift_ + 0;
	*jc = shift_ + n - i;
	ic += ic_ld;
	jc += jc_ld;      
      }

    if (degree_>=2)
      {
	status =  bms_ordering_quadrilateral(degree_ - 2,
					     c_storage_,
					     c_m_,
					     c_n_,
					     ic,
					     c_ld_,
					     shift_ + 1);
	WMESH_STATUS_CHECK(status);
      }
  
    return WMESH_STATUS_SUCCESS;
  }



  wmesh_status_t bms_ordering_edge(wmesh_int_t		degree_,
				   wmesh_int_t		c_storage_,
				   wmesh_int_t		c_m_,
				   wmesh_int_t		c_n_,
				   wmesh_int_p		c_,
				   wmesh_int_t		c_ld_)
  {
    WMESH_CHECK_POINTER(c_);
    wmesh_int_p c = c_;

    c[0] = 0;
    c += c_ld_;
    
    c[0] = degree_;
    c += c_ld_;
    
    for (wmesh_int_t i = 1;i < degree_;++i)
      {
	c[0] = i;
	c += c_ld_;
      }
    return WMESH_STATUS_SUCCESS;
  }

  wmesh_status_t bms_ordering_face(wmesh_int_t		element_,
				   wmesh_int_t		degree_,
				   wmesh_int_t		c_storage_,
				   wmesh_int_t		c_m_,
				   wmesh_int_t		c_n_,
				   wmesh_int_p		c_v_,
				   wmesh_int_t		c_ld_) 
  {  
    wmesh_status_t status;
    WMESH_CHECK_POINTER(c_v_);
    switch(element_)
      {      
      case WMESH_ELEMENT_TRIANGLE:
	{
	  static constexpr wmesh_int_t shift = 0;
	  status = bms_ordering_triangle(degree_,
					 c_storage_,
					 c_m_,
					 c_n_,
					 c_v_,
					 c_ld_,
					 shift);
	  WMESH_STATUS_CHECK(status);
	  return WMESH_STATUS_SUCCESS;
	}
      
      case WMESH_ELEMENT_QUADRILATERAL:
	{
	  static constexpr wmesh_int_t shift = 0;
	  status = bms_ordering_quadrilateral(degree_,
					      c_storage_,
					      c_m_,
					      c_n_,
					      c_v_,
					      c_ld_,
					      shift);
	  WMESH_STATUS_CHECK(status);
	  return WMESH_STATUS_SUCCESS;
	}

      case WMESH_ELEMENT_NODE:
      case WMESH_ELEMENT_EDGE:
      case WMESH_ELEMENT_TETRAHEDRON:
      case WMESH_ELEMENT_PYRAMID:
      case WMESH_ELEMENT_WEDGE:
      case WMESH_ELEMENT_HEXAHEDRON:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
	  break;
	}
      }  
    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
  };	  
 

  wmesh_status_t bms_ordering_volume(wmesh_int_t	element_,
				     wmesh_int_t	degree_,
				     wmesh_int_t	c_storage_,
				     wmesh_int_t	c_m_,
				     wmesh_int_t	c_n_,
				     wmesh_int_p	c_,
				     wmesh_int_t	c_ld_,
				     
				     wmesh_int_t 	s_e2n_m_,
				     wmesh_int_t 	s_e2n_n_,
				     const_wmesh_int_p 	s_e2n_v_,
				     wmesh_int_t 	s_e2n_ld_,
				     
				     wmesh_int_t 	s_t2n_m_,
				     wmesh_int_t 	s_t2n_n_,
				     const_wmesh_int_p 	s_t2n_v_,
				     wmesh_int_t 	s_t2n_ld_,
				     
				     wmesh_int_t 	s_q2n_m_,
				     wmesh_int_t 	s_q2n_n_,
				     const_wmesh_int_p 	s_q2n_v_,
				     wmesh_int_t 	s_q2n_ld_) 
  {

    WMESH_CHECK_POINTER(c_);

  
    static constexpr wmesh_int_t 	refpyr_c[] 	=  { 0,0,0,
							     1,0,0,
							     1,1,0,
							     0,1,0,
							     0,0,1};
  
    static constexpr  wmesh_int_t 	reftetra_c[] 	= { 0,0,0,
							    1,0,0,
							    0,1,0,
							    0,0,1};
  
    static constexpr wmesh_int_t 	refwedge_c[] 	= { 0,0,0,
							    1,0,0,
							    0,1,0,
							    0,0,1,
							    1,0,1,
							    0,1,1};
  
    static constexpr wmesh_int_t 	refhexa_c[] 	= { 0,0,0,
							    1,0,0,
							    1,1,0,
							    0,1,0,
							    0,0,1,
							    1,0,1,
							    1,1,1,
							    0,1,1};

  
#if 0
    static const wmesh_int_t wedgequadface_cnc[] = {0,2,5,3,
						    4,5,2,1,
						    0,3,4,1};
	
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
    
    static const wmesh_int_t hexaface_cnc[]  = {0,3,2,1,
						4,5,6,7,				 
						0,1,5,4,
						1,2,6,5,				 
						2,3,7,6,
						3,0,4,7 };
      

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
      
    static const wmesh_int_t s_q2n[] = {0,3,2,1};
    static const wmesh_int_t s_s_t2n[] = {0,1,4,
					1,2,4, 
					2,3,4, 
					3,0,4};
      
#endif
    
    wmesh_int_t reface[32];
    wmesh_int_t reedge[8];
    switch(element_)
      {

      case WMESH_ELEMENT_TETRAHEDRON:
	{
	  wmesh_int_t idx = 0;	
	  const wmesh_int_t n = (degree_>0)? degree_-1 : 0;
#if 0
	  static constexpr wmesh_int_t ref_icoo_n  = 4;
	  static constexpr wmesh_int_t ref_icoo_ld = 3;
	  for (wmesh_int_t i=0;i<ref_icoo_n;++i)
	    {
	      *ic = shift_ + reftetra_c[ref_icoo_ld*i+0] * degree_;
	      *jc = shift_ + reftetra_c[ref_icoo_ld*i+1] * degree_;
	      *kc = shift_ + reftetra_c[ref_icoo_ld*i+2] * degree_;
	      ic += ic_ld;
	      jc += jc_ld;
	      kc += jc_ld;
	    }
	
	  for (wmesh_int_t iedge=0;iedge<s_e2n_n_;++iedge)
	    {	    
	      reedge[0] = reftetra_c[3*s_e2n_v_[iedge * 2 + 0]+0];
	      reedge[1] = reftetra_c[3*s_e2n_v_[iedge * 2 + 0]+1];
	      reedge[2] = reftetra_c[3*s_e2n_v_[iedge * 2 + 0]+2];
	    
	      reedge[3] = reftetra_c[3*s_e2n_v_[iedge * 2 + 1]+0];
	      reedge[4] = reftetra_c[3*s_e2n_v_[iedge * 2 + 1]+1];
	      reedge[5] = reftetra_c[3*s_e2n_v_[iedge * 2 + 1]+2];
	      for (wmesh_int_t i=0;i<n;++i)
		{
		  const wmesh_int_t l1 = (i+1);
		  const wmesh_int_t l0 = degree_-(i+1);
		  *ic = shift_ + reedge[0] * l0 + reedge[3] * l1;
		  *jc = shift_ + reedge[1] * l0 + reedge[4] * l1;
		  *kc = shift_ + reedge[2] * l0 + reedge[5] * l1;			    
		  ic += ic_ld;
		  jc += jc_ld;
		  kc += jc_ld;
		} 
	    } 
	
	  if (degree_>2)
	    {
	      for (wmesh_int_t iface=0;iface<s_t2n_n;++iface)
		{
		  reface[0] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 0]+0];
		  reface[1] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 0]+1];
		  reface[2] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 0]+2];
		
		  reface[3] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 1]+0];
		  reface[4] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 1]+1];
		  reface[5] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 1]+2];
		
		  reface[6] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 2]+0];
		  reface[7] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 2]+1];
		  reface[8] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 2]+2];
		
		  for (wmesh_int_t i=0;i<n;++i)
		    {
		      for (wmesh_int_t j=0;j<n-(i+1);++j)
			{
			  const wmesh_int_t l0 		= degree_-( (i+1) + (j+1) );
			  const wmesh_int_t l1 		= (i+1);
			  const wmesh_int_t l2 		= (j+1);			
			  *ic = reface[0] * l0 + reface[3] * l1 + reface[6] * l2;
			  *jc = reface[1] * l0 + reface[4] * l1 + reface[7] * l2;
			  *kc = reface[2] * l0 + reface[5] * l1 + reface[8] * l2;			
			  ic += ic_ld;
			  jc += jc_ld;
			  kc += jc_ld;
			} 
		    }	      	
		} 
	    }
#endif	
	  c_[idx*c_ld_+0] = 0;	  
	  c_[idx*c_ld_+1] = 0;	
	  c_[idx*c_ld_+2] = 0;	
	  ++idx;
	  c_[idx*c_ld_+0] = degree_;	  
	  c_[idx*c_ld_+1] = 0;	
	  c_[idx*c_ld_+2] = 0;	
	  ++idx;
	  c_[idx*c_ld_+0] = 0;	  
	  c_[idx*c_ld_+1] = degree_;	
	  c_[idx*c_ld_+2] = 0;	
	  ++idx;
	  c_[idx*c_ld_+0] = 0;	  
	  c_[idx*c_ld_+1] = 0;
	  c_[idx*c_ld_+2] = degree_;	
	  ++idx;
	
	  for (wmesh_int_t iedge=0;iedge<s_e2n_n_;++iedge)
	    {
	      reedge[0] = reftetra_c[3*s_e2n_v_[iedge * 2 + 0]+0];
	      reedge[1] = reftetra_c[3*s_e2n_v_[iedge * 2 + 0]+1];
	      reedge[2] = reftetra_c[3*s_e2n_v_[iedge * 2 + 0]+2];
	      reedge[3] = reftetra_c[3*s_e2n_v_[iedge * 2 + 1]+0];
	      reedge[4] = reftetra_c[3*s_e2n_v_[iedge * 2 + 1]+1];
	      reedge[5] = reftetra_c[3*s_e2n_v_[iedge * 2 + 1]+2];
	      for (wmesh_int_t i=0;i<n;++i)
		{
		  const wmesh_int_t l1 = (i+1);
		  const wmesh_int_t l0 = degree_-(i+1);		  
		  c_[idx*c_ld_+0] = reedge[0] * l0 + reedge[3] * l1;
		  c_[idx*c_ld_+1] = reedge[1] * l0 + reedge[4] * l1;
		  c_[idx*c_ld_+2] = reedge[2] * l0 + reedge[5] * l1;			    
		  ++idx;
		} 
	    } 
      
	  if (degree_>2)
	    {
	      for (wmesh_int_t iface=0;iface<s_t2n_n_;++iface)
		{
		  reface[0] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 0]+0];
		  reface[1] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 0]+1];
		  reface[2] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 0]+2];

		  reface[3] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 1]+0];
		  reface[4] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 1]+1];
		  reface[5] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 1]+2];

		  reface[6] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 2]+0];
		  reface[7] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 2]+1];
		  reface[8] = reftetra_c[3*s_t2n_v_[iface * s_t2n_ld_ + 2]+2];
		  
		  for (wmesh_int_t i=0;i<n;++i)
		    {
		      for (wmesh_int_t j=0;j<n-(i+1);++j)
			{
			  const wmesh_int_t l0 		= degree_-( (i+1) + (j+1) );
			  const wmesh_int_t l1 		= (i+1);
			  const wmesh_int_t l2 		= (j+1);			
			  c_[idx*c_ld_+0] = reface[0] * l0 + reface[3] * l1 + reface[6] * l2;
			  c_[idx*c_ld_+1] = reface[1] * l0 + reface[4] * l1 + reface[7] * l2;
			  c_[idx*c_ld_+2] = reface[2] * l0 + reface[5] * l1 + reface[8] * l2;			
			  ++idx;
			} 
		    }
	      		  
		} 
	    }
	
	  /* INSIDE TETRA */
	  for (wmesh_int_t  k=0;k<n-1;++k)
	    {
	      for (wmesh_int_t  i=0;i<n-(k+1);++i)
		{		    
		  for (wmesh_int_t j=0;j<n-(i+1)-(k+1);j++)
		    {			  
		      c_[idx*c_ld_+0] = i+1;
		      c_[idx*c_ld_+1] = j+1;	
		      c_[idx*c_ld_+2] = k+1;	
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
	  c_[c_ld_*idx+0] = 0;
	  c_[c_ld_*idx+1] = 0;
	  c_[c_ld_*idx+2] = 0;
	  ++idx;
	  c_[c_ld_*idx+0] = degree_;
	  c_[c_ld_*idx+1] = 0;
	  c_[c_ld_*idx+2] = 0;
	  ++idx;
	  c_[c_ld_*idx+0] = degree_;
	  c_[c_ld_*idx+1] = degree_;
	  c_[c_ld_*idx+2] = 0;
	  ++idx;
	  c_[c_ld_*idx+0] = 0;
	  c_[c_ld_*idx+1] = degree_;
	  c_[c_ld_*idx+2] = 0;
	  ++idx;
	  c_[c_ld_*idx+0] = 0;
	  c_[c_ld_*idx+1] = 0;
	  c_[c_ld_*idx+2] = degree_;
	  ++idx;
      
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = (i+1);	  
	      c_[c_ld_*idx+1] = 0;	
	      c_[c_ld_*idx+2] = 0;	
	      ++idx;	      
	    }
      
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = degree_; 
	      c_[c_ld_*idx+1] = (i+1);	
	      c_[c_ld_*idx+2] = 0;	
	      ++idx;	      
	    } 
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = degree_-(i+1); 
	      c_[c_ld_*idx+1] = degree_;	
	      c_[c_ld_*idx+2] = 0;	
	      ++idx;	      
	    }
      
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = 0; 
	      c_[c_ld_*idx+1] = degree_-(i+1);	
	      c_[c_ld_*idx+2] = 0;	
	      ++idx;	      
	    }
      
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = 0; 
	      c_[c_ld_*idx+1] = 0;	
	      c_[c_ld_*idx+2] = (i+1);	
	      ++idx;	      
	    }
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = degree_-(i+1); 
	      c_[c_ld_*idx+1] = 0;	
	      c_[c_ld_*idx+2] = (i+1);	
	      ++idx;	      
	    }

	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = degree_-(i+1);
	      c_[c_ld_*idx+1] = degree_-(i+1);	
	      c_[c_ld_*idx+2] = (i+1);	
	      ++idx;	      
	    }

	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = 0;
	      c_[c_ld_*idx+1] = degree_-(i+1);	
	      c_[c_ld_*idx+2] = (i+1);	
	      ++idx;	      
	    }

	  //      std::cout << "idx " << idx << std::endl;           

	  for (wmesh_int_t iface=0;iface<s_t2n_n_;++iface)
	    {
	    
	      reface[0] = refpyr_c[3*s_t2n_v_[iface * s_t2n_ld_ + 0]+0];
	      reface[1] = refpyr_c[3*s_t2n_v_[iface * s_t2n_ld_ + 0]+1];
	      reface[2] = refpyr_c[3*s_t2n_v_[iface * s_t2n_ld_ + 0]+2];
		  
	      reface[3] = refpyr_c[3*s_t2n_v_[iface * s_t2n_ld_ + 1]+0];
	      reface[4] = refpyr_c[3*s_t2n_v_[iface * s_t2n_ld_ + 1]+1];
	      reface[5] = refpyr_c[3*s_t2n_v_[iface * s_t2n_ld_ + 1]+2];
		  
	      reface[6] = refpyr_c[3*s_t2n_v_[iface * s_t2n_ld_ + 2]+0];
	      reface[7] = refpyr_c[3*s_t2n_v_[iface * s_t2n_ld_ + 2]+1];
	      reface[8] = refpyr_c[3*s_t2n_v_[iface * s_t2n_ld_ + 2]+2];

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
		    
		      c_[c_ld_*idx+0] = xi;	  
		      c_[c_ld_*idx+1] = xj;	
		      c_[c_ld_*idx+2] = xk;
#if 0
		      std::cout << "----" << std::endl;
		      std::cout << "s " << xi << std::endl;
		      std::cout << "s " << xj << std::endl;
		      std::cout << "s " << xk << std::endl;
#endif
		      ++idx;
		    } 
		} 
	    } 


	  for (wmesh_int_t iface=0;iface<s_q2n_n_;++iface)
	    {
	    
	      reface[0] = refpyr_c[3*s_q2n_v_[iface * s_q2n_ld_ + 0]+0];
	      reface[1] = refpyr_c[3*s_q2n_v_[iface * s_q2n_ld_ + 0]+1];
	      reface[2] = refpyr_c[3*s_q2n_v_[iface * s_q2n_ld_ + 0]+2];
		  
	      reface[3] = refpyr_c[3*s_q2n_v_[iface * s_q2n_ld_ + 1]+0];
	      reface[4] = refpyr_c[3*s_q2n_v_[iface * s_q2n_ld_ + 1]+1];
	      reface[5] = refpyr_c[3*s_q2n_v_[iface * s_q2n_ld_ + 1]+2];
		  
	      reface[6] = refpyr_c[3*s_q2n_v_[iface * s_q2n_ld_ + 2]+0];
	      reface[7] = refpyr_c[3*s_q2n_v_[iface * s_q2n_ld_ + 2]+1];
	      reface[8] = refpyr_c[3*s_q2n_v_[iface * s_q2n_ld_ + 2]+2];
		
	      reface[9]  = refpyr_c[3*s_q2n_v_[iface * s_q2n_ld_ + 3]+0];
	      reface[10] = refpyr_c[3*s_q2n_v_[iface * s_q2n_ld_ + 3]+1];
	      reface[11] = refpyr_c[3*s_q2n_v_[iface * s_q2n_ld_ + 3]+2];

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
		      c_[c_ld_*idx+0] = xi / degree_;	  
		      c_[c_ld_*idx+1] = xj / degree_;	
		      c_[c_ld_*idx+2] = xk / degree_;

#if 0
		      std::cout << "idx " << idx << std::endl;
		      std::cout << c_[c_ld_*idx+0] << std::endl;
		      std::cout << c_[c_ld_*idx+1] << std::endl;
		      std::cout << c_[c_ld_*idx+2] << std::endl;
#endif
		      ++idx;
		    } 
		} 

	    }
      
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
			  c_[c_ld_*idx+0] = i;	  
			  c_[c_ld_*idx+1] = j;
			  c_[c_ld_*idx+2] = k;
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

	  c_[c_ld_*idx+0] = 0;
	  c_[c_ld_*idx+1] = 0;
	  c_[c_ld_*idx+2] = 0;
	  ++idx;
	  c_[c_ld_*idx+0] = degree_;
	  c_[c_ld_*idx+1] = 0;
	  c_[c_ld_*idx+2] = 0;
	  ++idx;
	  c_[c_ld_*idx+0] = 0;
	  c_[c_ld_*idx+1] = degree_;
	  c_[c_ld_*idx+2] = 0;
	  ++idx;

	  c_[c_ld_*idx+0] = 0;
	  c_[c_ld_*idx+1] = 0;
	  c_[c_ld_*idx+2] = degree_;
	  ++idx;
	  c_[c_ld_*idx+0] = degree_;
	  c_[c_ld_*idx+1] = 0;
	  c_[c_ld_*idx+2] = degree_;
	  ++idx;
	  c_[c_ld_*idx+0] = 0;
	  c_[c_ld_*idx+1] = degree_;
	  c_[c_ld_*idx+2] = degree_;
	  ++idx;

	  /* EDGE 0 */
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = (i+1);	  
	      c_[c_ld_*idx+1] = 0;	
	      c_[c_ld_*idx+2] = 0;	
	      ++idx;	      
	    } 
	  /* EDGE 1 */
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = degree_-(i+1);	  
	      c_[c_ld_*idx+1] = (i+1);	
	      c_[c_ld_*idx+2] = 0;	
	      ++idx;	      
	    } 
	  /* EDGE 2 */
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = 0;	  
	      c_[c_ld_*idx+1] = degree_-(i+1);
	      c_[c_ld_*idx+2] = 0;	
	      ++idx;	      
	    } 

	  /* EDGE 3 */
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = (i+1);	  
	      c_[c_ld_*idx+1] = 0;	
	      c_[c_ld_*idx+2] = degree_;	
	      ++idx;	      
	    } 
	  /* EDGE 4 */
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = degree_-(i+1);	  
	      c_[c_ld_*idx+1] = (i+1);	
	      c_[c_ld_*idx+2] = degree_;	
	      ++idx;	      
	    } 
	  /* EDGE 5 */
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = 0;	  
	      c_[c_ld_*idx+1] = degree_-(i+1);
	      c_[c_ld_*idx+2] = degree_;	
	      ++idx;	      
	    } 

	  /* EDGE 6 */
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = 0;	  
	      c_[c_ld_*idx+1] = 0;	
	      c_[c_ld_*idx+2] = (i+1);	
	      ++idx;	      
	    } 
	  /* EDGE 7 */
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = degree_;	  
	      c_[c_ld_*idx+1] = 0;	
	      c_[c_ld_*idx+2] = (i+1);	
	      ++idx;	      
	    } 
	  /* EDGE 8 */
	  for (wmesh_int_t i=0;i<n;++i)
	    {
	      c_[c_ld_*idx+0] = 0;	  
	      c_[c_ld_*idx+1] = degree_;
	      c_[c_ld_*idx+2] = (i+1);	
	      ++idx;	      
	    } 

	  for (wmesh_int_t iface=0;iface<s_t2n_n_;++iface)
	    {
	      for (wmesh_int_t i=0;i<n;++i)
		{
		  for (wmesh_int_t j=0;j<n-(i+1);++j)
		    {
		      c_[c_ld_*idx+0] = (i+1);	  
		      c_[c_ld_*idx+1] = (j+1);
		      c_[c_ld_*idx+2] = iface*degree_;
		      ++idx;	      
		    } 
		}       
	    } 

	  for (wmesh_int_t iface=0;iface<s_q2n_n_;++iface)
	    {
	      reface[0] = refwedge_c[3*s_q2n_v_[iface * s_q2n_ld_ + 0]+0];
	      reface[1] = refwedge_c[3*s_q2n_v_[iface * s_q2n_ld_ + 0]+1];
	      reface[2] = refwedge_c[3*s_q2n_v_[iface * s_q2n_ld_ + 0]+2];
	      
	      reface[3] = refwedge_c[3*s_q2n_v_[iface * s_q2n_ld_ + 1]+0];
	      reface[4] = refwedge_c[3*s_q2n_v_[iface * s_q2n_ld_ + 1]+1];
	      reface[5] = refwedge_c[3*s_q2n_v_[iface * s_q2n_ld_ + 1]+2];
	      
	      reface[6] = refwedge_c[3*s_q2n_v_[iface * s_q2n_ld_ + 2]+0];
	      reface[7] = refwedge_c[3*s_q2n_v_[iface * s_q2n_ld_ + 2]+1];
	      reface[8] = refwedge_c[3*s_q2n_v_[iface * s_q2n_ld_ + 2]+2];
	      
	      reface[9]  = refwedge_c[3*s_q2n_v_[iface * s_q2n_ld_ + 3]+0];
	      reface[10] = refwedge_c[3*s_q2n_v_[iface * s_q2n_ld_ + 3]+1];
	      reface[11] = refwedge_c[3*s_q2n_v_[iface * s_q2n_ld_ + 3]+2];

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
			    
		      c_[c_ld_*idx+0] 	= xi/degree_;	  
		      c_[c_ld_*idx+1] 	= xj/degree_;	
		      c_[c_ld_*idx+2] 	= xk/degree_;	
		      ++idx;
		    } 
		} 
	    } 

	  /* interior */
	  for (wmesh_int_t k=0;k<n;++k)
	    {
	      for (wmesh_int_t i=0;i<n;++i)
		{
		  for (wmesh_int_t j=0;j<n-(i+1);++j)
		    {
		      c_[c_ld_*idx+0] = (i+1);	  
		      c_[c_ld_*idx+1] = (j+1);
		      c_[c_ld_*idx+2] = (k+1);
		      ++idx;	      
		    } 
		} 
	    } 
	  return WMESH_STATUS_SUCCESS;
	}

      case WMESH_ELEMENT_HEXAHEDRON:
	{
	  wmesh_int_t idx = 0;
      
	  c_[idx*c_ld_+0] = 0;	  
	  c_[idx*c_ld_+1] = 0;	
	  c_[idx*c_ld_+2] = 0;	
	  ++idx;
      
	  c_[idx*c_ld_+0] = degree_;	  
	  c_[idx*c_ld_+1] = 0;	
	  c_[idx*c_ld_+2] = 0;	
	  ++idx;

	  c_[idx*c_ld_+0] = degree_;	  
	  c_[idx*c_ld_+1] = degree_;	
	  c_[idx*c_ld_+2] = 0;	
	  ++idx;
      
	  c_[idx*c_ld_+0] = 0;	  
	  c_[idx*c_ld_+1] = degree_;	
	  c_[idx*c_ld_+2] = 0;	
	  ++idx;
	
	  c_[idx*c_ld_+0] = 0;	  
	  c_[idx*c_ld_+1] = 0;	
	  c_[idx*c_ld_+2] = degree_;	
	  ++idx;

	  c_[idx*c_ld_+0] = degree_;	  
	  c_[idx*c_ld_+1] = 0;	
	  c_[idx*c_ld_+2] = degree_;	
	  ++idx;

	  c_[idx*c_ld_+0] = degree_;	  
	  c_[idx*c_ld_+1] = degree_;	
	  c_[idx*c_ld_+2] = degree_;	
	  ++idx;

	  c_[idx*c_ld_+0] = 0;	  
	  c_[idx*c_ld_+1] = degree_;	
	  c_[idx*c_ld_+2] = degree_;	
	  ++idx;


	  for (wmesh_int_t iedge=0;iedge<s_e2n_n_;++iedge)
	    {
	      reedge[0] = refhexa_c[3*s_e2n_v_[iedge * s_e2n_ld_ + 0]+0];
	      reedge[1] = refhexa_c[3*s_e2n_v_[iedge * s_e2n_ld_ + 0]+1];
	      reedge[2] = refhexa_c[3*s_e2n_v_[iedge * s_e2n_ld_ + 0]+2];

	      reedge[3] = refhexa_c[3*s_e2n_v_[iedge * s_e2n_ld_ + 1]+0];
	      reedge[4] = refhexa_c[3*s_e2n_v_[iedge * s_e2n_ld_ + 1]+1];
	      reedge[5] = refhexa_c[3*s_e2n_v_[iedge * s_e2n_ld_ + 1]+2];
		
	      for (wmesh_int_t i=0;i<degree_-1;++i)
		{
		  const wmesh_int_t l1 		= (i+1);
		  const wmesh_int_t l0 		= degree_-(i+1);
		  c_[idx*c_ld_+0] 	= reedge[0] * l0 + reedge[3] * l1;
		  c_[idx*c_ld_+1] 	= reedge[1] * l0 + reedge[4] * l1;
		  c_[idx*c_ld_+2] 	= reedge[2] * l0 + reedge[5] * l1;
		  ++idx;
		} 
	    }
	


	  for (wmesh_int_t iface=0;iface<s_q2n_n_;++iface)
	    {
	      reface[0] = refhexa_c[3*s_q2n_v_[iface * s_q2n_ld_ + 0]+0];
	      reface[1] = refhexa_c[3*s_q2n_v_[iface * s_q2n_ld_ + 0]+1];
	      reface[2] = refhexa_c[3*s_q2n_v_[iface * s_q2n_ld_ + 0]+2];

	      reface[3] = refhexa_c[3*s_q2n_v_[iface * s_q2n_ld_ + 1]+0];
	      reface[4] = refhexa_c[3*s_q2n_v_[iface * s_q2n_ld_ + 1]+1];
	      reface[5] = refhexa_c[3*s_q2n_v_[iface * s_q2n_ld_ + 1]+2];

	      reface[6] = refhexa_c[3*s_q2n_v_[iface * s_q2n_ld_ + 2]+0];
	      reface[7] = refhexa_c[3*s_q2n_v_[iface * s_q2n_ld_ + 2]+1];
	      reface[8] = refhexa_c[3*s_q2n_v_[iface * s_q2n_ld_ + 2]+2];

	      reface[9]  = refhexa_c[3*s_q2n_v_[iface * s_q2n_ld_ + 3]+0];
	      reface[10] = refhexa_c[3*s_q2n_v_[iface * s_q2n_ld_ + 3]+1];
	      reface[11] = refhexa_c[3*s_q2n_v_[iface * s_q2n_ld_ + 3]+2];

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
			
		      c_[idx*c_ld_+0] 	= xi/degree_;	  
		      c_[idx*c_ld_+1] 	= xj/degree_;	
		      c_[idx*c_ld_+2] 	= xk/degree_;
		      ++idx;
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
		      c_[idx*c_ld_+0] = i+1;
		      c_[idx*c_ld_+1] = j+1;	
		      c_[idx*c_ld_+2] = k+1;
		      ++idx;
		    }
		}
	    }
	  return WMESH_STATUS_SUCCESS;
	}
      
      case WMESH_ELEMENT_EDGE:
      case WMESH_ELEMENT_TRIANGLE:
      case WMESH_ELEMENT_QUADRILATERAL:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
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
  
    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);

    return WMESH_STATUS_INVALID_ARGUMENT;

  };	  

  wmesh_status_t bms_ordering(wmesh_int_t		element_,
			      wmesh_int_t		degree_,
			      wmesh_int_t		c_storage_,
			      wmesh_int_t		c_m_,
			      wmesh_int_t		c_n_,
			      wmesh_int_p		c_v_,
			      wmesh_int_t		c_ld_) 
  {

    
    wmesh_status_t status;
    WMESH_CHECK_POINTER(c_v_);
    switch(element_)
      {      
      case WMESH_ELEMENT_TRIANGLE:
      case WMESH_ELEMENT_QUADRILATERAL:
	{
	  status = bms_ordering_face(element_,
				     degree_,
				     c_storage_,
				     c_m_,
				     c_n_,
				     c_v_,
				     c_ld_);
	  WMESH_STATUS_CHECK(status);
	  return WMESH_STATUS_SUCCESS;
	}

      case WMESH_ELEMENT_EDGE:
	{
	  status = bms_ordering_edge(degree_,
				     c_storage_,
				     c_m_,
				     c_n_,
				     c_v_,
				     c_ld_);
	  WMESH_STATUS_CHECK(status);
	  return WMESH_STATUS_SUCCESS;
	}
	
      case WMESH_ELEMENT_TETRAHEDRON:
      case WMESH_ELEMENT_PYRAMID:
      case WMESH_ELEMENT_WEDGE:
      case WMESH_ELEMENT_HEXAHEDRON:
	{
	  wmesh_int_t s_e2n_size;
	  wmesh_int_t s_e2n_ptr[4+1];
	  wmesh_int_t s_e2n_m[4];
	  wmesh_int_t s_e2n_n[4];
	  wmesh_int_t s_e2n_v[128];
	  wmesh_int_t s_e2n_ld[4];

	  
	  status = bms_s_e2n(3,
			     &s_e2n_size,
			     s_e2n_ptr,
			     s_e2n_m,
			     s_e2n_n,
			     s_e2n_v,
			     s_e2n_ld);
	  WMESH_STATUS_CHECK(status);

	  wmesh_int_t s_t2n_size;
	  wmesh_int_t s_t2n_ptr[4+1];
	  wmesh_int_t s_t2n_m[4];
	  wmesh_int_t s_t2n_n[4];
	  wmesh_int_t s_t2n_v[128];
	  wmesh_int_t s_t2n_ld[4];
	  status = bms_s_t2n(&s_t2n_size,
			     s_t2n_ptr,
			     s_t2n_m,
			     s_t2n_n,
			     s_t2n_v,
			     s_t2n_ld);
	  WMESH_STATUS_CHECK(status);
	  
	  wmesh_int_t s_q2n_size;
	  wmesh_int_t s_q2n_ptr[4+1];
	  wmesh_int_t s_q2n_m[4];
	  wmesh_int_t s_q2n_n[4];
	  wmesh_int_t s_q2n_v[128];
	  wmesh_int_t s_q2n_ld[4];
	  
	  status = bms_s_q2n(&s_q2n_size,
			     s_q2n_ptr,
			     s_q2n_m,
			     s_q2n_n,
			     s_q2n_v,
			     s_q2n_ld);
	  WMESH_STATUS_CHECK(status);
	  
	  status = bms_ordering_volume(element_,
				       degree_,
				       c_storage_,
				       c_m_,
				       c_n_,
				       c_v_,
				       c_ld_,

				       s_e2n_m[element_-4],
				       s_e2n_n[element_-4],
				       s_e2n_v + s_e2n_ptr[element_-4],
				       s_e2n_ld[element_-4],

				       s_t2n_m[element_-4],
				       s_t2n_n[element_-4],
				       s_t2n_v + s_t2n_ptr[element_-4],
				       s_t2n_ld[element_-4],

				       s_q2n_m[element_-4],
				       s_q2n_n[element_-4],
				       s_q2n_v + s_q2n_ptr[element_-4],
				       s_q2n_ld[element_-4]);
	  
	  return WMESH_STATUS_SUCCESS;
	}
	
      case WMESH_ELEMENT_NODE:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
	  break;
	}

      }  
    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
  };	  

#if 0  
  wmesh_status_t bms_ordering_sub(wmesh_int_t		element_,
				  wmesh_int_t 		sub_element_,
				  wmesh_int_t 		sub_idx_,
				  wmesh_int_t		degree_,
				  wmesh_int_t		c_storage_,
				  wmesh_int_t		c_m_,
				  wmesh_int_t		c_n_,
				  const_wmesh_int_p	c_v_,
				  wmesh_int_t		c_ld_,
				  wmesh_int_t		sub_c_storage_,
				  wmesh_int_t		sub_c_m_,
				  wmesh_int_t		sub_c_n_,
				  const_wmesh_int_p	sub_c_v_,
				  wmesh_int_t		sub_c_ld_,
				  wmesh_int_t 		s_x2n_m_,
				  wmesh_int_t 		s_x2n_n_,
				  const_wmesh_int_p	s_x2n_v_,
				  wmesh_int_t 		s_x2n_ld_,
				  wmesh_int_p		sub_c2n_,
				  wmesh_int_t		sub_c2n_inc_) 
  {

    wmesh_int_t ref[32];
    {
      wmesh_status_t status = bms_ordering_vertices(element_,
						    ref);
      WMESH_STATUS_CHECK(status);
    }

    wmesh_int_t topodim 	= c_m_;
    wmesh_int_t sub_topodim 	= sub_c_m_;
    wmesh_int_t s_cnc[4];
 	  
    wmesh_int_t bez[3];
    wmesh_int_t sub_bez[3];
    
    //
    // Get local connectivity.
    // 
    for (wmesh_int_t i=0;i<s_x2n_m_;++i)
      {
	s_cnc[i] = s_x2n_v_[s_x2n_ld_*sub_idx+i];
      }
    
    wmesh_int_t sub_nodes[12];
    for (wmesh_int_t i=0;i<s_x2n_m_;++i)
      {
	for (wmesh_int_t k=0;k<topodim;++k)
	  {
	    sub_nodes[i * topodim + k] = ref[s_cnc[i] * topodim + k];
	  }
      }
    
    //
    // Get the reference coordinates.
    //
	  	  
    // ref[topodim * + i];
    //
    // Traverse. (1,0) * ( degree_ -  r) + r*  (0,1) 
    //
    wmesh_int_t lagr[4];
    for (wmesh_int_t j=0;j<sub_c_n_;++j)
      {
	for (wmesh_int_t i=0;i<sub_topodim;++i)
	  {
	    sub_bez[i] = sub_c_v_[sub_c_ld_*j+i];
	  }
	
	if (s_x2n_m_ == 2)
	  {
	    lagr[0] = (degree_ - sub_bez[0]);
	    lagr[1] = sub_bez[0];
	  }
	else if (s_x2n_m_ == 3)
	  {
	    lagr[0] = (degree_ - sub_bez[0]-sub_bez[1]);
	    lagr[1] = sub_bez[0];
	    lagr[2] = sub_bez[1];
	  }
	else if (s_x2n_m_ == 4)
	  {
	    lagr[0] = (degree_ - sub_bez[0])*(degree_ - sub_bez[1]);
	    lagr[1] = sub_bez[0] * (degree_ - sub_bez[1]);
	    lagr[2] = sub_bez[0] * sub_bez[1];
	    lagr[3] = (degree_ - sub_bez[0]) * sub_bez[1];
	  }
	
	//
	// Interpolate.
	//
	for (wmesh_int_t i=0;i<topodim;++i)
	  {
	    bez[i] = 0;
	    for (wmesh_int_t k = 0;k<s_x2n_m_;++k)
	      {
		bez[i] += sub_nodes[topodim * k + i] * lagr[k];
	      }
	  }
      }	      

    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
  };	  
#endif

  
  wmesh_status_t bms_ordering_topoid(wmesh_int_t		element_,
				     wmesh_int_t		degree_,
				     wmesh_int_t		topoid_n_,
				     wmesh_int_p		topoid_v_,
				     wmesh_int_t		topoid_inc_)
  {
    

    WMESH_CHECK_POINTER(topoid_v_);

    switch(element_)
      {      
      
#define TREAT_CASE(_c)							\
	case _c:							\
	  {								\
	    wmesh_int_t topodim;					\
	    wmesh_status_t status = bms_template_element2topodim<_c>(&topodim); \
	    WMESH_STATUS_CHECK(status);					\
									\
	    const wmesh_int_t num_nodes 		= bms_traits_element<_c>::s_num_nodes; \
	    const wmesh_int_t num_edges 		= bms_traits_element<_c>::s_num_edges; \
	    const wmesh_int_t num_triangles 		= bms_traits_element<_c>::s_num_triangles; \
	    const wmesh_int_t num_quadrilaterals 	= bms_traits_element<_c>::s_num_quadrilaterals; \
									\
	    const wmesh_int_t num_dofs_0  = (degree_ > 0) ? num_nodes : 0; \
	    const wmesh_int_t num_dofs_1  = num_edges * bms_template_ndofs_interior<WMESH_ELEMENT_EDGE>(degree_); \
	    const wmesh_int_t num_dofs_2t = (topodim > 1) ? num_triangles * bms_template_ndofs_interior<WMESH_ELEMENT_TRIANGLE>(degree_) : 0; \
	    const wmesh_int_t num_dofs_2q = (topodim > 1) ? num_quadrilaterals  * bms_template_ndofs_interior<WMESH_ELEMENT_QUADRILATERAL>(degree_) : 0; \
	    const wmesh_int_t num_dofs_3  = (topodim > 2) ? bms_template_ndofs_interior<_c>(degree_) : 0; \
									\
	    return  bms_ordering_topoid_calculate(num_dofs_0,		\
						  num_dofs_1,		\
						  num_dofs_2t,		\
						  num_dofs_2q,		\
						  num_dofs_3,		\
						  topoid_n_,		\
						  topoid_v_,		\
						  topoid_inc_);		\
	  }
	

	TREAT_CASE(WMESH_ELEMENT_EDGE);
	TREAT_CASE(WMESH_ELEMENT_TRIANGLE);
	TREAT_CASE(WMESH_ELEMENT_QUADRILATERAL);

	TREAT_CASE(WMESH_ELEMENT_TETRAHEDRON);
	TREAT_CASE(WMESH_ELEMENT_PYRAMID);
	TREAT_CASE(WMESH_ELEMENT_WEDGE);
	TREAT_CASE(WMESH_ELEMENT_HEXAHEDRON);


#undef TREAT_CASE
      case WMESH_ELEMENT_NODE:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
	  break;
	}

      }  
    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
  };	  

}






  template<typename T>
  wmesh_status_t bms_ordering_linear_shape(wmesh_int_t 			element_,
					   wmesh_int_t 			c_storage_,
					   wmesh_int_t 			c_m_,
					   wmesh_int_t 			c_n_,
					   const T * __restrict__ 	c_v_,
					   wmesh_int_t 			c_ld_,
					   wmesh_int_t 			ev_storage_,
					   wmesh_int_t 			ev_m_,
					   wmesh_int_t 			ev_n_,
					   T * __restrict__ 	ev_v_,
					   wmesh_int_t 			ev_ld_)
  {
    static constexpr T one =  static_cast<T>(1);
    static constexpr T zero =  static_cast<T>(0);
    const wmesh_int_t num_nodes = (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_n_ : c_m_;
    const wmesh_int_t r_inc = (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_ld_ : 1;
    const wmesh_int_t s_inc = (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_ld_ : 1;
    const wmesh_int_t t_inc = (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_ld_ : 1;

    const T * __restrict__ r = c_v_;
    const T * __restrict__ s = (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_v_+1 : c_v_+c_ld_*1;
    const T * __restrict__ t = (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_v_+2 : c_v_+c_ld_*2;


    switch(element_)
      {
      case WMESH_ELEMENT_EDGE:
	{
	  const wmesh_int_t l_inc = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_ld_ : 1;
	  T * __restrict__ l0 = ev_v_;
	  T * __restrict__ l1 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+1 : ev_v_+ev_ld_*1;

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l0[l_inc*i] = one - r[r_inc*i];		
	    }
	  
	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {
	      l1[l_inc*i] = r[r_inc*i];    		
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}

      case WMESH_ELEMENT_TRIANGLE:
	{
	  const wmesh_int_t l_inc = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_ld_ : 1;
	  T * __restrict__ l0 = ev_v_;
	  T * __restrict__ l1 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+1 : ev_v_+ev_ld_*1;
	  T * __restrict__ l2 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+2 : ev_v_+ev_ld_*2;

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l0[l_inc*i] = one - (r[r_inc*i]+s[s_inc*i]);		
	    }
	  
	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {
	      l1[l_inc*i] = r[r_inc*i];    		
	    }

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {
	      l2[l_inc*i] = s[s_inc*i];    		
	    }

	  return WMESH_STATUS_SUCCESS;
	}

      case WMESH_ELEMENT_QUADRILATERAL:
	{
	  const wmesh_int_t l_inc = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_ld_ : 1;
	  T * __restrict__ l0 = ev_v_;
	  T * __restrict__ l1 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+1 : ev_v_+ev_ld_*1;
	  T * __restrict__ l2 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+2 : ev_v_+ev_ld_*2;
	  T * __restrict__ l3 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+3 : ev_v_+ev_ld_*3;

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l0[l_inc*i] = (one - r[r_inc*i])*(one - s[s_inc*i]);	
	    }
	  
	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {
	      l1[l_inc*i] = r[r_inc*i]*(one - s[s_inc*i]);	
	    }

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {
	      l2[l_inc*i] = r[r_inc*i]*s[s_inc*i];	
	    }
	  
	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l3[l_inc*i] = (one - r[r_inc*i])*s[s_inc*i];	
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}

      case WMESH_ELEMENT_TETRAHEDRON:
	{
	  const wmesh_int_t l_inc = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_ld_ : 1;
	  T * __restrict__ l0 = ev_v_;
	  T * __restrict__ l1 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+1 : ev_v_+ev_ld_*1;
	  T * __restrict__ l2 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+2 : ev_v_+ev_ld_*2;
	  T * __restrict__ l3 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+3 : ev_v_+ev_ld_*3;

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l0[l_inc*i] = one - (r[r_inc*i] + s[s_inc*i] + t[t_inc*i] );	
	    }
	  
	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {
	      l1[l_inc*i] = r[r_inc*i];	
	    }

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {
	      l2[l_inc*i] = s[s_inc*i];	
	    }
	  
	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l3[l_inc*i] = t[t_inc*i];	
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}
	  
      case WMESH_ELEMENT_PYRAMID:
	{
	  const wmesh_int_t l_inc = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_ld_ : 1;
	  T * __restrict__ l0 = ev_v_;
	  T * __restrict__ l1 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+1 : ev_v_+ev_ld_*1;
	  T * __restrict__ l2 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+2 : ev_v_+ev_ld_*2;
	  T * __restrict__ l3 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+3 : ev_v_+ev_ld_*3;
	  T * __restrict__ l4 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+4 : ev_v_+ev_ld_*4;

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l0[l_inc*i] = (t[t_inc*i] < one) ?  (one - (t[t_inc*i]+r[r_inc*i]) ) * (one - (t[t_inc*i] + s[s_inc*i]) ) / (one - t[t_inc*i]) : zero;
	    }
	  
	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {
	      l1[l_inc*i] = (t[t_inc*i] < one) ?  (r[r_inc*i] * (one-t[t_inc*i]-s[s_inc*i])) / (one - t[t_inc*i]) : zero;
	    }

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {
	      l2[l_inc*i] = (t[t_inc*i] < one) ?  (r[r_inc*i] * s[s_inc*i]) / (one - t[t_inc*i]) : zero;
	    }

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {
	      l3[l_inc*i] = (t[t_inc*i] < one) ?  ( (one - (t[t_inc*i] + r[r_inc*i]) ) * s[s_inc*i] ) / (one - t[t_inc*i]) : zero;
	    }

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l4[l_inc*i] = t[t_inc*i];	
	    }
	  
	  return WMESH_STATUS_SUCCESS;
	}

	case WMESH_ELEMENT_WEDGE:
	{
	  const wmesh_int_t l_inc = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_ld_ : 1;
	  T * __restrict__ l0 = ev_v_;
	  T * __restrict__ l1 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+1 : ev_v_+ev_ld_*1;
	  T * __restrict__ l2 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+2 : ev_v_+ev_ld_*2;
	  T * __restrict__ l3 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+3 : ev_v_+ev_ld_*3;
	  T * __restrict__ l4 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+4 : ev_v_+ev_ld_*4;
	  T * __restrict__ l5 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+5 : ev_v_+ev_ld_*5;
	  
	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l0[l_inc*i] = ( one - (r[r_inc*i] + s[s_inc*i]) ) * ( one - t[t_inc*i] );
	    }
	  
	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {
	      l1[l_inc*i] = r[r_inc*i] * ( one - t[t_inc*i] );
	    }

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {
	      l2[l_inc*i] = s[s_inc*i] * ( one - t[t_inc*i] );	
	    }
	  
	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l3[l_inc*i] = ( one - (r[r_inc*i] + s[s_inc*i])) * t[t_inc*i] ;
	    }

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l4[l_inc*i] = r[r_inc*i] * t[t_inc*i] ;
	    }

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l5[l_inc*i] = s[s_inc*i] *  t[t_inc*i];
	    }

	  return WMESH_STATUS_SUCCESS;
	}

	case WMESH_ELEMENT_HEXAHEDRON:
	{
	  const wmesh_int_t l_inc = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_ld_ : 1;
	  T * __restrict__ l0 = ev_v_;
	  T * __restrict__ l1 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+1 : ev_v_+ev_ld_*1;
	  T * __restrict__ l2 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+2 : ev_v_+ev_ld_*2;
	  T * __restrict__ l3 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+3 : ev_v_+ev_ld_*3;
	  T * __restrict__ l4 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+4 : ev_v_+ev_ld_*4;
	  T * __restrict__ l5 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+5 : ev_v_+ev_ld_*5;
	  T * __restrict__ l6 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+6 : ev_v_+ev_ld_*6;
	  T * __restrict__ l7 = (ev_storage_ == WMESH_STORAGE_INTERLEAVE) ? ev_v_+7 : ev_v_+ev_ld_*7;
	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l0[l_inc*i] = ( one - r[r_inc*i] )* ( one - s[s_inc*i] ) * ( one - t[t_inc*i] );
	    }
	  
	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {
	      l1[l_inc*i] = ( ( r[r_inc*i] )* ( one - s[s_inc*i] ) * ( one - t[t_inc*i] ) );
	    }

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {
	      l2[l_inc*i] = ( ( r[r_inc*i] )* ( s[s_inc*i] ) * ( one - t[t_inc*i] ) );
	    }
	  
	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l3[l_inc*i] = ( one - r[r_inc*i] )* ( s[s_inc*i] ) * ( one - t[t_inc*i] );
	    }

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l4[l_inc*i] = ( one - r[r_inc*i] )* ( one - s[s_inc*i] ) * ( t[t_inc*i] );
	    }

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l5[l_inc*i] = ( ( r[r_inc*i] )* ( one - s[s_inc*i] ) * ( t[t_inc*i] ) );
	    }

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l6[l_inc*i] = ( ( r[r_inc*i] )* ( s[s_inc*i] ) * ( t[t_inc*i] ) );
	    }

	  for (wmesh_int_t i=0;i<num_nodes;++i)
	    {	      
	      l7[l_inc*i] = ( one - r[r_inc*i] )* ( s[s_inc*i] ) * ( t[t_inc*i] );
	    }

	  return WMESH_STATUS_SUCCESS;
	}

      }
    return WMESH_STATUS_INVALID_ENUM;
  }



template
wmesh_status_t bms_ordering_linear_shape<float>(wmesh_int_t 			element_,
						wmesh_int_t 			c_storage_,
						wmesh_int_t 			c_m_,
						wmesh_int_t 			c_n_,
						const float * __restrict__ 	c_v_,
						wmesh_int_t 			c_ld_,
						wmesh_int_t 			ev_storage_,
						wmesh_int_t 			ev_m_,
						wmesh_int_t 			ev_n_,
						 float * __restrict__ 	ev_v_,
						wmesh_int_t 			ev_ld_);

template
wmesh_status_t bms_ordering_linear_shape<double>(wmesh_int_t 			element_,
						 wmesh_int_t 			c_storage_,
						 wmesh_int_t 			c_m_,
						 wmesh_int_t 			c_n_,
						 const double * __restrict__ 	c_v_,
						 wmesh_int_t 			c_ld_,
						 wmesh_int_t 			ev_storage_,
						 wmesh_int_t 			ev_m_,
						 wmesh_int_t 			ev_n_,
						 double * __restrict__ 	ev_v_,
						 wmesh_int_t 			ev_ld_);
