#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"
#include "wmesh_utils.hpp"

#include <iostream>
#include "wmesh.hpp"
#include "GenericEncoding.hpp"

extern "C" wmesh_status_t wfe_ndofs(wmesh_int_t element_,
				    wmesh_int_t d_,
				    wmesh_int_p ndofs_)
{
#ifdef TREAT_CASE
#error TREAT_CASE is already defined
#else
#define TREAT_CASE(_f) case _f: ndofs_[0] = wfe_ndofs_template<_f>(d_); return WMESH_STATUS_SUCCESS
#endif
  switch(element_)
    {
    TREAT_CASE(WMESH_ELEMENT_NODE);
    TREAT_CASE(WMESH_ELEMENT_EDGE);
    TREAT_CASE(WMESH_ELEMENT_TRIANGLE);
    TREAT_CASE(WMESH_ELEMENT_QUADRILATERAL);
    TREAT_CASE(WMESH_ELEMENT_TETRAHEDRON);
    TREAT_CASE(WMESH_ELEMENT_PYRAMID);
    TREAT_CASE(WMESH_ELEMENT_WEDGE);
    TREAT_CASE(WMESH_ELEMENT_HEXAHEDRON);
    
    }
  return WMESH_STATUS_INVALID_ENUM;
#undef TREAT_CASE
}



template<typename T>
wmesh_status_t wmesh_permutate(wmesh_int_t 		x_m_,
			       wmesh_int_t 		x_n_,
			       T*__restrict__		x_v_,
			       wmesh_int_t 		x_ld_,
			       wmesh_int_t 		p_n_,
			       const_wmesh_int_p 	p_v_,
			       wmesh_int_t 		p_ld_,
			       wmesh_int_t              work_n_,
			       T*__restrict__		work_)
{
  if (work_n_ < x_m_*x_n_)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
    }
  if (x_n_ != p_n_)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
    }
  for (wmesh_int_t i=0;i<x_n_;++i)
    {
      wmesh_int_t j 	= p_v_[i * p_ld_];
      for (wmesh_int_t k=0;k<x_m_;++k)
	{
	  work_[x_m_*i+k] = x_v_[x_ld_*j+k];
	}
    }
  
  // 
  // copy back.
  //
  for (wmesh_int_t i=0;i<x_n_;++i)
    {
      for (wmesh_int_t k=0;k<x_m_;++k)
	{
	  x_v_[x_ld_*i+k] = work_[x_m_*i+k];
	}
    }
  return 0;
}

extern "C" const char * file_extension(const char * filename_)
{
  int i = -1;
  int len = 0;
  for (;;)
    {
      if (filename_[len] == '.')
	{
	  i = len;
	}
      if (filename_[len++] == '\0')
	{
	  break;
	}
    }
  return (i>=0) ? &filename_[i] : nullptr;
}

template<typename T>
wmesh_status_t wmesh_raw_bounding_box(wmesh_int_t 		coo_m_,
				      wmesh_int_t 		coo_n_,
				      const T *__restrict__ 	coo_v_,
				      wmesh_int_t 		coo_ld_,
				      T *__restrict__ 		box_)
{
  double c[3],p[3],q[3];
  
  {
    double initial_value 	= *coo_v_;    
    for (wmesh_int_t i=0;i<coo_m_;++i)
      {
	p[i] = initial_value;
      }    
    for (wmesh_int_t i=0;i<coo_m_;++i)
      {
	q[i] = initial_value;
      }
  }
    
  for (wmesh_int_t j=0;j<coo_n_;++j)
    {
      for (wmesh_int_t i=0;i<coo_m_;++i)
	{
	  c[i] = coo_v_[coo_ld_*j + i];
	}	
      for (wmesh_int_t i=0;i<coo_m_;++i)
	{
	  p[i] = (c[i] < p[i]) ? c[i] : p[i];
	}
      for (wmesh_int_t i=0;i<coo_m_;++i)
	{
	  q[i] = (c[i] > q[i]) ? c[i] : q[i];
	}      
    }
    
  for (wmesh_int_t i=0;i<coo_m_;++i)
    {
      box_[i] = p[i];
    }      
  for (wmesh_int_t i=0;i<coo_m_;++i)
    {
      box_[coo_m_ + i] = q[i];
    }      

  return WMESH_STATUS_SUCCESS;
}



template
wmesh_status_t wmesh_raw_bounding_box<double>(wmesh_int_t 		 coo_m_,
					      wmesh_int_t 		 coo_n_,
					      const double *__restrict__ coo_v_,
					      wmesh_int_t 		 coo_ld_,
					      double *__restrict__ 	 box_);



extern "C"
unsigned long long int hilbert_coordinate(double*__restrict__ 	crd,
					  const double*__restrict__ box,
					  int 	itr)
{
  unsigned long long int IntCrd[3], m=1LL<<63, cod;
  int b, GeoWrd, NewWrd, BitTab[3] = {1,2,4};
  int rot[8], GeoCod[8]={0,3,7,4,1,2,6,5};  /* Z curve = {5,4,7,6,1,0,3,2} */
  int HilCod[8][8] = {{0,7,6,1,2,5,4,3}, {0,3,4,7,6,5,2,1}, {0,3,4,7,6,5,2,1}, {2,3,0,1,6,7,4,5}, \
		      {2,3,0,1,6,7,4,5}, {6,5,2,1,0,3,4,7}, {6,5,2,1,0,3,4,7}, {4,3,2,5,6,1,0,7}};
  
  /* Convert double precision coordinates to integers */
  double tmp[3];
  unsigned long long int * h = (unsigned long long int *)&tmp[0];
  
  tmp[0] = (crd[0] - box[0]) * box[0+3];
  tmp[1] = (crd[1] - box[1]) * box[1+3];
  tmp[2] = (crd[2] - box[2]) * box[2+3];

  IntCrd[0] = h[0];
  IntCrd[1] = h[1];
  IntCrd[2] = h[2];

  /* Binary hilbert renumbering loop */
  
  cod = 0;
  
  rot[0] = GeoCod[0];
  rot[1] = GeoCod[1];
  rot[2] = GeoCod[2];
  rot[3] = GeoCod[3];
  rot[4] = GeoCod[4];
  rot[5] = GeoCod[5];
  rot[6] = GeoCod[6];
  rot[7] = GeoCod[7];
  
  for(b=0;b<itr;b++)
    {
      GeoWrd = 0;
      
      if(IntCrd[0] & m) GeoWrd |= BitTab[0];	  
      IntCrd[0] = IntCrd[0]<<1;

      if(IntCrd[1] & m) GeoWrd |= BitTab[1];	  
      IntCrd[1] = IntCrd[1]<<1;

      if(IntCrd[2] & m) GeoWrd |= BitTab[2];	  
      IntCrd[2] = IntCrd[2]<<1;

      NewWrd 	= rot[ GeoWrd ];      
      cod 	= cod<<3 | NewWrd;

      rot[0] = HilCod[ NewWrd ][ rot[0] ];
      rot[1] = HilCod[ NewWrd ][ rot[1] ];
      rot[2] = HilCod[ NewWrd ][ rot[2] ];
      rot[3] = HilCod[ NewWrd ][ rot[3] ];
      rot[4] = HilCod[ NewWrd ][ rot[4] ];
      rot[5] = HilCod[ NewWrd ][ rot[5] ];
      rot[6] = HilCod[ NewWrd ][ rot[6] ];
      rot[7] = HilCod[ NewWrd ][ rot[7] ];
    }
  
  return cod;
}

