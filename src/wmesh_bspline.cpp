#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"
#include <chrono>
#include <iostream>
#include "wmesh.hpp"
#include "GenericEncoding.hpp"
#include "wmesh_utils.hpp"
#include "wmesh_bspline.h"
#include <math.h>

#include "wmesh-blas.h"

#ifdef __cplusplus
extern "C" {
#endif




  
 #define hermite_interval_dtt(_n,_r,_roff,_p,_poff)			\
  { wmesh_int_t _i;								\
    for (_i=0;_i<(_n);++_i)						\
      { const double s  = ((_p)[_i*(_poff)]+1.0)*0.5;			\
	(_r)[(_roff)*_i+0] = (s*12.0-6.0);				\
	(_r)[(_roff)*_i+1] = (6.0-s*12.0);				\
	(_r)[(_roff)*_i+2] = (s*6.0-4.0);				\
	(_r)[(_roff)*_i+3] = (s*6.0-2.0);				\
      }									\
  }

#define hermite_interval_dt(_n,_r,_roff,_p,_poff)	\
  { wmesh_int_t _i;						\
    for (_i=0;_i<(_n);++_i)				\
      {							\
	const double s  = ((_p)[_i*(_poff)]+1.0)*0.5;	\
	(_r)[(_roff)*_i+0] = (s-1.0)*s*6.0;		\
	(_r)[(_roff)*_i+1] = (1.0-s)*s*6.0;		\
	(_r)[(_roff)*_i+2] = (s*(s*3.0-4.0)+1.0);	\
	(_r)[(_roff)*_i+3] = s*(s*3.0-2.0);		\
      }							\
  }


#define hermite_interval(_n,_r,_roff,_p,_poff)				\
  { wmesh_int_t _i;								\
    for (_i=0;_i<(_n);++_i)						\
      {									\
	const double s  		= ((_p)[_i*(_poff)]+1.0)*0.5;		\
	const double s2 		= s*s;					\
	(_r)[(_roff)*_i+0] 	= (s*2.0-3.0)*s2+1.0;			\
	(_r)[(_roff)*_i+1] 	= (s*(-2.0)+3.0)*s2;			\
	(_r)[(_roff)*_i+2] 	= s*(s-1.0)*(s-1.0);			\
	(_r)[(_roff)*_i+3] 	= s2*(s-1.0);				\
      }									\
  }
 

void FiniteElementSegmentHermite_Eval(const char*transpose_,
				      const_wmesh_int_p numPoints_,
				      const double*__restrict__ t_,
				      const_wmesh_int_p toff_,
				      double * __restrict__ r_,
				      const_wmesh_int_p roff_,
				      wmesh_int_p err_)
{
  static const 	double r1_2 = ((double)0.5);
  static const 	double r1 	= ((double)1.0);
  static const 	double r2 	= ((double)2.0);
  static const 	double r3 	= ((double)3.0);

  const wmesh_int_t numPoints 	= numPoints_[0];
  const wmesh_int_t toff		= toff_[0];
  const wmesh_int_t roff		= roff_[0];

  err_[0] = 0;
  if (transpose_[0] == 'N')
    {
      { wmesh_int_t pointIndex;
	for (pointIndex = 0;pointIndex<numPoints;++pointIndex)
	  {
	    const double t 	= t_[pointIndex * toff];
	    const double s  	= (r1+t)*r1_2;
	    const double sxs	= s*s;
	    const double sm1	= s-r1;
	    r_[pointIndex * roff + 1] = (r3-r2*s)*sxs;
	    r_[pointIndex * roff + 0] = r1-r_[pointIndex * roff + 1];
	    r_[pointIndex * roff + 2] = s*sm1*sm1;
	    r_[pointIndex * roff + 3] = sxs*sm1;
	  } }
    }
  else if (transpose_[0] == 'T')
    {
      { wmesh_int_t pointIndex;
	for (pointIndex = 0;pointIndex<numPoints;++pointIndex)
	  {
	    const double t 	= t_[pointIndex * toff];
	    const double s  	= (r1+t)*r1_2;
	    const double sxs	= s*s;
	    const double sm1	= s-r1;
	    r_[pointIndex + roff * 1] = (r3-r2*s)*sxs;
	    r_[pointIndex + roff * 0] = r1-r_[pointIndex + roff * 1];
	    r_[pointIndex + roff * 2] = s*sm1*sm1;
	    r_[pointIndex + roff * 3] = sxs*sm1;
	  } }
    }
  else
    {
      err_[0] = 7;
    }
}




void FiniteElementSegmentHermite_EvalDt(const char*transpose_,
					const_wmesh_int_p numPoints_,
					const double*__restrict__ t_,
					const_wmesh_int_p toff_,
					double * __restrict__ r_,
					const_wmesh_int_p roff_,
					wmesh_int_p err_)
{
  static const 	double r1_2 = ((double)0.5);
  static const 	double r1 	= ((double)1.0);
  static const 	double r2 	= ((double)2.0);
  static const 	double r3 	= ((double)3.0);
  static const 	double r4 	= ((double)4.0);
  static const 	double r6 	= ((double)6.0);
  
  const wmesh_int_t numPoints 	= numPoints_[0];
  const wmesh_int_t toff		= toff_[0];
  const wmesh_int_t roff		= roff_[0];

  err_[0] = 0;
  if (transpose_[0] == 'N')
    {
      { wmesh_int_t pointIndex;
	for (pointIndex = 0;pointIndex<numPoints;++pointIndex)
	  {
	    const double t 	= t_[pointIndex * toff];
	    const double s  	= (r1+t)*r1_2;
	    const double sm1	= s-r1;
	    r_[pointIndex * roff + 0] = sm1*s*r6;
	    r_[pointIndex * roff + 1] = -r_[pointIndex * roff + 0];
	    r_[pointIndex * roff + 2] = (s*r3-r4)*s+r1;
	    r_[pointIndex * roff + 3] = s*(s*r3-r2);
	  } }
    }
  else if (transpose_[0] == 'T')
    {
      { wmesh_int_t pointIndex;
	for (pointIndex = 0;pointIndex<numPoints;++pointIndex)
	  {
	    const double t 	= t_[pointIndex * toff];
	    const double s  	= (r1+t)*r1_2;
	    const double sm1	= s-r1;
	    r_[pointIndex + roff * 0] = sm1*s*r6;
	    r_[pointIndex + roff * 1] = -r_[pointIndex + roff * 0];
	    r_[pointIndex + roff * 2] = (s*r3-r4)*s+r1;
	    r_[pointIndex + roff * 3] = s*(s*r3-r2);
	  } }
    }
  else
    {
      err_[0] = 7;
    }
}

void FiniteElementSegmentHermite_EvalDDt(const char*			transpose_,
					 const_wmesh_int_p  		numPoints_,
					 const double * __restrict__  	t_,
					 const_wmesh_int_p 		toff_,
					 double * __restrict__ 		r_,
					 const_wmesh_int_p 		roff_,
					 wmesh_int_p 			err_)
{
  static const 	double r1_2 	= ((double)0.5);
  static const 	double r1 	= ((double)1.0);
  static const 	double r2 	= ((double)2.0);
  static const 	double r3 	= ((double)3.0);
  static const 	double r6 	= ((double)6.0);

  const wmesh_int_t numPoints 	= numPoints_[0];
  const wmesh_int_t toff		= toff_[0];
  const wmesh_int_t roff		= roff_[0];

  err_[0] = 0;
  if (transpose_[0] == 'N')
    {

	{ wmesh_int_t pointIndex;
	for (pointIndex = 0;pointIndex<numPoints;++pointIndex)
	  {
	    const double t 	= t_[pointIndex * toff];
	    const double s  	= (r1+t)*r1_2;
	    r_[pointIndex * roff + 0] = r6*(s*r2-r1);
	    r_[pointIndex * roff + 1] = -r_[pointIndex * roff + 0];
	    r_[pointIndex * roff + 2] = r2*(s*r3-r2);
	    r_[pointIndex * roff + 3] = r2*(s*r3-r1);
	  } }
    }
  else if (transpose_[0] == 'T')
    {
      { wmesh_int_t pointIndex;
	for (pointIndex = 0;pointIndex<numPoints;++pointIndex)
	  {
	    const double t 	= t_[pointIndex * toff];
	    const double s  	= (r1+t)*r1_2;
	    r_[pointIndex + roff * 0] = r6*(s*r2-r1);
	    r_[pointIndex + roff * 1] = -r_[pointIndex + roff * 0];
	    r_[pointIndex + roff * 2] = r2*(s*r3-r2);
	    r_[pointIndex + roff * 3] = r2*(s*r3-r1);
	  } }
    }
  else
    {
      err_[0] = 7;
    }
}

  struct wmesh_bspline_t
  {
    wmesh_int_t  	m_dim;
    wmesh_int_t 	m_numPoints;
    double		m_length;  
    double*		m_ctrlpts;  
    double*		m_pts;  
    double*		m_derivatives;  
    double*		m_normalized_length;  
  };


    void wmesh_bspline_eval_ddt	(const wmesh_bspline_t* 	self_,
				 const double * 		s_,
				 wmesh_int_t iedge_,
				 double*acc_)
  {
    if (iedge_>=self_->m_numPoints-1)
      {
	printf("OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n");
	exit(1);
      }
    const wmesh_int_t  dim = self_->m_dim;
    double hermite[4];
    hermite_interval_dtt(1,hermite,1,s_,1);   
    { wmesh_int_t  idim =0;
      for (idim =0;idim<dim;++idim)
	{
	  acc_[idim]  = 
	    hermite[0]   * self_->m_pts[iedge_*dim+idim]
	    + hermite[1] * self_->m_pts[(iedge_+1)*dim+idim]
	    + hermite[2] * self_->m_derivatives[iedge_*dim+idim]
	    + hermite[3] * self_->m_derivatives[(iedge_+1)*dim+idim];
	} } 
  }

  wmesh_status_t wmesh_bspline_dcoo	(const wmesh_bspline_t*__restrict__ 	self_,
					 wmesh_int_t 				idx_,
					 double * 				coo_)
  {
    WMESH_POINTER_CHECK(self_);
    WMESH_POINTER_CHECK(coo_);
    
    wmesh_int_t dim;
    
    wmesh_status_t  status = wmesh_bspline_dimension(self_,&dim);
    WMESH_STATUS_CHECK(status);
    
    for (wmesh_int_t dimIndex = 0;dimIndex<dim;++dimIndex)
      {	
	coo_[dimIndex] = self_->m_pts[idx_*dim+dimIndex];
      } 
    return WMESH_STATUS_SUCCESS;    
  }

  wmesh_status_t wmesh_bspline_ctrl_dcoo	(const wmesh_bspline_t*__restrict__ 	self_,
						 wmesh_int_t 				idx_,
						 double * 				coo_)
  {
    WMESH_POINTER_CHECK(self_);
    WMESH_POINTER_CHECK(coo_);
    
    wmesh_int_t dim;

    wmesh_status_t  status = wmesh_bspline_dimension(self_,&dim);
    WMESH_STATUS_CHECK(status);
    
    for (wmesh_int_t dimIndex = 0;dimIndex<dim;++dimIndex)
      {	
	coo_[dimIndex] = self_->m_ctrlpts[idx_*dim+dimIndex];
      } 
    return WMESH_STATUS_SUCCESS;    
  }

  wmesh_status_t  wmesh_bspline_write_medit(const wmesh_bspline_t*	bspline_,
					    const char * 			filename_)
  {
    

    wmesh_int_t dim  = bspline_->m_dim;
    const wmesh_int_t 	numPoints 	= bspline_->m_numPoints;
    const wmesh_int_t 	numEdges 	= numPoints-1;
    
    double tangent[3],
      ctrlPointCoordinates[3],
      pointCoordinates[3];
    
    
    { FILE * out = fopen(filename_,"w");
      fprintf(out,"MeshVersionFormatted\n1\n");
      fprintf(out,"Dimension\n" WMESH_INT_FORMAT "\n",dim);
      fprintf(out,"Vertices\n " WMESH_INT_FORMAT "\n",numPoints * 2);
      
      { wmesh_int_t ctrlPointIndex;
	for (ctrlPointIndex=0;ctrlPointIndex<numPoints;++ctrlPointIndex)
	  {
	    wmesh_bspline_ctrl_dcoo(bspline_,ctrlPointIndex,ctrlPointCoordinates);
	    { wmesh_int_t idim = 0;
	      for (idim = 0;idim<dim;++idim)
		{
		  fprintf(out," %e",ctrlPointCoordinates[idim]);
		} }
	  fprintf(out," 0\n");
	} }

    { wmesh_int_t pointIndex;
      for (pointIndex=0;pointIndex<numPoints;++pointIndex)
	{
	  wmesh_bspline_dcoo(bspline_,pointIndex,pointCoordinates);
	  { wmesh_int_t idim = 0;
	    for (idim =0;idim<dim;++idim)
	      {
		fprintf(out," " "%e" "",pointCoordinates[idim]);
	      } }
	  fprintf(out," 1\n");
	} }


    fprintf(out,"Edges\n" "" WMESH_INT_FORMAT "" "\n",numEdges * 2);
    { wmesh_int_t edgeIndex;
      for (edgeIndex=0;edgeIndex<numEdges;++edgeIndex)
	{	
	  fprintf(out,"" "" WMESH_INT_FORMAT "" " " "" WMESH_INT_FORMAT "" " 0\n",edgeIndex+1,edgeIndex+2);
	} }
    { wmesh_int_t edgeIndex;
      for (edgeIndex=0;edgeIndex<numEdges;++edgeIndex)
	{	
	  fprintf(out,"" "" WMESH_INT_FORMAT "" " " "" WMESH_INT_FORMAT "" " 1\n",numPoints+edgeIndex+1,numPoints+edgeIndex+2);
	} }


    fprintf(out,"Normals\n" "" WMESH_INT_FORMAT "" "\n",numPoints);

    { wmesh_int_t pointIndex;
      for (pointIndex=0;pointIndex<numPoints;++pointIndex)
	{
	  wmesh_bspline_dtan(bspline_,pointIndex,tangent);
	  {wmesh_int_t idim = 0;
	    for (idim = 0;idim<dim;++idim)
	      {
		fprintf(out," " "%e" "",tangent[idim]);
	      } }
	  fprintf(out,"\n");	  
	} }
    
    fprintf(out,"NormalAtVertices\n" "" WMESH_INT_FORMAT "" "\n",numPoints);
    
    { wmesh_int_t pointIndex;
      for (pointIndex=0;pointIndex<numPoints;++pointIndex)
	{
	  fprintf(out,"" "" WMESH_INT_FORMAT "" " " "" WMESH_INT_FORMAT "" "\n",numPoints+pointIndex+1,pointIndex+1);	  
	} }

    fprintf(out,"End\n");
    fclose(out); }
    return WMESH_STATUS_SUCCESS;
}

  
  void wmesh_bspline_eval_dt(const wmesh_bspline_t* 	const 	self_,
			     const double * 		const	s_,
			     wmesh_int_t 		iedge_,
			     double * __restrict__ 		const	normal_)
{

  wmesh_int_t  dim = self_->m_dim;
  double hermite[4];
  hermite_interval_dt(1,hermite,1,s_,1);   

  /* calcul de la tangente */
  { wmesh_int_t  idim =0;
    for (idim =0;idim<dim;++idim)
      {
	normal_[idim]  = 
	  hermite[0]   * self_->m_pts[iedge_*dim+idim]
	  + hermite[1] * self_->m_pts[(iedge_+1)*dim+idim]
	  + hermite[2] * self_->m_derivatives[iedge_*dim+idim]
	  + hermite[3] * self_->m_derivatives[(iedge_+1)*dim+idim];
      } } 

}

  void Vect_normalized(double * __restrict__ const self_,
		       wmesh_int_t  dim_)
  {
    double norm = ((double)0.0);
    { wmesh_int_t  idim =0;
      for (idim =0;idim<dim_;++idim)
	{
	  norm += self_[idim]*self_[idim];
	} } 
    norm = sqrt(norm);
    { wmesh_int_t  idim =0;
      for (idim =0;idim<dim_;++idim)
	{
	  self_[idim]/=norm;
	} } 
  }
#if 0  
  static void Blas_dcopy(const_wmesh_int_p len_,const double * src_,const_wmesh_int_p src_ld_,double * dest_,const_wmesh_int_p dest_ld_)
  {
    for (wmesh_int_t i=0;i<len_[0];++i)
      {
	dest_[dest_ld_[0]*i] = src_[src_ld_[0]*i];
      }
  }
#endif
  

  wmesh_bspline_t* wmesh_bspline_new2(wmesh_int_t  		dim_,
				      wmesh_int_t 		numPoints_,
				      const double *    	ctrlpts_,
				      wmesh_int_t 		ctrlptsoff_)
  {

  static const 	double 	oneOnTwo 			= ((double)0.5);
  
  const 	wmesh_int_t 	dimension 		= dim_;
  const 	wmesh_int_t 	numEdges			= numPoints_-1;





  wmesh_bspline_t*self		= (wmesh_bspline_t*)calloc(1,sizeof(wmesh_bspline_t));
  self->m_dim 					= dim_;
  self->m_numPoints 				= numPoints_;
  self->m_pts 					= (double * __restrict__)malloc(sizeof(double)*numPoints_*dim_);
  self->m_derivatives				= (double * __restrict__)malloc(sizeof(double)*numPoints_*dim_);
  self->m_ctrlpts					= (double * __restrict__)malloc(sizeof(double)*numPoints_*dim_);
  self->m_normalized_length			= (double * __restrict__)malloc(sizeof(double)*numPoints_);
  self->m_normalized_length[0] 			= ((double)0.0);
  self->m_normalized_length[numEdges] 		= ((double)1.0);
  self->m_length				= ((double)0.0);
  

  /*########################################################################*/

  { wmesh_int_t  idim = 0;
    for (idim = 0;idim<dim_;++idim)
      {
    
	BLAS_dcopy(&numPoints_,(double*)&ctrlpts_[idim],&ctrlptsoff_,&self->m_ctrlpts[idim],(wmesh_int_p)&dimension);
  } }
  
  /*########################################################################*/
  { wmesh_int_t  idim = 0;
    for (idim = 0;idim<dim_;++idim)
      {
	self->m_pts[idim] = ctrlpts_[idim];
      } }
  /**/
  { wmesh_int_t i;
    for (i=1;i<numEdges;++i)
      {
	{ wmesh_int_t  idim = 0;
	  for (idim = 0;idim<dim_;++idim)
	    {
	      self->m_pts[i*dim_+idim] 		= ctrlpts_[i*ctrlptsoff_+idim];
	    } }
      } }
  
  { wmesh_int_t  idim = 0;
    for (idim = 0;idim<dim_;++idim)
      {
	self->m_pts[numEdges*dim_+idim] = ctrlpts_[numEdges*ctrlptsoff_+idim];
      } }





  { wmesh_int_t  idim = 0;
    for (idim = 0;idim<dim_;++idim)
      {
	self->m_derivatives[idim] = (ctrlpts_[1*ctrlptsoff_+idim]-ctrlpts_[idim] );
      } }

  { wmesh_int_t i;
    for (i=1;i<numEdges;++i)
      {
	{ wmesh_int_t  idim = 0;
	  for (idim = 0;idim<dim_;++idim)
	    {
	      self->m_derivatives[i*dim_+idim] 	= (ctrlpts_[(i+1)*ctrlptsoff_+idim]-ctrlpts_[(i-1)*ctrlptsoff_+idim])*oneOnTwo;
	    } }
      } }

  { wmesh_int_t  idim = 0;
    for (idim = 0;idim<dim_;++idim)
      {
	self->m_derivatives[numEdges*dim_+idim] = (ctrlpts_[(numPoints_-1)*ctrlptsoff_+idim]-ctrlpts_[(numPoints_-2)*ctrlptsoff_+idim]);
      } }


  return self;
}


  
  wmesh_status_t wmesh_bspline_dlength	(const wmesh_bspline_t*__restrict__ self_,double * __restrict__ out_length_)
  {
    WMESH_POINTER_CHECK(self_);
    WMESH_POINTER_CHECK(out_length_);
    out_length_[0] =  self_->m_length;
    return WMESH_STATUS_SUCCESS;    
  }

  wmesh_status_t wmesh_bspline_num_nodes	(const wmesh_bspline_t*__restrict__ 	self_,
						 wmesh_int_p 				out_num_nodes_)
  {
    WMESH_POINTER_CHECK(self_);
    WMESH_POINTER_CHECK(out_num_nodes_);
    out_num_nodes_[0] =  self_->m_numPoints;
    return WMESH_STATUS_SUCCESS;    
  }
  
  wmesh_status_t wmesh_bspline_num_edges	(const wmesh_bspline_t*__restrict__ 	self_,
						 wmesh_int_p 				out_num_edges_)
  {
    WMESH_POINTER_CHECK(self_);
    WMESH_POINTER_CHECK(out_num_edges_);
    out_num_edges_[0] =  self_->m_numPoints-1;
    return WMESH_STATUS_SUCCESS;    
  }

    wmesh_status_t wmesh_bspline_dimension	(const wmesh_bspline_t*__restrict__ 	self_,
						 wmesh_int_p 				out_dimension_)
  {
    WMESH_POINTER_CHECK(self_);
    WMESH_POINTER_CHECK(out_dimension_);
    out_dimension_[0] =  self_->m_dim;
    return WMESH_STATUS_SUCCESS;    
  }

  
  

  wmesh_status_t wmesh_bspline_dtan	(const wmesh_bspline_t*__restrict__ 	self_,
					 wmesh_int_t 				idx_,
					 double * 				tan_)
  {
    WMESH_POINTER_CHECK(self_);
    WMESH_POINTER_CHECK(tan_);
    
    wmesh_int_t dim;
    
    wmesh_status_t  status = wmesh_bspline_dimension(self_,&dim);
    WMESH_STATUS_CHECK(status);
    
    for (wmesh_int_t dimIndex = 0;dimIndex<dim;++dimIndex)
      {	
	tan_[dimIndex] = self_->m_derivatives[idx_*dim+dimIndex];
      } 
    return WMESH_STATUS_SUCCESS;    
  }

wmesh_status_t 	wmesh_bspline_free(wmesh_bspline_t * self_)
{
  if (self_)
    {
      if (self_->m_ctrlpts)
	{
	  free(self_->m_ctrlpts);
	  self_->m_ctrlpts = nullptr;
	}
      if (self_->m_pts)
	{
	  free(self_->m_pts);
	  self_->m_pts = nullptr;
	}
      if (self_->m_derivatives)
	{
	  free(self_->m_derivatives);
	  self_->m_derivatives = nullptr;
	}
      if (self_->m_normalized_length)
	{
	  free(self_->m_normalized_length);
	  self_->m_normalized_length = nullptr;
	}
    }
    return WMESH_STATUS_SUCCESS;
}



bool wmesh_bspline_findEdge(const wmesh_bspline_t*  self_,
			    double			t_,
			    wmesh_int_p iedge_,
			    double*	s_)
{
  wmesh_int_t num_nodes = 0;
  wmesh_bspline_num_nodes(self_,&num_nodes);
  
  { wmesh_int_t k;
    for (k=1;k<num_nodes;++k)
      {
	if (self_->m_normalized_length[k]>=t_)
	  {
	    break;	    
	  }
      }
    if (k<num_nodes)
      {
	iedge_[0] 	= k-1;
	s_[0]		= ( (t_ - self_->m_normalized_length[k-1])/(self_->m_normalized_length[k]-self_->m_normalized_length[k-1]));	
	return true;
      }
    else
      {
	return false;
      } }
}




  wmesh_status_t wmesh_bspline_deval	(const wmesh_bspline_t* 	__restrict__	self_,
					 const double * 		__restrict__	s_,
					 wmesh_int_t			iedge_,
					 double * __restrict__	position_)
  {
    wmesh_int_t dim;    
    wmesh_status_t  status = wmesh_bspline_dimension(self_,&dim);
    WMESH_STATUS_CHECK(status);
    
    double hermite[4];
    hermite_interval(1,hermite,1,s_,1);   
    
    for (wmesh_int_t idim = 0;idim<dim;++idim)
      {
	position_[idim]  = 
	  hermite[0]   * self_->m_pts[iedge_*dim+idim]
	  + hermite[1] * self_->m_pts[(iedge_+1)*dim+idim]
	  + hermite[2] * self_->m_derivatives[iedge_*dim+idim]
	  + hermite[3] * self_->m_derivatives[(iedge_+1)*dim+idim];
      } 
    return WMESH_STATUS_SUCCESS;
  }

  wmesh_status_t  wmesh_bspline_transform_coordinates(const wmesh_bspline_t*	self_,
						      wmesh_int_t		nv_,
						      double * 			coo_,
						      wmesh_int_t 		coo_ld_,
						      wmesh_int_p		ids_,
						      const_wmesh_int_p	 	offsetNode_)
  {
    wmesh_int_t dim;    

    //
    // Get dimension.
    //
    wmesh_status_t  status = wmesh_bspline_dimension(self_,&dim);
    WMESH_STATUS_CHECK(status);

    double coo0[3];
    double coo[3];
    double tangent[3];
    wmesh_int_t N = nv_;
    wmesh_int_t Nz = N / offsetNode_[0];
    if (N%offsetNode_[0]!=0)
      {
	printf("wmesh_bspline_transform_coordinates::wrong offset\n");
	exit(1);
      }
    
    double * __restrict__ approximationCourbure 	= (double*)malloc(sizeof(double)*(self_->m_numPoints-1)*3);
    double * __restrict__ maxCourbureEdge 	= (double*)malloc(sizeof(double)*(self_->m_numPoints-1));
    double * __restrict__ maxCourburePoint 	= (double*)malloc(sizeof(double)*(self_->m_numPoints));    
    double * __restrict__ basis = (double*)malloc(sizeof(double)*12*Nz);
    if (dim!=3)
      {
	printf("temporary dim invalid\n");
	exit(1);
      }

  for (wmesh_int_t iv=0;iv<offsetNode_[0];++iv)
    {
      wmesh_int_t idx 	= (Nz - 1) * offsetNode_[0] + iv;
      wmesh_int_t idTopology 	= ids_[idx];
      for (wmesh_int_t i=0;i<dim;++i)
	{
	  coo0[i] = coo_[idx*coo_ld_ + i];
	}
      
      if (coo0[2]>1.0)
	{
	  coo0[2] = ((double)1.0);
	}
      
      for (wmesh_int_t i=0;i<dim;++i)
	{
	  coo_[idx*coo_ld_ + i] = coo0[i];
	}	
    } 
  
  for (wmesh_int_t iv=0;iv<offsetNode_[0];++iv)
    {
      wmesh_int_t idx = iv;
      wmesh_int_t idTopology = ids_[idx];
      for (wmesh_int_t i=0;i<dim;++i)
	{
	  coo0[i] = coo_[idx*coo_ld_ + i];
	}
      
      if (coo0[2]<0.0)
	{
	  coo0[2] = ((double)0.0);
	}
      
      for (wmesh_int_t i=0;i<dim;++i)
	{
	  coo_[idx*coo_ld_ + i] = coo0[i];
	}      
    } 
  
  
  for (wmesh_int_t i=0;i<self_->m_numPoints-1;++i)
    {
      double nn[3];
      double s1[1];	
      s1[0] = ((double)-1.0);
      wmesh_bspline_eval_ddt(self_,s1,i,nn);
      approximationCourbure[3*i+0] = nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2];
      
      s1[0] = ((double)0.0);
      wmesh_bspline_eval_ddt(self_,s1,i,nn);
      approximationCourbure[3*i+1] = nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2];
      
      s1[0] = ((double)1.0);
      wmesh_bspline_eval_ddt(self_,s1,i,nn);
      approximationCourbure[3*i+2] = nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2];
#if 0
      printf("BBB %f %f %f\n",nn[0],nn[1],nn[2]);
      printf("AAA "ifmt" "rfmt" "rfmt" "rfmt"\n",i,approximationCourbure[3*i+0],approximationCourbure[3*i+1],approximationCourbure[3*i+2]);
#endif
    } 



  for (wmesh_int_t i=0;i<self_->m_numPoints-1;++i)
    {
      /**/
	const double a0 = approximationCourbure[3*i+0];
	const double a1 = approximationCourbure[3*i+1];
	const double a2 = approximationCourbure[3*i+2];
	double nn[3];
	double s1[1];
	
	double r = -1.0;
	const double d0 = (2.0*r-1.0)*0.5 * a0 + (-2.0 * r)*a1 + (2.0*r+1.0)*0.5*a2;
	r = 1.0;
	const double d1 = (2.0*r-1.0)*0.5 * a0 + (-2.0 * r)*a1 + (2.0*r+1.0)*0.5*a2;
	const double dd = d1*d0;
	if ( (dd<0.0) )
	  {
	    /*
	      a+b=d1
	      b-a=d0

	      a+b=d1
	      b-a=d1
	      
	      b=d1
	      a=0
	      

	      b = (d1+d0)/2.0
	      a = (d1-d0)/2.0
	      x = -b/a
	      
	      x = (d1+d0)/(d0-d1)
	    */
	    s1[0] = (d1+d0)/(d0-d1);
	  }
	else  
	  {
	    if (a0<a2)
	      {
		s1[0] = 1.0;
	      }
	    else
	      {
		s1[0] = -1.0;
	      }
	  }
	wmesh_bspline_eval_ddt(self_,s1,i,nn);
	maxCourbureEdge[i] = (nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2]);
      } 


  for (wmesh_int_t i=0;i<self_->m_numPoints-1;++i)
    {
      double nn[3];
      double s1[1];	
      
      double mx = 0.0;
      double xx=-1.0;
      while (xx<=1.0)
	{
	  s1[0] = xx;
	  wmesh_bspline_eval_ddt(self_,s1,i,nn);
	  double s3 = sqrt(nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2]);
	  if (s3>mx)
	    mx=s3;
	  
	  xx+=1.0e-2;
	}
      maxCourbureEdge[i] = mx;
#if 0
      printf("BBB %f %f %f\n",nn[0],nn[1],nn[2]);
      printf("AAA "ifmt" "rfmt" "rfmt" "rfmt"\n",i,approximationCourbure[3*i+0],approximationCourbure[3*i+1],approximationCourbure[3*i+2]);
#endif
    }


  maxCourburePoint[0] = maxCourbureEdge[0];
    for (wmesh_int_t i=1;i<self_->m_numPoints-1;++i)
      {
	maxCourburePoint[i] = std::max(maxCourbureEdge[i-1],maxCourbureEdge[i]);
      }
  maxCourburePoint[self_->m_numPoints-1] = maxCourbureEdge[self_->m_numPoints-2];


 


  double * __restrict__ ss = (double*)malloc(sizeof(double)*self_->m_numPoints*3);
  
    for (wmesh_int_t i=0;i<self_->m_numPoints;++i)
      {
	ss[3*i+0] = ((double)i);
	ss[3*i+1] = maxCourburePoint[i];
	ss[3*i+2] = ((double)0.0);
      } 
  wmesh_bspline_t* splineRayon = wmesh_bspline_new2(3,
						    self_->m_numPoints,
						    ss,
						    3);
  
  double * __restrict__ rayon = (double*)malloc(sizeof(double)*Nz);
  double minrayon=1e+30;
  { wmesh_int_t iz;
    for (iz=0;iz<Nz;++iz)
      {
	double * __restrict__ local_basis 	= &basis[12*iz];	


	wmesh_int_t idx = iz*offsetNode_[0];
	wmesh_int_t idTopology = ids_[idx];
	for (wmesh_int_t i=0;i<dim;++i)
	  {
	    coo0[i] = coo_[idx*coo_ld_ + i];
	  }

	const double s0 = coo0[ 2  ];
	double s1[1];
	s1[0]= s0;

	wmesh_int_t iedge = -1;
	bool edgeHasBeenFound = wmesh_bspline_findEdge(self_,s0,&iedge,s1);
	if (! edgeHasBeenFound)
	  {
	    printf("! edgeHasBeenFound\n",s0);
	    exit(1);
	  }
	double posrayon[3];
	s1[0] *= ((double)2.0);
	s1[0] -= ((double)1.0);

	wmesh_bspline_deval		(self_,s1,iedge,coo);
	wmesh_bspline_eval_dt		(self_,s1,iedge,tangent);
	wmesh_bspline_deval		(splineRayon,s1,iedge,posrayon);

	Vect_normalized	(tangent,dim);

	const double c0 = maxCourburePoint[iedge];
	const double c1 = maxCourburePoint[iedge+1];

	rayon[iz] = (c0+c1)*0.5 + (c1-c0)*0.5 * s1[0];
	rayon[iz] = ((double)1.0)/rayon[iz];
	rayon[iz] = ((double)1.0)/posrayon[1];

#if 0
	printf("radius %e\n",((double)1.0)/posrayon[1]);
#endif
	if (rayon[iz]>1.0)
	  {
	    rayon[iz] = ((double)1.0);
	  }
	
	double vi[3];
	double vj[3];
	double t[3];
	
	t[0] = tangent[0];
	t[1] = tangent[1];
	t[2] = tangent[2];
	
	/* on prend une tangente perturbee qu on orthonormalise */
	vi[0] = tangent[0] + ((double)random())/((double)RAND_MAX);
	vi[1] = tangent[1] + ((double)random())/((double)RAND_MAX);
	vi[2] = tangent[2] + ((double)random())/((double)RAND_MAX);
	const double alpha = vi[0] *tangent[0] +vi[1] *tangent[1] +vi[2] *tangent[2];
	vi[0] = vi[0] - tangent[0]*alpha;
	vi[1] = vi[1] - tangent[1]*alpha;
	vi[2] = vi[2] - tangent[2]*alpha;
	Vect_normalized	(vi,dim);
	/* on calcule le produit vectoriel */
	vj[0] = vi[1] * t[2] - t[1] * vi[2];
	vj[1] = vi[2] * t[0] - t[2] * vi[0];
	vj[2] = vi[0] * t[1] - t[0] * vi[1];
	
	local_basis[0] = vi[0];
	local_basis[1] = vi[1];
	local_basis[2] = vi[2];
	local_basis[3] = vj[0];
	local_basis[4] = vj[1];
	local_basis[5] = vj[2];

	local_basis[6] = coo[0];
	local_basis[7] = coo[1];
	local_basis[8] = coo[2];

	local_basis[9] = t[0];
	local_basis[10] = t[1];
	local_basis[11] = t[2];

      } }



    for (wmesh_int_t iz=1;iz<Nz;++iz)
      {

	double * __restrict__ local_basis 			= &basis[12*iz];
	double * __restrict__ previous_local_basis 	= &basis[12*(iz-1)];
	
	const double a0 = previous_local_basis[0]*local_basis[0] + previous_local_basis[1]*local_basis[1] + previous_local_basis[2]*local_basis[2];
	const double b0 = previous_local_basis[0]*local_basis[3] + previous_local_basis[1]*local_basis[4] + previous_local_basis[2]*local_basis[5];
	double projection_of_previous_vi[3];
	double projection_of_previous_vj[3];
	projection_of_previous_vi[0]  = a0 * local_basis[0] + b0 * local_basis[3];
	projection_of_previous_vi[1]  = a0 * local_basis[1] + b0 * local_basis[4];
	projection_of_previous_vi[2]  = a0 * local_basis[2] + b0 * local_basis[5];


	Vect_normalized(projection_of_previous_vi,dim);

	double t[3];
	t[0] = local_basis[9];
	t[1] = local_basis[10];
	t[2] = local_basis[11];

	projection_of_previous_vj[0] = projection_of_previous_vi[1] * t[2] - t[1] * projection_of_previous_vi[2];
	projection_of_previous_vj[1] = projection_of_previous_vi[2] * t[0] - t[2] * projection_of_previous_vi[0];
	projection_of_previous_vj[2] = projection_of_previous_vi[0] * t[1] - t[0] * projection_of_previous_vi[1];
	Vect_normalized(projection_of_previous_vj,dim);	

	/* maintenant on change la base */
	
	local_basis[0] = projection_of_previous_vi[0];
	local_basis[1] = projection_of_previous_vi[1];
	local_basis[2] = projection_of_previous_vi[2];
	local_basis[3] = projection_of_previous_vj[0];
	local_basis[4] = projection_of_previous_vj[1];
	local_basis[5] = projection_of_previous_vj[2];
      }

  for (wmesh_int_t iz=0;iz<Nz;++iz)
    {
      double * __restrict__ local_basis = &basis[12*iz];
      /* 
	 on traverse tous les noeuds de la frame 
      */	
      { double aa[3];
	double bb[3];
	aa[0] = local_basis[0];
	aa[1] = local_basis[1];
	aa[2] = local_basis[2];
	
	bb[0] = local_basis[3];
	bb[1] = local_basis[4];
	bb[2] = local_basis[5];

	const double nx = local_basis[9];
	const double ny = local_basis[10];
	const double nz = local_basis[11];	
	/* 
	   on effectue une rotation autour de la normale d'angle PI/2 de vj
	   si on tombe sur vi on change de signe
	   on veut un repere correcte ou vj est obtenu par rotation de PI	   
	*/
	const double sx = aa[0] + nx*nx * bb[0] + (nx*ny - nz)*bb[1] + (nx*nz + ny)*bb[2];
	const double sy = aa[1] + (nx*ny+nz) * bb[0] + (ny*ny)*bb[1] + (ny*nz - nx)*bb[2];
	const double sz = aa[2] + (nx*nz-ny) * bb[0] + (ny*nz + nx)*bb[1] + (nz*nz)*bb[2];
	if (sx*sx+sy*sy+sz*sz>1.0)
	  {
	    local_basis[3]*=((double)-1.0);
	    local_basis[4]*=((double)-1.0);
	    local_basis[5]*=((double)-1.0);
	  } }
      double rr = ((double)2.0)*rayon[iz]*((double)0.5);
#if 0
      printf("minrayon %e\n",minrayon);
#endif
      { wmesh_int_t iv;
	for (iv=0;iv<offsetNode_[0];++iv)
	  {


	    
	    wmesh_int_t idx = iz*offsetNode_[0] + iv;
	    wmesh_int_t idTopology = ids_[idx];
	    for (wmesh_int_t i=0;i<dim;++i)
	      {
		coo0[i] = coo_[idx*coo_ld_ + i];
	      }
	    coo0[0] = ( (coo0[0] + ((double)1.0))*((double)0.5) - ((double)0.5))*rr;
	    if (dim > 2)
	      {
		coo0[1] = ( (coo0[1] + ((double)1.0))*((double)0.5) -((double)0.5))*rr;
	      }
	    const double xx = coo0[0];
	    const double yy = coo0[1];			
	    coo0[0]    = local_basis[6] + (xx * local_basis[0] + yy * local_basis[3]);
	    coo0[1]    = local_basis[7] + (xx * local_basis[1] + yy * local_basis[4]);
	    coo0[2]    = local_basis[8] + (xx * local_basis[2] + yy * local_basis[5]);

	    
	    for (wmesh_int_t i=0;i<dim;++i)
	      {
		coo_[idx*coo_ld_ + i] = coo0[i];
	      }
  } }
  }
  return WMESH_STATUS_SUCCESS;
}



wmesh_status_t wmesh_bspline_ddef(wmesh_bspline_t**		self__,
				  wmesh_int_t  			dim_,
				  wmesh_int_t  			numPoints_,
				  const double *  __restrict__	ctrlpts_,
				  wmesh_int_t			ctrlptsoff_)
  {
    self__[0] = (wmesh_bspline_t*)calloc(1,sizeof(wmesh_bspline_t));
    wmesh_bspline_t*self_ = self__[0];
    
    static const 	double 	oneOnSix 		= ((double)1.0)/((double)6.0);
    static const 	double	twoOnThree 		= ((double)2.0)/((double)3.0);
    static const 	double 	oneOnTwo		= ((double)0.5);
    static const 	wmesh_int_t 	nbCubaturePoints 	= ((wmesh_int_t)10);
    static const 	wmesh_int_t 	nequal1			= ((wmesh_int_t)1);
    static const 	wmesh_int_t 	nequal3 		= ((wmesh_int_t)3);
    static const 	wmesh_int_t 	nequal4 		= ((wmesh_int_t)4);
    static const 	double 	requal1 		= ((double)1.0);
    static const 	double 	requal0 		= ((double)0.0);
    const 		wmesh_int_t 	dimension 		= dim_;
    const 		wmesh_int_t 	numEdges		= numPoints_-1;
    const 		wmesh_int_t 	numEdges_X_dim 		= numEdges * dim_;
    double 				coo[3];  
    double 				tmp[nbCubaturePoints];
    
    self_->m_dim 				= dim_;
    self_->m_numPoints				= numPoints_;
    self_->m_pts 				= (double * __restrict__)malloc(sizeof(double)*numPoints_*dim_);
    self_->m_derivatives			= (double * __restrict__)malloc(sizeof(double)*numPoints_*dim_);
    self_->m_ctrlpts				= (double * __restrict__)malloc(sizeof(double)*numPoints_*dim_);
    self_->m_normalized_length			= (double * __restrict__)malloc(sizeof(double)*numPoints_);
    self_->m_normalized_length[0]		= ((double)0.0);
    self_->m_normalized_length[numEdges]	= ((double)1.0);
    self_->m_length				= ((double)0.0);

    
    double* dofValues 	= (double*)malloc		(nequal4* ((numEdges+1)*dim_ )* sizeof(double));
    
    double* evalCurve 	= (double*)malloc		(nbCubaturePoints*
						 numEdges_X_dim* sizeof(double));
    
    double* evalDtHermite = (double*)malloc		(nbCubaturePoints*
							 nequal4 * sizeof(double));
    
    
    /*########################################################################*/
    
    { wmesh_int_t  idim =0;
      for (idim =0;idim<dim_;++idim)
	{
	  BLAS_dcopy(&numPoints_,(double*)&ctrlpts_[idim],&ctrlptsoff_,&self_->m_ctrlpts[idim],(wmesh_int_p)&dimension);
	} }
    
    /*########################################################################*/
    { wmesh_int_t  idim =0;
      for (idim =0;idim<dim_;++idim)
	{
	  self_->m_pts[idim] = ctrlpts_[idim];
	} }
    /**/
    { wmesh_int_t  idim =0;
      for (idim =0;idim<dim_;++idim)
	{
	  self_->m_derivatives[idim] = (ctrlpts_[1*ctrlptsoff_+idim]-ctrlpts_[idim] );
	} }
    /**/
    { wmesh_int_t i;
      for (i=1;i<numEdges;++i)
	{
	  { wmesh_int_t  idim =0;
	    for (idim =0;idim<dim_;++idim)
	      {
		self_->m_pts[i*dim_+idim] 		= oneOnSix * ( ctrlpts_[(i-1)*ctrlptsoff_+idim]  + ctrlpts_[(i+1)*ctrlptsoff_+idim] ) + twoOnThree * ctrlpts_[i*ctrlptsoff_+idim];
	      } }
	} }
    { wmesh_int_t i;
      for (i=1;i<numEdges;++i)
	{
	  { wmesh_int_t  idim =0;
	    for (idim =0;idim<dim_;++idim)
	      {
		self_->m_derivatives[i*dim_+idim] 	= (ctrlpts_[(i+1)*ctrlptsoff_+idim]-ctrlpts_[(i-1)*ctrlptsoff_+idim])*oneOnTwo;
	      } }
	} }
    
    { wmesh_int_t  idim =0;
      for (idim =0;idim<dim_;++idim)
	{
	  self_->m_pts[numEdges*dim_+idim] = ctrlpts_[numEdges*ctrlptsoff_+idim];
	} }
    
    { wmesh_int_t  idim =0;
      for (idim =0;idim<dim_;++idim)
	{
	  self_->m_derivatives[numEdges*dim_+idim] = (ctrlpts_[(numPoints_-1)*ctrlptsoff_+idim]-ctrlpts_[(numPoints_-2)*ctrlptsoff_+idim]);
	} }
    /*########################################################################*/
    
    static double w1d[]={
      0.295524224714752870173892994651338329421046717026853601354308029755995938217152329270356595793754216722717164401252558386818490789552005826001936342494186966609562718648884168043231305061535867409083051270663865287483901746874726597515954450775158914556548308329986393605934912382356670244,
      0.295524224714752870173892994651338329421046717026853601354308029755995938217152329270356595793754216722717164401252558386818490789552005826001936342494186966609562718648884168043231305061535867409083051270663865287483901746874726597515954450775158914556548308329986393605934912382356670244,
      0.2692667193099963550912269215694693528597599384608837958005632762421534323191792767642266367092527607555958114503686983086929234693811452415564658846634423711656014432259960141729044528030344411297902977067142537534806284608399276575006911686749842814086288868533208042150419508881916391898,
      0.2692667193099963550912269215694693528597599384608837958005632762421534323191792767642266367092527607555958114503686983086929234693811452415564658846634423711656014432259960141729044528030344411297902977067142537534806284608399276575006911686749842814086288868533208042150419508881916391898,
      0.2190863625159820439955349342281631924587718705226770898809565436351999106529512812426839931772021927865912168728128876347666269080669475688309211843316656677105269915322077536772652826671027878246851010208832173320064273483254756250668415885349420711613410227291565477768928313300688702802,
      0.2190863625159820439955349342281631924587718705226770898809565436351999106529512812426839931772021927865912168728128876347666269080669475688309211843316656677105269915322077536772652826671027878246851010208832173320064273483254756250668415885349420711613410227291565477768928313300688702802,
      0.1494513491505805931457763396576973324025566396694273678354772687532386547266300109459472646347319519140057525610454363382344517067454976014713716011937109528798134828865118770953566439639333773939909201690204649083815618779157522578300343427785361756927642128792412282970150172590842897331,
      0.1494513491505805931457763396576973324025566396694273678354772687532386547266300109459472646347319519140057525610454363382344517067454976014713716011937109528798134828865118770953566439639333773939909201690204649083815618779157522578300343427785361756927642128792412282970150172590842897331,
      0.066671344308688137593568809893331792857864834320158145128694881613412064084087101776785509685058877821090054714520419331487507126254403762139304987316994041634495363706400187011242315504393526242450629832718198718647480566044117862086478449236378557180717569208295026105115288152794421677,
      0.066671344308688137593568809893331792857864834320158145128694881613412064084087101776785509685058877821090054714520419331487507126254403762139304987316994041634495363706400187011242315504393526242450629832718198718647480566044117862086478449236378557180717569208295026105115288152794421677};

    static double p1d[]={
      -0.1488743389816312108848260011297199846175648594206916957079892535159036173556685213711776297994636912300311608052553388261028901818643765402316761969968090913050737827720371059070942475859422743249837177174247346216914852902942929003193466659082433838094355075996833570230005003837280634351,
      0.1488743389816312108848260011297199846175648594206916957079892535159036173556685213711776297994636912300311608052553388261028901818643765402316761969968090913050737827720371059070942475859422743249837177174247346216914852902942929003193466659082433838094355075996833570230005003837280634351,
      -0.4333953941292471907992659431657841622000718376562464965027015131437669890777035012251027579501177212236829350409989379472742247577232492051267741032822086200952319270933462032011328320387691584063411149801129823141488787443204324766414421576788807708483879452488118549797039287926964254222,
      0.4333953941292471907992659431657841622000718376562464965027015131437669890777035012251027579501177212236829350409989379472742247577232492051267741032822086200952319270933462032011328320387691584063411149801129823141488787443204324766414421576788807708483879452488118549797039287926964254222,
      -0.6794095682990244062343273651148735757692947118348094676648171889525585753950749246150785735704803794998339020473993150608367408425766300907682741718202923543197852846977409718369143712013552962837733153108679126932544954854729341324727211680274268486617121011712030227181051010718804444161,
      0.6794095682990244062343273651148735757692947118348094676648171889525585753950749246150785735704803794998339020473993150608367408425766300907682741718202923543197852846977409718369143712013552962837733153108679126932544954854729341324727211680274268486617121011712030227181051010718804444161,
      -0.8650633666889845107320966884234930485275430149653304525219597318453747551380555613567907289460457706944046310864117651686783001614934535637392729396890950011571349689893051612072435760480900979725923317923795535739290595879776956832427702236942765911483643714816923781701572597289139322313,
      0.8650633666889845107320966884234930485275430149653304525219597318453747551380555613567907289460457706944046310864117651686783001614934535637392729396890950011571349689893051612072435760480900979725923317923795535739290595879776956832427702236942765911483643714816923781701572597289139322313,
      -0.9739065285171717200779640120844520534282699466923821192312120666965952032346361596257235649562685562582330425187742112150221686014344777799205409587259942436704413695764881258799146633143510758737119877875210567067452435368713683033860909388311646653581707125686970668737259229449284383797,
      0.9739065285171717200779640120844520534282699466923821192312120666965952032346361596257235649562685562582330425187742112150221686014344777799205409587259942436704413695764881258799146633143510758737119877875210567067452435368713683033860909388311646653581707125686970668737259229449284383797};

  { wmesh_int_t ierr;
    FiniteElementSegmentHermite_EvalDt("N",
				       &nbCubaturePoints,
				       p1d,
				       &nequal1,
				       evalDtHermite,
				       &nequal4,
				       &ierr); }

  //  hermite_interval_dt(nbCubaturePoints,Matrix_at(evalDtHermite,0,0),4,p1d,1);     
  
  { wmesh_int_t  idim =0;
    for (idim =0;idim<dim_;++idim)
      {
	BLAS_dcopy((wmesh_int_p)&numEdges,&self_->m_pts[idim],(wmesh_int_p)&dimension,dofValues+0+idim*4,(wmesh_int_p)&nequal4);
      } }
  
  { wmesh_int_t  idim =0;
    for (idim =0;idim<dim_;++idim)
      {
	BLAS_dcopy((wmesh_int_p)&numEdges,&self_->m_pts[dim_+idim],(wmesh_int_p)&dimension,dofValues+1+idim*4,(wmesh_int_p)&nequal4);
      } }
  
  { wmesh_int_t  idim =0;
    for (idim =0;idim<dim_;++idim)
      {
	BLAS_dcopy((wmesh_int_p)&numEdges,&self_->m_derivatives[idim],(wmesh_int_p)&dimension,dofValues+2+idim*4,(wmesh_int_p)&nequal4);
      } }

  { wmesh_int_t  idim =0;
    for (idim =0;idim<dim_;++idim)
      {
	BLAS_dcopy((wmesh_int_p)&numEdges,&self_->m_derivatives[dim_+idim],(wmesh_int_p)&dimension,dofValues+3+idim*4,(wmesh_int_p)&nequal4);
      } }

  
  BLAS_dgemm("N",
	     "N",
	     (wmesh_int_p)&nbCubaturePoints,
	     (wmesh_int_p)&numEdges_X_dim,
	     (wmesh_int_p)&nequal4,
	     (double*)&requal1,
	     evalDtHermite,
	     (wmesh_int_p)&nbCubaturePoints,
	     dofValues,
	     (wmesh_int_p)&nequal4,
	     (double*)&requal0,
	     evalCurve,
	     (wmesh_int_p)&nbCubaturePoints);
#if 0
  Matrix_MM		(evalDtHermite,
			 &requal1,
			 dofValues,
			 &requal0,
			 evalCurve);
#endif
  
  double length = ((double)0.0);

  { wmesh_int_t iedge;
    for (iedge=0;iedge<numEdges;++iedge)
      {
	{ wmesh_int_t i;
	  for (i=0;i<nbCubaturePoints;++i)
	    {
	      { wmesh_int_t  idim =0;
		for (idim =0;idim<dim_;++idim)
		  {
	            const double *  c = evalCurve + (idim *nbCubaturePoints  + i) + iedge*nbCubaturePoints;
		    coo[idim] = c[0];
		  } } 


	      { double s = ((double)0.0);
		{ wmesh_int_t  idim =0;
		  for (idim =0;idim<dim_;++idim)
		    {
		      s+=coo[idim]*coo[idim];
		    } }
		tmp[i] = sqrt(s);
	      }
	    } 
	  const double iedge_length 			= BLAS_ddot((wmesh_int_p)&nbCubaturePoints,w1d,(wmesh_int_p)&nequal1,tmp,(wmesh_int_p)&nequal1);
	  length 				+= iedge_length*((double)0.5);
	  self_->m_normalized_length[iedge+1] 	= length;
	}
      } }

  free(dofValues);
  free(evalCurve);

  std::cout << "length " << length << std::endl;//printf("length = "rfmt"\n",length);
  self_->m_length = length;
  { wmesh_int_t iedge;
    for (iedge=1;iedge<numPoints_;++iedge)
      {
	self_->m_normalized_length[iedge] /= self_->m_length;
      } }

#if 1
  /* 
     POUdouble UNwmesh_int_tFOdoubleMwmesh_int_tSEdouble 
  */
  { wmesh_int_t iedge;
    for (iedge=1;iedge<numPoints_;++iedge)
      {
	self_->m_normalized_length[iedge] = ((double)iedge)/((double)numPoints_-1);
      } }
#endif
  
#if 1

  /*
    -1  -1/2 0   0   0   0
    -1  -1   1   0   0   0
    0  -1  -1   1   0   0
    0   0  -1  -1   1   0
    0   0   0  -1  -1   1 
    0   0   0   0  1/2  1 


    P_{i-1} * B0(1) + P_{i} * B1(1) + T_{i-1} * B2(1) + T_{i} * B3(1)
    = 
    P_{i} * B0(0) + P_{i+1} * B1(0) + T_{i} * B2(0) + T_{i+1} * B3(0)
    
    T_{i-1} * B2(1) + T_{i} * (B3(1) -  B2(0)) - T_{i+1} * B3(0)
    = 
    - P_{i-1} * B0(1)  + P_{i} * ( B0(0) - B1(1) ) + P_{i+1} * B1(0) 

    T_{0} * B2(0) + T_{1} * B3(0) = - P_{0} * B0(0) - P_{1} * B1(0) 
    T_{n-1} * B2(1) + T_{n} * B3(1) = - P_{n-1} * B0(1) - P_{n} * B1(1) 
  */

  
  double ee[1];
  double coeff0[4];
  double coeff1[4];
  ee[0] = -1.0;


  { wmesh_int_t ierr;
    FiniteElementSegmentHermite_EvalDDt("N",
					&nequal1,
					ee,
					&nequal1,
					coeff0,
					&nequal4,
					&ierr); }
  

  ee[0] = 1.0;

  { wmesh_int_t ierr;
    FiniteElementSegmentHermite_EvalDDt("N",
					&nequal1,
					ee,
					&nequal1,
					coeff1,
					&nequal4,
					&ierr); }

  printf("coeff0 %e %e %e %e\n",coeff0[0],coeff0[1],coeff0[2],coeff0[3]);
  printf("coeff1 %e %e %e %e\n",coeff1[0],coeff1[1],coeff1[2],coeff1[3]);
  
  double * __restrict__ A 		= (double*)calloc(numPoints_*numPoints_,sizeof(double));
  double * __restrict__ B 		= (double*)calloc(numPoints_*3,sizeof(double));
  
  /* T_{0} * B2(0) + T_{1} * B3(0) = - P_{0} * B0(0) - P_{1} * B1(0) */
  A[0] 			= coeff0[2];
  A[numPoints_] 	= coeff0[3];
  B[0] 			= -self_->m_pts[0] * coeff0[0] - self_->m_pts[dim_+0] * coeff0[1];
  B[numPoints_+0] 	= -self_->m_pts[1] * coeff0[0] - self_->m_pts[dim_+1] * coeff0[1];
  B[2*numPoints_+0]	= -self_->m_pts[2] * coeff0[0] - self_->m_pts[dim_+2] * coeff0[1];
  
  { wmesh_int_t i;
    for (i=1;i<numPoints_-1;++i)
      {
	/*T_{i-1} * B2(1) + T_{i} * (B3(1) -  B2(0)) - T_{i+1} * B3(0)*/
	A[(i-1)*numPoints_+i] 	= coeff1[2];
	A[i*numPoints_+i] 	= coeff1[3]-coeff0[2];
	A[(i+1)*numPoints_+i] 	= -coeff0[3];
	/* - P_{i-1} * B0(1)  + P_{i} * ( B0(0) - B1(1) ) + P_{i+1} * B1(0) */
	B[i] 		= -self_->m_pts[ (i-1)*dim_ + 0 ] * coeff1[0] + self_->m_pts[ (i+1)*dim_ + 0 ] * coeff0[1] + self_->m_pts[ i*dim_ + 0 ] * (coeff0[0]-coeff1[1]);
	B[numPoints_+i] 	= -self_->m_pts[ (i-1)*dim_ + 1 ] * coeff1[0] + self_->m_pts[ (i+1)*dim_ + 1 ] * coeff0[1] + self_->m_pts[ i*dim_ + 1 ] * (coeff0[0]-coeff1[1]);
	B[2*numPoints_+i]	= -self_->m_pts[ (i-1)*dim_ + 2 ] * coeff1[0] + self_->m_pts[ (i+1)*dim_ + 2 ] * coeff0[1] + self_->m_pts[ i*dim_ + 2 ] * (coeff0[0]-coeff1[1]);
      } }

  /*T_{n-1} * B2(1) + T_{n} * B3(1) = - P_{n-1} * B0(1) - P_{n} * B1(1) */  
  A[(numPoints_-2)*numPoints_+numPoints_-1] 	= coeff1[2];
  A[(numPoints_-1)*numPoints_+numPoints_-1] 	= coeff1[3];  

  B[numPoints_-1] 	= -self_->m_pts[(numPoints_-2)*dim_+0] * coeff1[0] - self_->m_pts[(numPoints_-1)*dim_+0] * coeff1[1];
  B[numPoints_+numPoints_-1] 	= -self_->m_pts[(numPoints_-2)*dim_+1] * coeff1[0] - self_->m_pts[(numPoints_-1)*dim_+1] * coeff1[1];
  B[2*numPoints_+numPoints_-1]	= -self_->m_pts[(numPoints_-2)*dim_+2] * coeff1[0] - self_->m_pts[(numPoints_-1)*dim_+2] * coeff1[1];


  { wmesh_int_t i;
    for (i=0;i<numPoints_;++i)
      {
	{ wmesh_int_t j;
	  for (j=0;j<numPoints_;++j)
	    {
	      printf(" %1.1f",A[j*numPoints_+i]);
	    } }
	printf("\n");
      } }


  wmesh_int_t info_lapack  = 0;
  wmesh_int_p perm = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*numPoints_);
  { wmesh_int_t i;
    for (i=0;i<numPoints_;++i)
      {
	printf("%e %e %e\n",B[i],B[numPoints_+i],B[2*numPoints_+i]);
      } }
  
  { wmesh_int_t i;
    for (i=0;i<numPoints_-1;++i)
      {
	double c1;
	double c2;
	c1 = coeff0[0]*self_->m_pts[i*dim_+0]
	  +coeff0[1]*self_->m_pts[(i+1)*dim_+0]
	  +coeff0[2]*self_->m_derivatives[(i)*dim_+0]
	  +coeff0[3]*self_->m_derivatives[(i+1)*dim_+0];
	c2 = coeff1[0]*self_->m_pts[i*dim_+0]
	  +coeff1[1]*self_->m_pts[(i+1)*dim_+0]
	  +coeff1[2]*self_->m_derivatives[(i)*dim_+0]
	  +coeff1[3]*self_->m_derivatives[(i+1)*dim_+0];
	printf("AVANT %e %e\n",c1,c2);
      } }

  LAPACK_dgesv((wmesh_int_p)&numPoints_,
	       (wmesh_int_p)&nequal3,
	       A,
	       (wmesh_int_p)&numPoints_,
	       perm,
	       B,
	       (wmesh_int_p)&numPoints_,
	       (wmesh_int_p)&info_lapack);   
  
  //  printf("info lapack "ifmt"\n",info_lapack);

  { wmesh_int_t i;
    for (i=0;i<numPoints_;++i)
      {
	self_->m_derivatives[i*dim_+0] = B[i];
	self_->m_derivatives[i*dim_+1] = B[numPoints_+i];
	self_->m_derivatives[i*dim_+2] = B[2*numPoints_+i];
	printf("%e %e %e\n",B[i],B[numPoints_+i],B[2*numPoints_+i]);
      } }

  { wmesh_int_t i;
    for (i=0;i<numPoints_-1;++i)
      {
	double c1;
	double c2;
	c1 = coeff0[0]*self_->m_pts[i*dim_+0]
	  +coeff0[1]*self_->m_pts[(i+1)*dim_+0]
	  +coeff0[2]*self_->m_derivatives[(i)*dim_+0]
	  +coeff0[3]*self_->m_derivatives[(i+1)*dim_+0];
	c2 = coeff1[0]*self_->m_pts[i*dim_+0]
	  +coeff1[1]*self_->m_pts[(i+1)*dim_+0]
	  +coeff1[2]*self_->m_derivatives[(i)*dim_+0]
	  +coeff1[3]*self_->m_derivatives[(i+1)*dim_+0];
	printf("APdoubleES %e %e\n",c1,c2);
      } }

  free(perm);
  free(A);
  free(B);
#endif  
  return WMESH_STATUS_SUCCESS;
}



#ifdef __cplusplus
}
#endif












    
#if 0

#include "wmesh_bspline_t.h"
#include "Blas.h"
#include <math.h>
#include "Matrix.h"
#include "FiniteElementSegmentHermite.h"
#include "SparseSymbolic.h"
extern long int random();


   
  void wmesh_bspline_get_dofGeometry		(const wmesh_bspline_t* 	const 	self_,
						 cst_pwmesh_int_t 		const	edgeIndex_,
						 double * __restrict__ 			const	dofGeometry_,
						 cst_pI			const 	dofGeometryOff_)
  {
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(edgeIndex_);
  DebugVerif(0 < edgeIndex_[0]+1);
  DebugVerif(wmesh_bspline_get_numEdges(self_)>edgeIndex_[0]);
  DebugVerif(dofGeometryOff_ >= 4);
#endif
  const I 	dofGeometryOff 	= dofGeometryOff_[0];
  const wmesh_int_t  	dim 		= wmesh_bspline_get_spatialDimension	(self_);
  const I 	edgeIndex	= edgeIndex_[0];
  { wmesh_int_t  idim =0;
    for (idim =0;idim<dim;++idim)
      {
	dofGeometry_[dofGeometryOff * idim + 0] = self_->m_pts[edgeIndex*dim+idim];
	dofGeometry_[dofGeometryOff * idim + 1] = self_->m_pts[(edgeIndex+1)*dim+idim];
	dofGeometry_[dofGeometryOff * idim + 2] = self_->m_derivatives[edgeIndex*dim+idim];
	dofGeometry_[dofGeometryOff * idim + 3] = self_->m_derivatives[(edgeIndex+1)*dim+idim];
      } }
}






  

  



void wmesh_bspline_t_defInterpolatory(pwmesh_bspline_t 	const	self_,
			      cst_wmesh_int_t  		dim_,
			      const I 		numPoints_,
			      const double *  	const 	ctrlpts_,
			      const I 		ctrlptsoff_)
{
#ifndef NDEBUG
  DebugVerif(numPoints_>1);
  DebugVerif(ctrlpts_);
  DebugVerif(ctrlptsoff_>=dim_);
#endif
  static const 	I 	nequal1			= ((I)1);
  static const 	I 	nequal3 		= ((I)3);
  const 	I 	dimension 		= dim_;
  const 	I 	numEdges		= numPoints_-1;


  wmesh_bspline_t_clear(self_);

  self_->m_dim 					= dim_;
  self_->m_numPoints 				= numPoints_;
  self_->m_pts 					= (double * __restrict__)malloc(sizeof(double)*numPoints_*dim_);
  self_->m_derivatives				= (double * __restrict__)malloc(sizeof(double)*numPoints_*dim_);
  self_->m_ctrlpts				= (double * __restrict__)malloc(sizeof(double)*numPoints_*dim_);
  self_->m_normalized_length			= (double * __restrict__)malloc(sizeof(double)*numPoints_);
  self_->m_normalized_length[0] 			= ((double)0.0);
  self_->m_normalized_length[numEdges] 		= ((double)1.0);
  self_->m_length				= ((double)0.0);

  /*########################################################################*/

  { wmesh_int_t  idim =0;
    for (idim =0;idim<dim_;++idim)
      {
	BLAS_dcopy((wmesh_int_p)&numPoints_,&ctrlpts_[idim],&ctrlptsoff_,&self_->m_ctrlpts[idim],(wmesh_int_p)&dimension);
      } }
  
  { I arrayLength = numPoints_ * dim_;
    BLAS_dcopy(&arrayLength,self_->m_ctrlpts,&nequal1,self_->m_pts,&nequal1); } 
  
  /*########################################################################*/
  double * __restrict__ A = calloc(numPoints_*numPoints_,sizeof(R));
  A[0] = ((R)2.0);
  A[numPoints_] = ((R)1.0);
  { I pointIndex = 0;
    const I numPoints_minus1 = numPoints_-1;
    for (pointIndex=1;pointIndex < numPoints_minus1;++pointIndex)
      {
	A[(pointIndex-1) * numPoints_ + pointIndex] = ((R)1.0);
	A[pointIndex * numPoints_ + pointIndex] = ((R)4.0);
	A[(pointIndex+1) * numPoints_ + pointIndex] = ((R)1.0);
      } }
  A[(numPoints_-2) * numPoints_ + numPoints_-1] = ((R)1.0);
  A[(numPoints_-1) * numPoints_ + numPoints_-1] = ((R)2.0);

  pI perm = malloc(sizeof(I)*numPoints_);
  double * __restrict__ B 	  = calloc(numPoints_*3,sizeof(R));

  { I pointIndex = 0;
    const I numPoints_minus1 = numPoints_-1;
    I coordinateIndex;
    for (coordinateIndex=0;coordinateIndex<dim_;++coordinateIndex)
      {
	pointIndex = 0;
	B[numPoints_*coordinateIndex + pointIndex] = 3.0*( self_->m_pts[dim_*(pointIndex+1) + coordinateIndex] - self_->m_pts[dim_*pointIndex + coordinateIndex]  );

	for (pointIndex=1;pointIndex < numPoints_minus1;++pointIndex)
	  {
	    B[numPoints_*coordinateIndex + pointIndex] = 3.0*( self_->m_pts[dim_*(pointIndex+1) + coordinateIndex] - self_->m_pts[dim_*(pointIndex-1) + coordinateIndex]  );
	  }
	
	pointIndex = numPoints_-1;
	B[numPoints_*coordinateIndex + pointIndex] = 3.0*( self_->m_pts[dim_*(pointIndex) + coordinateIndex] - self_->m_pts[dim_*(pointIndex-1) + coordinateIndex]  );
	
      } }

  I info_lapack;
  dgesv_((wmesh_int_p)&numPoints_,
	&nequal3,
	A,
	&numPoints_,
	perm,
	B,
	&numPoints_,
	&info_lapack);   

  {
    
    I coordinateIndex,pointIndex;
    for (pointIndex=0;pointIndex < numPoints_;++pointIndex)
      {
	for (coordinateIndex=0;coordinateIndex<dim_;++coordinateIndex)
	  {
	    self_->m_derivatives[dim_*pointIndex + coordinateIndex] = B[numPoints_*coordinateIndex + pointIndex];
	  }
      }
  }
  
}




pwmesh_bspline_t wmesh_bspline_t_new(cst_wmesh_int_t  		dim_,
		     const I 		numPoints_,
		     const double *  	const 	ctrlpts_,
		     const I 		ctrlptsoff_)
{
  pwmesh_bspline_t self					= malloc(sizeof(wmesh_bspline_t));

  wmesh_bspline_t_defInterpolatory(self,dim_,numPoints_,ctrlpts_,ctrlptsoff_);
//  wmesh_bspline_t_def(self,dim_,numPoints_,ctrlpts_,ctrlptsoff_);
  return self;
}












void Vect_normalized(double * __restrict__ const self_,
		     cst_wmesh_int_t  dim_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  R norm = ((R)0.0);
  { wmesh_int_t  idim =0;
    for (idim =0;idim<dim_;++idim)
      {
	norm += self_[idim]*self_[idim];
      } } 
  norm = nsSQRT(norm);
  { wmesh_int_t  idim =0;
    for (idim =0;idim<dim_;++idim)
      {
	self_[idim]/=norm;
      } } 
}


void  wmesh_bspline_transform_coordinates(const wmesh_bspline_t* 	const	self_,
						    pMeshNode 		const 	meshNode_,
						    cst_pI		const 	offsetNode_)
{
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(meshNode_);
  DebugVerif(offsetNode_);
#endif

  cst_wmesh_int_t  dim = wmesh_bspline_get_spatialDimension(self_);
  R coo0[3];
  R coo[3];
  R tangent[3];
  const I N = MeshNode_get_nv(meshNode_);

  /* 
     calcul de la normale pour la premiere frame 
   */

  /* 
     calcul de la normale pour la frame suivante
  */
  /* 
     on ajuste la normale pour minimiser l'angle avec la normale precedente
  */
  const I Nz = N / offsetNode_[0];
  if (N%offsetNode_[0]!=0)
    {
      printf("wmesh_bspline_transform_coordinates::wrong offset\n");
      exit(1);
    }


  double * __restrict__ basis = malloc(sizeof(R)*12*Nz);
  if (dim!=__wmesh_int_t _3)
    {
      printf("temporary dim invalid\n");
      exit(1);
    }


  { I iv;
    for (iv=0;iv<offsetNode_[0];++iv)
      {
	I idTopology;
	MeshNode_get_vertex(meshNode_,
			    (Nz-1)*offsetNode_[0] + iv,
			    &idTopology,
			    coo0);
	if (coo0[2]>1.0)
	  {
	    coo0[2] = ((R)1.0);
	  }
	MeshNode_set_vertex(meshNode_,
			    (Nz-1)*offsetNode_[0] + iv,
			    idTopology,
			    coo0);
      } }


  { I iv;
    for (iv=0;iv<offsetNode_[0];++iv)
      {
	I idTopology;
	MeshNode_get_vertex(meshNode_,
			    iv,
			    &idTopology,
			    coo0);
	if (coo0[2]<0.0)
	  {
	    coo0[2] = ((R)0.0);
	  }
	MeshNode_set_vertex(meshNode_,
			    iv,
			    idTopology,
			    coo0);	
      } }



  double * __restrict__ approximationCourbure 	= malloc(sizeof(R)*(self_->m_numPoints-1)*3);
  double * __restrict__ maxCourbureEdge 		= malloc(sizeof(R)*(self_->m_numPoints-1));
  double * __restrict__ maxCourburePoint 		= malloc(sizeof(R)*(self_->m_numPoints));

  { I i;
    for (i=0;i<self_->m_numPoints-1;++i)
      {
	R nn[3];
	R s1[1];	
	s1[0] = ((R)-1.0);
	wmesh_bspline_eval_ddt(self_,s1,i,nn);
	approximationCourbure[3*i+0] = nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2];

	s1[0] = ((R)0.0);
	wmesh_bspline_eval_ddt(self_,s1,i,nn);
	approximationCourbure[3*i+1] = nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2];

	s1[0] = ((R)1.0);
	wmesh_bspline_eval_ddt(self_,s1,i,nn);
	approximationCourbure[3*i+2] = nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2];
#if 0
	printf("BBB %f %f %f\n",nn[0],nn[1],nn[2]);
	printf("AAA "ifmt" "rfmt" "rfmt" "rfmt"\n",i,approximationCourbure[3*i+0],approximationCourbure[3*i+1],approximationCourbure[3*i+2]);
#endif
      } }


  { I i;
    for (i=0;i<self_->m_numPoints-1;++i)
      {
	/**/
	const R a0 = approximationCourbure[3*i+0];
	const R a1 = approximationCourbure[3*i+1];
	const R a2 = approximationCourbure[3*i+2];
	R nn[3];
	R s1[1];
#if 0
	*l0++ = half*r*(r-r1);
	*l1++ = r1 - r*r;
	*l2++ = half*r*(r+r1);
#endif
	
	R r = -1.0;
	const R d0 = (2.0*r-1.0)*0.5 * a0 + (-2.0 * r)*a1 + (2.0*r+1.0)*0.5*a2;
	r = 1.0;
	const R d1 = (2.0*r-1.0)*0.5 * a0 + (-2.0 * r)*a1 + (2.0*r+1.0)*0.5*a2;
	const R dd = d1*d0;
	if ( (dd<0.0) )
	  {
	    /*
	      a+b=d1
	      b-a=d0

	      a+b=d1
	      b-a=d1
	      
	      b=d1
	      a=0
	      

	      b = (d1+d0)/2.0
	      a = (d1-d0)/2.0
	      x = -b/a
	      
	      x = (d1+d0)/(d0-d1)
	    */
	    s1[0] = (d1+d0)/(d0-d1);
	  }
	else  
	  {
	    if (a0<a2)
	      {
		s1[0] = 1.0;
	      }
	    else
	      {
		s1[0] = -1.0;
	      }
	  }
	wmesh_bspline_eval_ddt(self_,s1,i,nn);
	maxCourbureEdge[i] = (nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2]);
#if 0
	const R tt = MAX(maxCourbure[i],1.0e-4);
	maxCourbure[i] = ((R)1.0)/tt;
#endif
	printf("BB "ifmt" "rfmt" \n",i,maxCourbureEdge[i]);
      } }

  { I i;
    for (i=0;i<self_->m_numPoints-1;++i)
      {
	R nn[3];
	R s1[1];	
	
	R mx = 0.0;
	R xx=-1.0;
	while (xx<=1.0)
	  {
	    s1[0] = xx;
	    wmesh_bspline_eval_ddt(self_,s1,i,nn);
	    R s3 = nsSQRT(nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2]);
	    if (s3>mx)
	      mx=s3;
	    
	    xx+=1.0e-2;
	  }
	maxCourbureEdge[i] = mx;
#if 0
	printf("BBB %f %f %f\n",nn[0],nn[1],nn[2]);
	printf("AAA "ifmt" "rfmt" "rfmt" "rfmt"\n",i,approximationCourbure[3*i+0],approximationCourbure[3*i+1],approximationCourbure[3*i+2]);
#endif
      } }


  maxCourburePoint[0] = maxCourbureEdge[0];
  { I i;
    for (i=1;i<self_->m_numPoints-1;++i)
      {
	maxCourburePoint[i] = MAX(maxCourbureEdge[i-1],maxCourbureEdge[i]);
      } }
  maxCourburePoint[self_->m_numPoints-1] = maxCourbureEdge[self_->m_numPoints-2];


  { I i;
    for (i=0;i<self_->m_numPoints;++i)
      {
#if 0
	maxCourburePoint[i] = MAX(maxCourburePoint[i],1.0e-4);
#endif
	printf("AAA "ifmt" "rfmt" \n",i,maxCourburePoint[i]);
#if 0
	printf("BBB %f %f %f\n",nn[0],nn[1],nn[2]);
#endif
#if 0

	maxCourburePoint[i] = ((R)1.0)/maxCourburePoint[i];
#endif
      } }


  double * __restrict__ ss = malloc(sizeof(R)*self_->m_numPoints*3);
  { I i;
    for (i=0;i<self_->m_numPoints;++i)
      {
	ss[3*i+0] = ((R)i);
	ss[3*i+1] = maxCourburePoint[i];
	ss[3*i+2] = ((R)0.0);
      } }
  pwmesh_bspline_t splineRayon = wmesh_bspline_t_new2(__wmesh_int_t _3,
				      self_->m_numPoints,
				      ss,
				      3);

  double * __restrict__ rayon = malloc(sizeof(R)*Nz);
  R minrayon=1e+30;
  { I iz;
    for (iz=0;iz<Nz;++iz)
      {
	double * __restrict__ local_basis 	= &basis[12*iz];	
	I idTopology 	= 0;
	MeshNode_get_vertex(meshNode_,
			    iz*offsetNode_[0],
			    &idTopology,
			    coo0);
	const R s0 = coo0[ 2  ];
	R s1[1];
	s1[0]= s0;

	
	I iedge = -1;
	bool edgeHasBeenFound = wmesh_bspline_findEdge(self_,s0,&iedge,s1);
	if (! edgeHasBeenFound)
	  {
	    printf("! edgeHasBeenFound "rfmt"\n",s0);
	    exit(1);
	  }
	R posrayon[3];
	s1[0] *= ((R)2.0);
	s1[0] -= ((R)1.0);

	wmesh_bspline_eval		(self_,s1,iedge,coo);
	wmesh_bspline_eval_dt		(self_,s1,iedge,tangent);
	wmesh_bspline_eval		(splineRayon,s1,iedge,posrayon);


#if 0
	printf("rayon courbure %e\n",1.0/(MAX(nsSQRT(nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2]),0.001)));
#endif

	Vect_normalized	(tangent,dim);
	const R c0 = maxCourburePoint[iedge];
	const R c1 = maxCourburePoint[iedge+1];

	rayon[iz] = (c0+c1)*0.5 + (c1-c0)*0.5 * s1[0];
	rayon[iz] = ((R)1.0)/rayon[iz];
	rayon[iz] = ((R)1.0)/posrayon[1];


	//	printf("radius %e\n",((R)1.0)/posrayon[1]);
	if (rayon[iz]>1.0)
	  {
	    rayon[iz] = ((R)1.0);
	  }

	
	R vi[3];
	R vj[3];
	R t[3];
	
	t[0] = tangent[0];
	t[1] = tangent[1];
	t[2] = tangent[2];
	
	/* on prend une tangente perturbee qu on orthonormalise */
	vi[0] = tangent[0] + ((R)random())/((R)RAND_MAX);
	vi[1] = tangent[1] + ((R)random())/((R)RAND_MAX);
	vi[2] = tangent[2] + ((R)random())/((R)RAND_MAX);
	const R alpha = vi[0] *tangent[0] +vi[1] *tangent[1] +vi[2] *tangent[2];
	vi[0] = vi[0] - tangent[0]*alpha;
	vi[1] = vi[1] - tangent[1]*alpha;
	vi[2] = vi[2] - tangent[2]*alpha;
	Vect_normalized	(vi,dim);
	/* on calcule le produit vectoriel */
	vj[0] = vi[1] * t[2] - t[1] * vi[2];
	vj[1] = vi[2] * t[0] - t[2] * vi[0];
	vj[2] = vi[0] * t[1] - t[0] * vi[1];
	
	local_basis[0] = vi[0];
	local_basis[1] = vi[1];
	local_basis[2] = vi[2];
	local_basis[3] = vj[0];
	local_basis[4] = vj[1];
	local_basis[5] = vj[2];

	local_basis[6] = coo[0];
	local_basis[7] = coo[1];
	local_basis[8] = coo[2];

	local_basis[9] = t[0];
	local_basis[10] = t[1];
	local_basis[11] = t[2];

      } }


  { I iz;
    for (iz=1;iz<Nz;++iz)
      {


	double * __restrict__ local_basis 			= &basis[12*iz];
	double * __restrict__ previous_local_basis 	= &basis[12*(iz-1)];
	
	const R a0 = previous_local_basis[0]*local_basis[0] + previous_local_basis[1]*local_basis[1] + previous_local_basis[2]*local_basis[2];
	const R b0 = previous_local_basis[0]*local_basis[3] + previous_local_basis[1]*local_basis[4] + previous_local_basis[2]*local_basis[5];
	R projection_of_previous_vi[3];
	R projection_of_previous_vj[3];
	projection_of_previous_vi[0]  = a0 * local_basis[0] + b0 * local_basis[3];
	projection_of_previous_vi[1]  = a0 * local_basis[1] + b0 * local_basis[4];
	projection_of_previous_vi[2]  = a0 * local_basis[2] + b0 * local_basis[5];


	Vect_normalized(projection_of_previous_vi,dim);

	R t[3];
	t[0] = local_basis[9];
	t[1] = local_basis[10];
	t[2] = local_basis[11];

	projection_of_previous_vj[0] = projection_of_previous_vi[1] * t[2] - t[1] * projection_of_previous_vi[2];
	projection_of_previous_vj[1] = projection_of_previous_vi[2] * t[0] - t[2] * projection_of_previous_vi[0];
	projection_of_previous_vj[2] = projection_of_previous_vi[0] * t[1] - t[0] * projection_of_previous_vi[1];
	Vect_normalized(projection_of_previous_vj,dim);	

	/* maintenant on change la base */
	
	local_basis[0] = projection_of_previous_vi[0];
	local_basis[1] = projection_of_previous_vi[1];
	local_basis[2] = projection_of_previous_vi[2];
	local_basis[3] = projection_of_previous_vj[0];
	local_basis[4] = projection_of_previous_vj[1];
	local_basis[5] = projection_of_previous_vj[2];

      } }
  


  { I iz;
    for (iz=0;iz<Nz;++iz)
      {
	double * __restrict__ local_basis = &basis[12*iz];
	/* 
	   on traverse tous les noeuds de la frame 
	*/	
	{ R aa[3];
	  R bb[3];
	  aa[0] = local_basis[0];
	  aa[1] = local_basis[1];
	  aa[2] = local_basis[2];

	  bb[0] = local_basis[3];
	  bb[1] = local_basis[4];
	  bb[2] = local_basis[5];

	  const R nx = local_basis[9];
	  const R ny = local_basis[10];
	  const R nz = local_basis[11];	
	  /* 
	     on effectue une rotation autour de la normale d'angle PI/2 de vj
	     si on tombe sur vi on change de signe
	     on veut un repere correcte ou vj est obtenu par rotation de PI	   
	  */
	  const R sx = aa[0] + nx*nx * bb[0] + (nx*ny - nz)*bb[1] + (nx*nz + ny)*bb[2];
	  const R sy = aa[1] + (nx*ny+nz) * bb[0] + (ny*ny)*bb[1] + (ny*nz - nx)*bb[2];
	  const R sz = aa[2] + (nx*nz-ny) * bb[0] + (ny*nz + nx)*bb[1] + (nz*nz)*bb[2];
	  if (sx*sx+sy*sy+sz*sz>1.0)
	    {
	      local_basis[3]*=((R)-1.0);
	      local_basis[4]*=((R)-1.0);
	      local_basis[5]*=((R)-1.0);
	    } }
#if 0
	R rr = ((R)2.0)*rayon[iz];
	rr = ((R)2.0)*minrayon;
#endif
	R rr = ((R)2.0)*rayon[iz]*((R)0.5);
#if 1
	printf("minrayon %e\n",minrayon);
#endif
	{ I iv;
	  for (iv=0;iv<offsetNode_[0];++iv)
	    {
	      I idTopology 	= 0;
	      MeshNode_get_vertex(meshNode_,
				  iz*offsetNode_[0] + iv,
				  &idTopology,
				  coo0);
#if 0
	      coo0[0] = ( (coo0[0] + ((R)1.0))*((R)0.5) - ((R)0.5)) *0.325;
	      if (dim > __wmesh_int_t _2)
		{
		  coo0[1] = ( (coo0[1] + ((R)1.0))*((R)0.5) -((R)0.5)) *0.325;
		}
#endif
	      coo0[0] = ( (coo0[0] + ((R)1.0))*((R)0.5) - ((R)0.5))*rr;
	      if (dim > __wmesh_int_t _2)
		{
		  coo0[1] = ( (coo0[1] + ((R)1.0))*((R)0.5) -((R)0.5))*rr;
		}
	      const R xx = coo0[0];
	      const R yy = coo0[1];			
	      coo0[0]    = local_basis[6] + (xx * local_basis[0] + yy * local_basis[3]);
	      coo0[1]    = local_basis[7] + (xx * local_basis[1] + yy * local_basis[4]);
	      coo0[2]    = local_basis[8] + (xx * local_basis[2] + yy * local_basis[5]);	      
	      MeshNode_set_vertex(meshNode_,
				  iz*offsetNode_[0] + iv,
				  idTopology,
				  coo0);
	    } }

      } }
}






#ifdef __cplusplus
}
#endif
#endif
