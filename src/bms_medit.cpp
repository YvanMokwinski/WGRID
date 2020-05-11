#include <array>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <valarray>
#include <iostream>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"
#include "libmeshb7.h"
#include <chrono>


namespace Medit
{
  template<typename _int_t> static int GmfType();
  template<typename _int_t> static int GmfVecType();
  template<> int GmfType<int>()		{return GmfInt;};
  template<> int GmfType<long>()	{return GmfLong;};
  template<> int GmfType<long long>()	{return GmfLong;};
  template<> int GmfType<double>()	{return GmfDouble;};
  template<> int GmfType<float>()	{return GmfFloat;};

  template<> int GmfVecType<int>()		{return GmfIntVec;};
  template<> int GmfVecType<long>()		{return GmfLongVec;};
  template<> int GmfVecType<long long>()	{return GmfLongVec;};
  template<> int GmfVecType<double>()		{return GmfDoubleVec;};
  template<> int GmfVecType<float>()		{return GmfFloatVec;};
};

#ifdef GMF_CHECK_STATUS
#error GMF_CHECK_STATUS already defined.
#endif
#define GMF_CHECK_STATUS(_status) if (_status != 1) return 1



#if 0
template<typename T>
static wmesh_status_t bms_read_medit_geometry(int64_t 				inm_,
					      wmesh_int_t			coo_m_,
					      wmesh_int_t			coo_n_,
					      T *__restrict__ 			coo_,
					      wmesh_int_t 			coo_ld_,
					      wmesh_int_p  			nflags_,
					      wmesh_int_t 			nflags_ld_)
{
  WMESH_CHECK_POINTER(coo_);
  WMESH_CHECK_POINTER(nflags_);

  int gmf_status;

  T*__restrict__ coo_first	= coo_;
  T*__restrict__ coo_last	= coo_ + (coo_n_-1) * coo_ld_;    
  const_wmesh_int_p nflags_first 	= nflags_;
  const_wmesh_int_p nflags_last 	= &nflags_[ (coo_n_-1) * nflags_ld_ ];
  
  gmf_status = GmfSetKwd(inm_,
			 GmfVertices,
			 coo_n_);
  GMF_CHECK_STATUS(gmf_status);
  switch (coo_m_)
    {
    case 3:
      {
	gmf_status = GmfSetBlock(inm_,
				 GmfVertices,
				 1,
				 coo_n_,
				 0,
				 NULL,
				 NULL,
				 Medit::GmfType<T>(), coo_first + 0, coo_last + 0,
				 Medit::GmfType<T>(), coo_first + 1, coo_last + 1,
				 Medit::GmfType<T>(), coo_first + 2, coo_last + 2,
				 Medit::GmfType<wmesh_int_t>(), nflags_first, nflags_last);
	GMF_CHECK_STATUS(gmf_status);
	break;
      }
    case 2:
      {
	gmf_status = GmfSetBlock(inm_,
				 GmfVertices,
				 1,
				 coo_n_,
				 0,
				 NULL,
				 NULL,
				 Medit::GmfType<T>(), coo_first + 0, coo_last + 0,
				 Medit::GmfType<T>(), coo_first + 1, coo_last + 1,
				 Medit::GmfType<wmesh_int_t>(), nflags_first, nflags_last);
	GMF_CHECK_STATUS(gmf_status);
      }
    default:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
	break;
      }
    }
  
  return WMESH_STATUS_SUCCESS;  
}
#endif

extern "C"
{
  static int s_gmf_types[5][4] = { {0,0,0,0},
				 {GmfEdges,0,0,0},
				 {GmfTriangles,GmfQuadrilaterals,0,0},
				 {0,0,0,0},
				 {GmfTetrahedra,GmfPyramids,GmfPrisms,GmfHexahedra}};
  wmesh_status_t
  bms_write_medit_open
  (int64_t*		inm_,
   wmesh_str_t 		filename_,
   wmesh_int_t          precision_,
   wmesh_int_t          dimension_)
  {
    inm_[0] = 0;
    int64_t inm;      
    inm = GmfOpenMesh(filename_,
		      GmfWrite,
		      precision_,
		      dimension_);    
    if (!inm)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
      }
    inm_[0] = inm;
    return WMESH_STATUS_SUCCESS;  
  }

  wmesh_status_t bms_read_medit_stat(int64_t inm_,
				     wmesh_int_p num_nodes_,
				     wmesh_int_p num_edges_,
				     wmesh_int_p num_triangles_,
				     wmesh_int_p num_quadrilaterals_,
				     wmesh_int_p num_tetrahedra_,
				     wmesh_int_p num_pyramids_,
				     wmesh_int_p num_wedges_,
				     wmesh_int_p num_hexahedra_)
  {
    WMESH_CHECK((inm_!=0));
    WMESH_CHECK_POINTER( num_nodes_);
    WMESH_CHECK_POINTER( num_edges_);
    WMESH_CHECK_POINTER( num_triangles_);
    WMESH_CHECK_POINTER( num_quadrilaterals_);
    WMESH_CHECK_POINTER( num_tetrahedra_);
    WMESH_CHECK_POINTER( num_pyramids_);
    WMESH_CHECK_POINTER( num_wedges_);
    WMESH_CHECK_POINTER( num_hexahedra_);

    num_nodes_[0] = GmfStatKwd(inm_,GmfVertices);
    num_edges_[0] = GmfStatKwd(inm_,GmfEdges);
    num_triangles_[0] = GmfStatKwd(inm_,GmfTriangles);
    num_quadrilaterals_[0] = GmfStatKwd(inm_,GmfQuadrilaterals);
    num_tetrahedra_[0] = GmfStatKwd(inm_,GmfTetrahedra);
    num_pyramids_[0] = GmfStatKwd(inm_,GmfPyramids);
    num_wedges_[0] = GmfStatKwd(inm_,GmfPrisms);
    num_hexahedra_[0] = GmfStatKwd(inm_,GmfHexahedra);
    return WMESH_STATUS_SUCCESS;  

  }
				     
  wmesh_status_t bms_read_medit_open(wmesh_str_t 	filename_,
				     int64_t*		inm_,
				     int32_t * 		version_,
				     int32_t * 		dim_)
  {
    inm_[0] = 0;
    int64_t inm;
    inm = GmfOpenMesh(filename_,
		      GmfRead,
		      version_,
		      dim_);    
    if (!inm)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
      }
    inm_[0] = inm;
    return WMESH_STATUS_SUCCESS;  
  }

  wmesh_status_t bms_medit_close(int64_t		inm_)
  {
    int gmf_status = GmfCloseMesh(inm_);
    GMF_CHECK_STATUS(gmf_status);
    return WMESH_STATUS_SUCCESS;  
  }
  
  wmesh_status_t bms_write_medit_topology(int64_t				inm_,
					  
					  const wmesh_int_t 			c2n_size,
					  const_wmesh_int_p 			c2n_ptr,
					  const_wmesh_int_p 			c2n_m,
					  const_wmesh_int_p 			c2n_n,
					  const_wmesh_int_p 			c2n_v,
					  const_wmesh_int_p 			c2n_ld,
					  
					  const wmesh_int_t 			c_c_size,
					  const_wmesh_int_p 			c_c_ptr,
					  const_wmesh_int_p 			c_c_m,
					  const_wmesh_int_p 			c_c_n,
					  const_wmesh_int_p 			c_c_v,
					  const_wmesh_int_p 			c_c_ld)
  {
    WMESH_CHECK_POINTER(c2n_v);
    WMESH_CHECK_POINTER(c_c_v);
    
    for (wmesh_int_t itype=0;itype<c2n_size;++itype)
      {
	const wmesh_int_t num_cells = c2n_n[itype];

	if (num_cells > 0)
	  {
	
	    const_wmesh_int_p c2n_first 	= &c2n_v[ c2n_ptr[itype] ];
	    const_wmesh_int_p c2n_last	= &c2n_v[ c2n_ptr[itype] + (num_cells-1) * c2n_ld[itype] ];
	    
	    const_wmesh_int_p c_c_first 	= &c_c_v[ c_c_ptr[itype] ];
	    const_wmesh_int_p c_c_last 	= &c_c_v[ c_c_ptr[itype] + (num_cells-1) * c_c_ld[itype] ];
	    
	    const int gmf_type = s_gmf_types[c2n_size][itype];	
	    
	    GmfSetKwd(inm_,
		      gmf_type,
		      num_cells);

	    GmfSetBlock(inm_,
			gmf_type,
			1,
			num_cells,
			  0,
			NULL,
			NULL,
			
			Medit::GmfVecType<wmesh_int_t>(), c2n_m[itype],c2n_first + 0, c2n_last + 0,
			Medit::GmfType<wmesh_int_t>(), c_c_first, c_c_last);
	  }
      }
    return WMESH_STATUS_SUCCESS;  
  };

  wmesh_status_t bms_write_medit_geometry(int64_t 			inm_,
						 wmesh_int_t			coo_m_,
						 wmesh_int_t			coo_n_,
						 const double *__restrict__ 	coo_,
						 wmesh_int_t 			coo_ld_,
						 const_wmesh_int_p  		nflags_,
						 wmesh_int_t 			nflags_ld_)
  {
    int gmf_status;
    WMESH_CHECK( (inm_!=0) );
    WMESH_CHECK( (coo_m_ >= 2) );
    WMESH_CHECK( (coo_n_ >= 2) );
    WMESH_CHECK( (coo_ld_ >= coo_m_) );
    WMESH_CHECK_POINTER(coo_);
    WMESH_CHECK_POINTER(nflags_);
    WMESH_CHECK( (nflags_ld_ >= 1) );
    
    const double*__restrict__ coo_first	= coo_;
    const double*__restrict__ coo_last	= coo_ + (coo_n_-1) * coo_ld_;    
    const_wmesh_int_p nflags_first 	= nflags_;
    const_wmesh_int_p nflags_last 	= &nflags_[ (coo_n_-1) * nflags_ld_ ];

    gmf_status = GmfSetKwd(inm_,GmfVertices,coo_n_);
    switch (coo_m_)
      {
      case 3:
	{
	  gmf_status = GmfSetBlock(inm_,
				   GmfVertices,
				   1,
				   coo_n_,
				   0,
				   NULL,
				   NULL,
				   Medit::GmfType<double>(), coo_first + 0, coo_last + 0,
				   Medit::GmfType<double>(), coo_first + 1, coo_last + 1,
				   Medit::GmfType<double>(), coo_first + 2, coo_last + 2,
				   Medit::GmfType<wmesh_int_t>(), nflags_first, nflags_last);
	  break;
	}
      case 2:
	{
	  gmf_status = GmfSetBlock(inm_,
				   GmfVertices,
				   1,
				   coo_n_,
				   0,
				   NULL,
				   NULL,
				   Medit::GmfType<double>(), coo_first + 0, coo_last + 0,
				   Medit::GmfType<double>(), coo_first + 1, coo_last + 1,
				   Medit::GmfType<wmesh_int_t>(), nflags_first, nflags_last);
	  break;
	}
      default:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
	  break;
	}
      }
    GMF_CHECK_STATUS(gmf_status);    
    return WMESH_STATUS_SUCCESS;  
  }
  

   wmesh_status_t bms_read_medit_topology(int64_t			inm_,
						
					  wmesh_int_t 			c2n_size,
					  const_wmesh_int_p 		c2n_ptr,
					  const_wmesh_int_p 		c2n_m,
					  const_wmesh_int_p 		c2n_n,
					  wmesh_int_p 			c2n_v,
					  const_wmesh_int_p 		c2n_ld,
						
					  wmesh_int_t 			c_c_size,
					  const_wmesh_int_p 		c_c_ptr,
					  const_wmesh_int_p 		c_c_m,
					  const_wmesh_int_p 		c_c_n,
					  wmesh_int_p 			c_c_v,
					  const_wmesh_int_p 		c_c_ld)
  
  {
    WMESH_CHECK( (inm_!=0) );

    WMESH_CHECK( (c2n_size > 0) );
    WMESH_CHECK_POINTER( c2n_ptr );
    WMESH_CHECK_POINTER( c2n_m );
    WMESH_CHECK_POINTER( c2n_n );
    WMESH_CHECK_POINTER( c2n_v );
    WMESH_CHECK_POINTER( c2n_ld );

    WMESH_CHECK( (c2n_size == c_c_size) );
    WMESH_CHECK_POINTER( c_c_ptr );
    WMESH_CHECK_POINTER( c_c_m );
    WMESH_CHECK_POINTER( c_c_n );
    WMESH_CHECK_POINTER( c_c_v );
    WMESH_CHECK_POINTER( c_c_ld );
    
    int gmf_status;    
    for (wmesh_int_t itype=0;itype<c2n_size;++itype)
      {
	const wmesh_int_t num_cells 	= c2n_n[itype];
	if (num_cells > 0)
	  {
	    const int gmf_type 		= s_gmf_types[c2n_size][itype];	
	    
	    wmesh_int_p c2n_v_first 	= c2n_v + c2n_ptr[itype];
	    wmesh_int_p c_c_v_first 	= c_c_v + c_c_ptr[itype];
	    wmesh_int_p c2n_v_last 	= c2n_v_first + (num_cells-1) * c2n_ld[itype];
	    wmesh_int_p c_c_v_last 	= c_c_v_first + (num_cells-1) * c_c_ld[itype];

	    GmfGotoKwd(inm_,gmf_type);
	    switch(c2n_m[itype])
	      {
	      case 2:
		{
		  gmf_status = GmfGetBlock(inm_,
					   gmf_type,
					   1,
					   num_cells,
					   0,
					   NULL,
					   NULL,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 0, c2n_v_last + 0,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 1, c2n_v_last + 1,
					   Medit::GmfType<wmesh_int_t>(), c_c_v_first + 0, c_c_v_last + 0);
		  break;
		}
	      case 3:
		{
		  gmf_status = GmfGetBlock(inm_,
					   gmf_type,
					   1,
					   num_cells,
					   0,
					   NULL,
					   NULL,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 0, c2n_v_last + 0,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 1, c2n_v_last + 1,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 2, c2n_v_last + 2,
					   Medit::GmfType<wmesh_int_t>(), c_c_v_first + 0, c_c_v_last + 0);
		  break;
		}
	      case 4:
		{
		  gmf_status = GmfGetBlock(inm_,
					   gmf_type,
					   1,
					   num_cells,
					   0,
					   NULL,
					   NULL,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 0, c2n_v_last + 0,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 1, c2n_v_last + 1,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 2, c2n_v_last + 2,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 3, c2n_v_last + 3,
					   Medit::GmfType<wmesh_int_t>(), c_c_v_first + 0, c_c_v_last + 0);
		  break;
		}
	  
	      case 5:
		{
		  gmf_status = GmfGetBlock(inm_,
					   gmf_type,
					   1,
					   num_cells,
					   0,
					   NULL,
					   NULL,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 0, c2n_v_last + 0,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 1, c2n_v_last + 1,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 2, c2n_v_last + 2,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 3, c2n_v_last + 3,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 4, c2n_v_last + 4,
					   Medit::GmfType<wmesh_int_t>(), c_c_v_first + 0, c_c_v_last + 0);
		  break;
		}

	      case 6:
		{
		  gmf_status = GmfGetBlock(inm_,
					   gmf_type,
					   1,
					   num_cells,
					   0,
					   NULL,
					   NULL,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 0, c2n_v_last + 0,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 1, c2n_v_last + 1,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 2, c2n_v_last + 2,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 3, c2n_v_last + 3,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 4, c2n_v_last + 4,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 5, c2n_v_last + 5,
					   Medit::GmfType<wmesh_int_t>(), c_c_v_first + 0, c_c_v_last + 0);
		  break;
		}
	  
	      case 8:
		{
		  gmf_status = GmfGetBlock(inm_,
					   gmf_type,
					   1,
					   num_cells,
					   0,
					   NULL,
					   NULL,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 0, c2n_v_last + 0,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 1, c2n_v_last + 1,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 2, c2n_v_last + 2,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 3, c2n_v_last + 3,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 4, c2n_v_last + 4,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 5, c2n_v_last + 5,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 6, c2n_v_last + 6,
					   Medit::GmfType<wmesh_int_t>(), c2n_v_first + 7, c2n_v_last + 7,
					   Medit::GmfType<wmesh_int_t>(), c_c_v_first + 0, c_c_v_last + 0);
		  break;
		}
		
	      }
	    GMF_CHECK_STATUS(gmf_status);
	  }
      }
    
    return WMESH_STATUS_SUCCESS;
  }
  
  wmesh_status_t bms_read_medit_geometry(int64_t 			inm_,
					 wmesh_int_t			coo_m_,
					 wmesh_int_t			coo_n_,
					 double *__restrict__ 		coo_,
					 wmesh_int_t 			coo_ld_,
					 wmesh_int_p  			nflags_,
					 wmesh_int_t 			nflags_ld_)
  {
    WMESH_CHECK( (inm_!=0) );
    WMESH_CHECK( (coo_m_ >= 2) );
    WMESH_CHECK( (coo_n_ >= 2) );
    WMESH_CHECK( (coo_ld_ >= coo_m_) );
    WMESH_CHECK_POINTER(coo_);
    WMESH_CHECK_POINTER(nflags_);
    WMESH_CHECK( (nflags_ld_ >= 1) );
    
    int gmf_status;
    
    double*__restrict__ coo_first	= coo_;
    double*__restrict__ coo_last	= coo_ + (coo_n_-1) * coo_ld_;    
    const_wmesh_int_p nflags_first 	= nflags_;
    const_wmesh_int_p nflags_last 	= &nflags_[ (coo_n_-1) * nflags_ld_ ];
    
    GmfGotoKwd(inm_,GmfVertices);
    switch (coo_m_)
      {
      case 2:
	{
	  gmf_status =  GmfGetBlock(inm_,
				    GmfVertices,
				    1,
				    coo_n_,
				    0,
				    NULL,
				    NULL,
				    Medit::GmfType<double>(), coo_first + 0, coo_last + 0,
				    Medit::GmfType<double>(), coo_first + 1, coo_last + 1,
				    Medit::GmfType<wmesh_int_t>(), nflags_first, nflags_last);
	  break;
	}
	
      case 3:
	{
	  gmf_status = GmfGetBlock(inm_,
				   GmfVertices,
				   1,
				   coo_n_,
				   0,
				   NULL,
				   NULL,
				   Medit::GmfType<double>(), coo_first + 0, coo_last + 0,
				   Medit::GmfType<double>(), coo_first + 1, coo_last + 1,
				   Medit::GmfType<double>(), coo_first + 2, coo_last + 2,
				   Medit::GmfType<wmesh_int_t>(), nflags_first, nflags_last);
	  break;
	}
      default:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
	  break;
	}
      }
    GMF_CHECK_STATUS(gmf_status);    
    return WMESH_STATUS_SUCCESS;  
  }

};
#undef GMF_CHECK_STATUS
