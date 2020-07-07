#include <array>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <valarray>
#include <iostream>
#include "wmesh-types.hpp"
#include "wmesh-status.h"

#include <chrono>
//#include "wfe_element.h"
#include "bms.h"
#include <array>
#include <string.h>
#include <stdarg.h>
#include <valarray>
#include <iostream>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_t.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std::chrono;

using namespace std::chrono;

extern "C"
{

  wmesh_status_t wmesh_write_medit(const wmesh_t* 		self_,
				   bool 			is_binary_,
				   const char * 		filename_,
				   ...)
  {
    WMESH_CHECK_POINTER(self_);
    WMESH_CHECK_POINTER(filename_);

    wmesh_status_t status;    
    wmesh_str_t filename;
    
    { va_list args;
      va_start(args,filename_);
      vsprintf(filename,filename_,args);
      va_end(args); }

    int64_t inm;
    int32_t version;
    if (is_binary_)
      {
#if WMESH_ILP64
	version = 4;
#else
	version = 3;
#endif
      }
    else
      {
	version = 1;
      }
    
    status = bms_write_medit_open(&inm,filename,version,self_->m_coo_m);
    WMESH_STATUS_CHECK(status);
    
    //
    // Write the geometry.
    //
    wmesh_int_t 		coo_m  	= self_->m_coo_m;
    wmesh_int_t 		coo_n  	= self_->m_coo_n;
    wmesh_int_t 		coo_ld 	= self_->m_coo_ld;
    const double*__restrict__ 	coo_v 	= self_->m_coo;
    
    status = bms_write_medit_geometry(inm,
				      coo_m,
				      coo_n,
				      coo_v,
				      coo_ld,
				      self_->m_n_c.v,
				      self_->m_n_c.ld);
    WMESH_STATUS_CHECK(status);

    status = bms_write_medit_topology(inm,

				      self_->m_c2n.m_size,
				      self_->m_c2n.m_ptr,
				      self_->m_c2n.m_m,
				      self_->m_c2n.m_n,
				      self_->m_c2n.m_data,
				      self_->m_c2n.m_ld,
				      
				      self_->m_c_c.m_size,
				      self_->m_c_c.m_ptr,
				      self_->m_c_c.m_m,
				      self_->m_c_c.m_n,
				      self_->m_c_c.m_data,
				      self_->m_c_c.m_ld);    
    WMESH_STATUS_CHECK(status);
    
#if 0    
    if (self_->m_bf2n.m_size > 0)
      {
	status = bms_write_medit_topology(inm,					
					  self_->m_bf2n.m_size,
					  self_->m_bf2n.m_ptr,
					  self_->m_bf2n.m_m,
					  self_->m_bf2n.m_n,
					  self_->m_bf2n.m_data,
					  self_->m_bf2n.m_ld,
					  self_->m_bf_c.m_size,
					  self_->m_bf_c.m_ptr,
					  self_->m_bf_c.m_m,
					  self_->m_bf_c.m_n,
					  self_->m_bf_c.m_data,
					  self_->m_bf_c.m_ld);
	
	WMESH_STATUS_CHECK(status);
      }
#endif
    //
    // Close.
    //
    status = bms_medit_close(inm);
    WMESH_STATUS_CHECK(status);
    
    return WMESH_STATUS_SUCCESS;  
  };

  wmesh_status_t wmesh_read_medit(wmesh_t** 		self__,
				  bool 			is_binary_,
				  const char * 		filename_,
				  ...)
  {

    wmesh_status_t status;
    char filename[256];    
    { va_list args;
      va_start(args,filename_);
      vsprintf(filename,filename_,args);
      va_end(args); }
    
    int32_t dim;
    int32_t version;
    int64_t inm;
    status =  bms_read_medit_open(filename,
				  &inm,
				  &version,
				  &dim);
    WMESH_STATUS_CHECK(status);
    int32_t ref_version = 1;
    if (is_binary_)
      {
#if WMESH_ILP64
	ref_version = 4;
#else
	ref_version = 3;
#endif
      }
    else
      {
	ref_version = 1;
      }

    if (version != ref_version)
      {
	std::cerr << "// ERROR::wmesh_read_medit: the version " << version << " is incompatible with the installation one (=" << ref_version << ")" << std::endl;
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
      }    

    wmesh_int_t
      topology_dimension = 0,
      num_nodes = 0,
      num_edges = 0,
      num_triangles = 0,
      num_quadrilaterals = 0,
      num_tetrahedra = 0,
      num_pyramids = 0,
      num_wedges = 0,
      num_hexahedra = 0;

    status =  bms_read_medit_stat(inm,
				  &num_nodes,
				  &num_edges,
				  &num_triangles,
				  &num_quadrilaterals,
				  &num_tetrahedra,
				  &num_pyramids,
				  &num_wedges,
				  &num_hexahedra);

    
    WMESH_STATUS_CHECK(status);

    //    std::cout << "numnodes " << num_nodes << std::endl;

    wmesh_int_t 	coo_m 		= dim;
    wmesh_int_t 	coo_n 		= num_nodes;
    wmesh_int_t 	coo_ld 		= coo_m;
    double * 		coo_v   	= nullptr;
  
    wmesh_int_t 	nflags_ld 	= 1;  
    wmesh_int_p   	nflags_v;

    wmesh_int_t c2n_size = 0;
    wmesh_int_t c2n_ptr[4+1];
    wmesh_int_t c2n_m[4]{4,5,6,8};
    wmesh_int_t c2n_n[4];
    wmesh_int_t c2n_ld[4];
    wmesh_int_p c2n_v = nullptr;

    wmesh_int_t c_c_size = 0;
    wmesh_int_t c_c_ptr[4+1];
    wmesh_int_t c_c_m[4]{1,1,1,1};
    wmesh_int_t c_c_n[4];
    wmesh_int_t c_c_ld[4];
    wmesh_int_p c_c_v = nullptr;
  
    c2n_n[0] = num_tetrahedra;
    c2n_n[1] = num_pyramids;
    c2n_n[2] = num_wedges;
    c2n_n[3] = num_hexahedra;

    wmesh_int_t bf2n_size = 0;
    wmesh_int_t bf2n_ptr[2+1];
    wmesh_int_t bf2n_m[2]{3,4};
    wmesh_int_t bf2n_n[2];
    wmesh_int_t bf2n_ld[2];
    wmesh_int_p bf2n_v = nullptr;


    wmesh_int_t bf_c_size = 0;
    wmesh_int_t bf_c_ptr[4+1];
    wmesh_int_t bf_c_m[4]{1,1,1,1};
    wmesh_int_t bf_c_n[4];
    wmesh_int_t bf_c_ld[4];
    wmesh_int_p bf_c_v = nullptr;

    wmesh_int_t total_num_bfaces = 0;
    wmesh_int_t total_num_cells = c2n_n[0] + c2n_n[1] + c2n_n[2] + c2n_n[3];
    if (total_num_cells > 0)
      {
	topology_dimension = 3;
	c2n_size = 4;
  
	bf2n_n[0] = num_triangles;
	bf2n_n[1] = num_quadrilaterals;
	total_num_bfaces = bf2n_n[0] + bf2n_n[1];
	if (total_num_bfaces > 0)
	  {
	    bf2n_size = 2;
	  }      
      }
    else
      {
	c2n_m[0] = 3;
	c2n_m[1] = 4;
      
	c2n_n[0] = num_triangles;
	c2n_n[1] = num_quadrilaterals;
	
	total_num_cells = c2n_n[0] + c2n_n[1];
	if (total_num_cells > 0)
	  {
	    topology_dimension = 2;
	    c2n_size = 2;
	    bf2n_n[0] = num_edges;
	    total_num_bfaces = bf2n_n[0];
	    if (total_num_bfaces > 0)
	      {
		bf2n_m[0] = 2;
		bf2n_size = 1;
	      }
	  }
	else
	  {
	    c2n_n[0] = num_edges;
	    total_num_cells = c2n_n[0];
	    if (total_num_cells > 0)
	      {
		topology_dimension = 1;
		c2n_size = 1;	  
	      }
	    else
	      {
		return WMESH_STATUS_INVALID_CONFIG;
	      }
	  }
      }

    c_c_size = c2n_size;
    for (wmesh_int_t i=0;i<c2n_size;++i) c_c_n[i] = c2n_n[i];
    for (wmesh_int_t i=0;i<c2n_size;++i) c_c_ld[i] = c_c_m[i];
    c_c_ptr[0] = 0;
    for (wmesh_int_t i=0;i<c2n_size;++i) c_c_ptr[i+1] = c_c_ptr[i] + c_c_n[i] * c_c_ld[i];

    bf_c_size = bf2n_size;
    for (wmesh_int_t i=0;i<bf2n_size;++i) bf_c_n[i] = bf2n_n[i];
    for (wmesh_int_t i=0;i<bf2n_size;++i) bf_c_ld[i] = bf_c_m[i];
    bf_c_ptr[0] = 0;
    for (wmesh_int_t i=0;i<bf2n_size;++i) bf_c_ptr[i+1] = bf_c_ptr[i] + bf_c_n[i] * bf_c_ld[i];

  
    for (wmesh_int_t i=0;i<c2n_size;++i)
      {
	c2n_ld[i] = c2n_m[i];
      }

    c2n_ptr[0] = 0;
    for (wmesh_int_t i=0;i<c2n_size;++i)
      {
	c2n_ptr[i+1] = c2n_ptr[i] + c2n_ld[i] * c2n_n[i];
      }
  
    for (wmesh_int_t i=0;i<c2n_size;++i)
      {
	bf2n_ld[i] = bf2n_m[i];
      }
  
    bf2n_ptr[0] = 0;
    for (wmesh_int_t i=0;i<bf2n_size;++i)
      {
	bf2n_ptr[i+1] = bf2n_ptr[i] + bf2n_ld[i] * bf2n_n[i];
      }

    if (bf2n_size > 0)
      {
	bf2n_v 	= (wmesh_int_p)malloc(bf2n_ptr[bf2n_size] * sizeof(wmesh_int_t));
	bf_c_v 	= (wmesh_int_p)malloc(bf_c_ptr[bf_c_size] * sizeof(wmesh_int_t));
      }
  
    c2n_v 	= (wmesh_int_p)malloc(c2n_ptr[c2n_size] * sizeof(wmesh_int_t));
    c_c_v 	= (wmesh_int_p)malloc(c_c_ptr[c_c_size]*sizeof(wmesh_int_t));
  
    nflags_v 	= (wmesh_int_p)malloc(num_nodes*sizeof(wmesh_int_t));
    coo_v   	= (double*)malloc((3*num_nodes)*sizeof(double));
    if (!coo_v)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }
    
    status = bms_read_medit_geometry(inm,
				     coo_m,
				     coo_n,
				     coo_v,
				     coo_ld,
				     nflags_v,
				     nflags_ld);
    WMESH_STATUS_CHECK(status);

    status = bms_read_medit_topology(inm,
				   
				     c2n_size,
				     c2n_ptr,
				     c2n_m,
				     c2n_n,
				     c2n_v,
				     c2n_ld,
				   
				     c_c_size,
				     c_c_ptr,
				     c_c_m,
				     c_c_n,
				     c_c_v,
				     c_c_ld);
    WMESH_STATUS_CHECK(status);

    if (bf2n_size > 0)
      {      
	status = bms_read_medit_topology(inm,
				       
					 bf2n_size,
					 bf2n_ptr,
					 bf2n_m,
					 bf2n_n,
					 bf2n_v,
					 bf2n_ld,
					 
					 bf_c_size,
					 bf_c_ptr,
					 bf_c_m,
					 bf_c_n,
					 bf_c_v,
					 bf_c_ld);
	WMESH_STATUS_CHECK(status);      
      }
  
    status = bms_medit_close(inm);
    WMESH_STATUS_CHECK(status);      
    
    status =  wmesh_factory(self__,
			    topology_dimension,

			    c2n_size,
			    c2n_ptr,
			    c2n_m,
			    c2n_n,
			    c2n_v,
			    c2n_ld,
			    
			    c_c_size,
			    c_c_ptr,
			    c_c_m,
			    c_c_n,
			    c_c_v,
			    c_c_ld,
				 
			    bf2n_size,
			    bf2n_ptr,
			    bf2n_m,
			    bf2n_n,
			    bf2n_v,
			    bf2n_ld,

			    bf_c_size,
			    bf_c_ptr,
			    bf_c_m,
			    bf_c_n,
			    bf_c_v,
			    bf_c_ld,				 
			    
			    dim,
			    num_nodes,
			    coo_v,
			    coo_ld,

			    nflags_v,
			    1);
    
    return WMESH_STATUS_SUCCESS;  
  };





  typedef enum __eVtkElement{ ERROR=0,
			      VERTEX=1,
			      POLYVERTEX=2,
			      LINE=3,
			      POLYLINE=4,
			      TRIANGLE=5,
			      TRIANGLE_STRIP=6,
			      POLYGON=7,
			      PIXEL=8,
			      QUAD=9,
			      TETRA=10,
			      VOXEL=11,
			      HEXA=12,
			      WEDGE=13,
			      PYRAMID=14,
			      QUADRATICEDGE=21,
			      QUADRATICTRIANGLE=22,
			      QUADRATICQUAD=23,
			      QUADRATICTETRA=24,
			      QUADRATICHEXA=25,
			      ALL} eVtkElement;

  
  
    wmesh_status_t wmesh_write_vtk(const wmesh_t* 	self__,
				   const char * 		filename_,
				   ...)
    {
      wmesh_t*self_ = (wmesh_t*)self__;

      wmesh_str_t filename;
    
      { va_list args;
	va_start(args,filename_);
	vsprintf(filename,filename_,args);
	va_end(args); }
    
      std::ofstream out(filename);
      out.precision(15);
      out.setf(std::ios::scientific);
    
      out << "# vtk DataFile Version 3.0"
	  << std::endl
	  << "vtk output"
	  << std::endl
	  << "ASCII"
	  << std::endl
	  << "DATASET UNSTRUCTURED_GRID"
	  << std::endl
	  << "POINTS "
	  << self_->m_num_nodes
	  << " double"
	  << std::endl;

      // geometry
      for (wmesh_int_t i=0;i<self_->m_num_nodes;++i)
	{
	  for (wmesh_int_t j=0;j<self_->m_coo_m;++j)
	    {
	      out << " " << self_->m_coo[self_->m_coo_ld*i+j];
	    }
	  out << std::endl;
	}
    
      unsigned long long int total_numCellsData = 0;
      wmesh_int_t num_cells = 0;
      for (wmesh_int_t i=0;i<self_->m_c2n.m_size;++i)
	{
	  num_cells += self_->m_c2n.m_n[i];
	}
      // topology
      for (wmesh_int_t i=0;i<self_->m_c2n.m_size;++i)
	{
	  wmesh_int_t m = self_->m_c2n.m_m[i];
	  wmesh_int_t n = self_->m_c2n.m_n[i];
	  total_numCellsData += n * ( 1 + m );
	}
      
      out << "CELLS "
	  << num_cells
	  << " "
	  << total_numCellsData
	  << std::endl;
    
      for (wmesh_int_t i=0;i<self_->m_c2n.m_size;++i)
	{
	  wmesh_int_t m = self_->m_c2n.m_m[i];
	  wmesh_int_t n = self_->m_c2n.m_n[i];
	  wmesh_int_t ld = self_->m_c2n.m_ld[i];
	  wmesh_int_t ptr = self_->m_c2n.m_ptr[i];
	  const_wmesh_int_p x = self_->m_c2n.m_data + ptr;
	  for (wmesh_int_t idx =0;idx<n;++idx)
	    {
	      out << m;
	      for (wmesh_int_t j=0;j<m;++j)
		{
		  out << " " << x[idx * ld + j]-1;
		}
	      out << std::endl;
	    }	
	}


    
      out << "CELL_TYPES " << num_cells << std::endl;

    
      if (self_->m_topology_dimension == 3)
	{    
	  if (self_->m_c2n.m_n[0]>0)
	    {
	      for (wmesh_int_t idx =0;idx<self_->m_c2n.m_n[0];++idx)
		{
		  out << TETRA << std::endl;
		}
	    }
	  if (self_->m_c2n.m_n[1]>0)
	    {
	      for (wmesh_int_t idx =0;idx<self_->m_c2n.m_n[1];++idx)
		{
		  out << PYRAMID << std::endl;
		}
	    }
	  if (self_->m_c2n.m_n[2]>0)
	    {
	      for (wmesh_int_t idx =0;idx<self_->m_c2n.m_n[2];++idx)
		{
		  out << WEDGE << std::endl;
		}
	    }
	  if (self_->m_c2n.m_n[3]>0)
	    {
	      for (wmesh_int_t idx =0;idx<self_->m_c2n.m_n[3];++idx)
		{
		  out << HEXA << std::endl;
		}
	    }
	}
      else if (self_->m_topology_dimension == 2)
	{
	  if (self_->m_c2n.m_n[0]>0)
	    {
	      for (wmesh_int_t idx =0;idx<self_->m_c2n.m_n[0];++idx)
		{
		  out << TRIANGLE << std::endl;
		}
	    }
	  if (self_->m_c2n.m_n[1]>0)
	    {
	      for (wmesh_int_t idx =0;idx<self_->m_c2n.m_n[1];++idx)
		{
		  out << QUAD << std::endl;
		}
	    }
	}
      else if (self_->m_topology_dimension == 1)
	{
	  if (self_->m_c2n.m_n[0]>0)
	    {
	      for (wmesh_int_t idx =0;idx<self_->m_c2n.m_n[0];++idx)
		{
		  out << LINE << std::endl;
		}
	    }
	}
    
      out << "CELL_DATA " << num_cells << std::endl;
      out << "SCALARS scalars int 1"  << std::endl;
      out << "LOOKUP_TABLE default"  << std::endl;

      for (wmesh_int_t s = 0;s<self_->m_c_c.m_size ;++s)
	{
	  for (wmesh_int_t j = 0;j < self_->m_c_c.m_n[ s ];++j)
	    {
	      for (wmesh_int_t i = 0;i < self_->m_c_c.m_m[ s ];++i)
		{
		  out << self_->m_c_c.m_data[self_->m_c_c.m_ptr[s] + self_->m_c_c.m_ld[s] * j + i ] << std::endl;
		}
	    }
	}

      out.close();
      return WMESH_STATUS_SUCCESS;  
    }

  

}
