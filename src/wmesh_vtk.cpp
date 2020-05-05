#include <array>
#include <string.h>
#include <stdarg.h>
#include <valarray>
#include <iostream>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"
#include "libmeshb7.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std::chrono;


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

extern "C"
{
  
  wmesh_status_t wmesh_write_vtk(const wmesh_t* 		self__,
				 const char * 		filename_,
				 ...)
  {
    wmesh_t*self_ = (wmesh_t*)self__;

#if 0
    {
      double p0[3] {1.0,1.0,0.2};
      double p1[3] {1.0,1.0,1.0};
      double p2[3] {1.0,0.0,1.0};
      double p3[3] {0.7,0.0,-0.3};
      double p4[3] {2.0,0.9,0.7};
      
      for (wmesh_int_t i=0;i<self_->m_num_nodes;++i)
	{
	  double r = self_->m_coo[3*i+0];
	  double s = self_->m_coo[3*i+1];
	  double t = self_->m_coo[3*i+2];
	  
	  double a0 = (t < 1.0) ? (1.0 - t  - r) * (1.0-t-s) / (1.0-t) : 0;
	  double a1 = (t < 1.0) ? (r * (1.0-t-s)) / (1.0-t) : 0;
	  double a2 = (t < 1.0) ? (r * s) / (1.0-t) : 0;
	  double a3 = (t < 1.0) ? ((1.0 - t - r) * s) / (1.0-t) : 0;
	  double a4 = t;
	  
	  double x = a0 * p0[0] + a1 * p1[0] + a2 * p2[0] + a3 * p3[0] + a4 * p4[0];
	  double y = a0 * p0[1] + a1 * p1[1] + a2 * p2[1] + a3 * p3[1] + a4 * p4[1];
	  double z = a0 * p0[2] + a1 * p1[2] + a2 * p2[2] + a3 * p3[2] + a4 * p4[2];
	  
	  self_->m_coo[3*i+0] = x;
	  self_->m_coo[3*i+1] = y;
	  self_->m_coo[3*i+2] = z;
	}

    }
#endif
    char filename[256];    
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
	out << self_->m_coo[3*i+0] << " " << self_->m_coo[3*i+1] << " " << self_->m_coo[3*i+2] << std::endl;
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
