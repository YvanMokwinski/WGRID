#include <iostream>

#include "wmesh-types.hpp"
#include "wmesh.hpp"
#include "wmesh_medit.hpp"

extern "C"
{

};


#if 0


#include <valarray>

#include <iostream>
#include <ostream>
#include "Program.hpp"
#include <limits>
#include <array>

#include "DimensionType.hpp"
#include "ReferenceShape.hpp"

#include "wmesh-types.hpp"
#include "Point.hpp"
#include "CellToNodes.hpp"
#include "Output/Medit.hpp"

#include "Input/Medit.hpp"
#include "WCOMMON/MeasureSteady.hpp"
//#include "Mesh.hpp"
#include "Output/Medit.hpp"
#include "Output/Vtk.hpp"
#include "wmesh_status.h"
#include "AHF/Mesh.hpp"

#include "EdgeCompare.hpp"
#include "FaceCompare.hpp"
#include "generate_finite_element_space.hpp"
#include "CellsToCellsCalculator.hpp"
using namespace std::chrono;

extern "C"
{
  void wint_sparsemat_def1(wint_sparsemat_t*	self_,
			   int_t 	        size_,
			   const int_t 	*	m_,
			   int_t             	ncells_)
  {
    self_->m_size = size_;
    self_->m_m = (int_t*)malloc(sizeof(int_t)*self_->m_size);
    self_->m_n = (int_t*)malloc(sizeof(int_t)*self_->m_size);
    self_->m_ld = (int_t*)malloc(sizeof(int_t)*self_->m_size);
    self_->m_ptr = (int_t*)malloc(sizeof(int_t)*(self_->m_size+1));
    self_->m_data = NULL;
    for (int i=0;i<self_->m_size;++i)
      {
	self_->m_m[i] = m_[i];
      }
    for (int i=0;i<self_->m_size;++i)
      {
	self_->m_n[i] = ncells_;
      }
    for (int i=0;i<self_->m_size;++i)
      {
	self_->m_ld[i] = self_->m_m[i];
      }

    self_->m_ptr[0] = 0;
    for (int i=1;i<=self_->m_size;++i)
      {
	self_->m_ptr[i] = self_->m_ptr[i-1] + self_->m_n[i-1] * self_->m_ld[i-1];
      }      
  }
  
  void wint_sparsemat_info(const wint_sparsemat_t* self_,FILE * out_)
  {

    fprintf(out_,"--- wint_sparsemat_info ---\n");
    fprintf(out_,"--- num segments %Ld\n",self_->m_size);
    for (int i=0;i<self_->m_size;++i)
      {
	fprintf(out_,"--- m %Ld\n",self_->m_m[i]);
      }
    for (int i=0;i<self_->m_size;++i)
      {
	fprintf(out_,"--- n %Ld\n",self_->m_n[i]);
      }
    for (int i=0;i<self_->m_size;++i)
      {
	fprintf(out_,"--- (m,n,ld,offset)  (%Ld,%Ld,%Ld,%Ld)\n",
		self_->m_m[i],
		self_->m_n[i],
		self_->m_ld[i],
		self_->m_ptr[i]);
      }
  }


  void wspace_c2n_info(wspace_c2n_t*	self_,
		      FILE*		out_)
  {
    fprintf(out_,"---- zone info ----\n");
    fprintf(out_,"M  = %Ld\n", self_->m_c2n.m);
    fprintf(out_,"N  = %Ld\n", self_->m_c2n.n);
    fprintf(out_,"LD = %Ld\n", self_->m_c2n.ld);
    fprintf(out_,"-------------------\n");
  }
  
  void wspace_c2n_def(wspace_c2n_t * self_,int_t celltype,int_t c2n_m_,int_t c2n_n_,int_t*c2n_v_,int_t c2n_ld_)
  {
    self_->m_cell_type = celltype;
    wint_mat_def(&self_->m_c2n,
		 c2n_m_,
		 c2n_n_,
		 c2n_v_,
		 c2n_ld_);
  };

  
  void wspace_sparse_c2n_info(wspace_sparse_c2n_t*self_,FILE*out_)
  {
    fprintf(stdout,"---- mesh info ----\n");
    wint_sparsemat_info(&self_->m_data,out_);
    for (int i=0;i<self_->m_nzones;++i)
      {
	wspace_c2n_info(&self_->m_zones[i],out_);
      }
    fprintf(stdout,"-------------------\n");
  }

  void wspace_sparse_c2n_print(wspace_sparse_c2n_t*self_,FILE*out_)
  {    
    fprintf(stdout,"\n");    
  }

  void wspace_sparse_c2n_def(wspace_sparse_c2n_t*	self_,
			     int_t 			nzones_,
			     const int_t*		data_ptr_,
			     const int_t*		data_m_,
			     const int_t*		data_n_,
			     int_t*			data_v_,
			     const int_t*		data_ld_)
    
  {
    wint_sparsemat_def(&self_->m_data,
		       nzones_);
    
    wint_sparsemat_set(&self_->m_data,
		    data_m_,
		    data_n_,
		    data_ld_,
		    data_ptr_,
		    data_v_);

    self_->m_nzones = nzones_;
    self_->m_zones = (wspace_c2n_t*)calloc(nzones_,sizeof(wspace_c2n_t));
    for (int i=0;i<nzones_;++i)
      {	
	wint_sparsemat_get(&self_->m_data,
			   i,
			   &self_->m_zones[i].m_c2n);
	self_->m_zones[i].m_cell_type = i;
      }
    
#if 0

	if (self_->ncells[i]>0)
	  {
	    wspace_c2n_def(&self_->zones[i],i,data_m[i],data_n[i],&data_v[data_ptr[i]],data_ld[i]);
	  }
    self_->ncells = (int_t*)calloc( nzones_,sizeof(int_t) );
    fprintf(stdout,"wspace_sparse_c2n_def\n");
    self_->data.m = 1;
    self_->data.n = 0;
    for (int i=0;i<nzones;++i)
      {
	self_->ncells[i] = data_n[i];
      }

    
    for (int i=0;i<nzones;++i)
      {	
	self_->data.n += data_ld[i] * data_n[i];
      }
    self_->data.v = data_v;
#endif    

  };
  
  void wspace_sparse_c2n_calculate(const wspace_sparse_c2n_t*self_)
  {
    static constexpr int_t degree = 3;
    const int_t nzones = self_->m_nzones;
    int_t size = 0;
    int_t space_c2d_m[4];
    int_t space_c2d_n[4];
    int_t space_c2d_ld[4];
    int_t space_c2d_ptr[5];
    
    {
      space_c2d_ptr[0] = 0;
      for (int_t izone=0;izone<self_->m_nzones;++izone)
	{
	  space_c2d_m[izone] = shape_get_num_dofs<degree>(izone);
	  space_c2d_n[izone] = self_->m_zones[izone].m_c2n.n;
	  space_c2d_ld[izone] = space_c2d_m[izone];
	  size += shape_get_num_dofs<degree>(izone) * self_->m_zones[izone].m_c2n.n;	
	  space_c2d_ptr[izone+1] = size;
	}
      std::cout << "total_numDofs " << size << std::endl;
    }

    printf("--------------------- ALLOC SPACE ------------------------\n");
    int_t * c2d = (int_t*)malloc(sizeof(int_t) * size);    
    wint_sparsemat_t space_c2d;
    wint_sparsemat_def(&space_c2d, nzones);
    wint_sparsemat_set(&space_c2d,
		    space_c2d_m,
		    space_c2d_n,
		    space_c2d_ld,
		    space_c2d_ptr,
		    c2d);
    printf("--------------------- ALLOC SPACE DONE ------------------------\n");


    printf("--------------------- CALCULATE ------------------------\n");

      for (int_t izone=0;izone<self_->m_nzones;++izone)
	{
	  if (izone == 3)
	    {
	      int_t numDofs=0;
	      if (self_->m_zones[3].m_c2n.n>0)
		{
		  wspace_c2d_calculate<degree,VolumeType::Hexahedron,int_t>(self_->m_zones[3].m_c2n.n,
									    self_->m_zones[3].m_c2n.v,
									    self_->m_zones[3].m_c2n.ld,
									    &numDofs,
									    c2d + space_c2d_ptr[3],
									    space_c2d_ld[3]);
		  printf("---------- numdofs hex %Ld\n",numDofs);
		}
	    }
	  
	  if (izone == 0)
	    {
	      int_t numDofs=0;
	      if (self_->m_zones[0].m_c2n.n>0)
		{
		  wspace_c2d_calculate<degree,VolumeType::Tetrahedron,int_t>(self_->m_zones[0].m_c2n.n,
									     self_->m_zones[0].m_c2n.v,
									     self_->m_zones[0].m_c2n.ld,
									     &numDofs,
									     c2d + space_c2d_ptr[0],
									     space_c2d_ld[0]);
		  printf("---------- numdofs tet %Ld\n",numDofs);
		}
	    }
	  
	  if (izone == 2)
	    {
	      int_t numDofs=0;
	      if (self_->m_zones[2].m_c2n.n>0)
		{
		  wspace_c2d_calculate<degree,VolumeType::Wedge,int_t>(self_->m_zones[2].m_c2n.n,
									     self_->m_zones[2].m_c2n.v,
									     self_->m_zones[2].m_c2n.ld,
									     &numDofs,
									     c2d + space_c2d_ptr[2],
									     space_c2d_ld[2]);
		  printf("---------- numdofs wedge %Ld\n",numDofs);
		}
	    }
	}
    printf("---------------------------------------------\n");
  }
  
  void wspace_sparse_c2n_calculate_c2c(const wspace_sparse_c2n_t*self_)
  {

    unsigned long long int 	numFaces = 0;
    unsigned long long int 	numFacesInterior = 0;
    unsigned long long int 	numFacesBoundary = 0;
    int_t size = 0;

    int_t nfaces[4] = {4,5,5,6};
    int_t * c2cz[4];
   for (int i=0;i<self_->m_nzones;++i)
      {
	if (self_->m_zones[i].m_c2n.n > 0)
	  {
	    size += self_->m_zones[i].m_c2n.n * nfaces[i];
	  }
      }  
    int_t*c2c = (int_t*)malloc(sizeof(int_t)*size);
    size = 0;
    for (int i=0;i<self_->m_nzones;++i)
      {
	if (self_->m_zones[i].m_c2n.n > 0)
	  {
	    c2cz[i] = c2c + size;
	    size += self_->m_zones[i].m_c2n.n * nfaces[i];
	  }
      }

    if (c2cz[0]){}
    printf("%Ld %Ld %p nface %Ld\n",self_->m_zones[3].m_c2n.m,self_->m_zones[3].m_c2n.n,c2cz[3],nfaces[3]);
    for (int i=0;i<4;++i)
      {
	for (int j=0;j<8;++j)
	  std::cout << " " << self_->m_zones[3].m_c2n.v[8*i+j];
	std::cout << std::endl;
      }

    

    std::cout << "Calculate c2c hexa" << std::endl;
    if (self_->m_zones[3].m_c2n.n>0)
      {
	calculator_cells_to_cells_t<VolumeType::Hexahedron,FaceType::Quadrilateral,int_t> calc;
	calc.calculate(self_->m_zones[3].m_c2n.n,
		       self_->m_zones[3].m_c2n.v,
		       self_->m_zones[3].m_c2n.ld,
		       c2cz[3],
		       nfaces[3],
		       &numFaces,
		       &numFacesInterior,
		       &numFacesBoundary);		
	std::cout << "Calculate c2c hexa numFaces         " << numFaces << std::endl;
	std::cout << "Calculate c2c hexa numFacesInterior " << numFacesInterior << std::endl;
	std::cout << "Calculate c2c hexa numFacesBoundary " << numFacesBoundary << std::endl;
      }
#if 0
    if (self_->m_zones[0].m_c2n.n>0)
      {
	calculator_cells_to_cells_t<VolumeType::Tetrahedron,FaceType::Triangle,int_t> calc;
	calc.calculate(self_->m_zones[0].m_c2n.n,
		       self_->m_zones[0].m_c2n.v,
		       self_->m_zones[0].m_c2n.ld,
		       c2cz[0],
		       nfaces[0],
		       &numFaces,
		       &numFacesInterior,
		       &numFacesBoundary);		
	std::cout << "Calculate c2c tetra numFaces         " << numFaces << std::endl;
	std::cout << "Calculate c2c tetra numFacesInterior " << numFacesInterior << std::endl;
	std::cout << "Calculate c2c tetra numFacesBoundary " << numFacesBoundary << std::endl;
      }
    
    if (self_->m_zones[2].m_c2n.n>0)
      {
	{
	  calculator_cells_to_cells_t<VolumeType::Wedge,FaceType::Quadrilateral,int_t> calc;
	  calc.calculate(self_->m_zones[2].m_c2n.n,
			 self_->m_zones[2].m_c2n.v,
			 self_->m_zones[2].m_c2n.ld,
			 c2cz[2],
			 nfaces[2],
			 &numFaces,
			 &numFacesInterior,
			 &numFacesBoundary);		
	  std::cout << "Calculate c2c numFaces         " << numFaces << std::endl;
	  std::cout << "Calculate c2c numFacesInterior " << numFacesInterior << std::endl;
	  std::cout << "Calculate c2c numFacesBoundary " << numFacesBoundary << std::endl;
	}
	{
	  calculator_cells_to_cells_t<VolumeType::Wedge,FaceType::Triangle,int_t> calc;
	  calc.calculate(self_->m_zones[2].m_c2n.n,
			 self_->m_zones[2].m_c2n.v,
			 self_->m_zones[2].m_c2n.ld,
			 c2cz[2],
			 nfaces[2],
			 &numFaces,
			 &numFacesInterior,
			 &numFacesBoundary);		
	  std::cout << "Calculate c2c numFaces         " << numFaces << std::endl;
	  std::cout << "Calculate c2c numFacesInterior " << numFacesInterior << std::endl;
	  std::cout << "Calculate c2c numFacesBoundary " << numFacesBoundary << std::endl;
	}
      }
#endif
    int_t N=0;
    for (int_t i=0;i<size;++i)
      {
	if (c2c[i]==0)
	  {	    
	    ++N;
	  }
      }
    std::cout << " xxxxxxxxxxxxxxxxx " << N << std::endl;
    for (int i=0;i<20;++i)
      {
	for (int j=0;j<6;++j)
	  std::cout << " " << c2cz[3][6*i+j];
	std::cout << std::endl;
      }
    
#if 0

    for (int i=0;i<self_->m_nzones;++i)
      {

	if (self_->m_zones[i].m_c2n.n > 0)
	  {
	    
	    if (self_->m_zones[2].m_c2n.n>0)
	      {
		calculator_cells_to_cells_t<VolumeType::Wedge,FaceType::Quadrilateral,int_t>::calculate(self_->m_zones[i].m_c2n.n,
													self_->m_zones[i].m_c2n.v,
													self_->m_zones[i].m_c2n.ld,
													c2c,
													6,
													&numFaces,
													&numFacesInterior,
													&numFacesBoundary);		
	      }
	  }
      }
#endif
#if 0
    wint_sparsemat_t  c2c_data;
    const wspace_c2n_t*zone = &self_->m_zones[3];

    wint_sparsemat_def1(&c2c_data,
		     4,
		     m_,
		     self_->data);
    int_t size = 0;
    if (self_->ncells[3]>0)
      {
	size+=self_->ncells[3]*6;
      }
    if (self_->ncells[2]>0)
      {
	size+=self_->ncells[2]*5;
      }
    if (self_->ncells[1]>0)
      {
	size+=self_->ncells[1]*5;
      }
    if (self_->ncells[0]>0)
      {
	size+=self_->ncells[0]*4;
      }
    
    int_t * c2c_global = (int_t*)malloc(sizeof(int_t)*size);
    int_t * c2c = c2c_global;
    if (self_->ncells[3]>0)
      {
	std::cout << "Calculate c2c" << std::endl;
	calculator_cells_to_cells_t<VolumeType::Hexahedron,FaceType::Quadrilateral,int_t>::calculate(zone->data.n,
												     zone->data.v,
												     zone->data.ld,
												     c2c,
												     6,
												     &numFaces,
												     &numFacesInterior,
												     &numFacesBoundary);
	c2c += 
      }
    if (self_->ncells[2]>0)
      {
	size+=self_->ncells[2]*5;
      }
#endif
    
#if 0
        calculator_cells_to_cells_t<VolumeType::Hexahedron,FaceType::Quadrilateral,int_t>::calculate(zone->data.n,
												 zone->data.v,
												 zone->data.ld,
												 c2c,
												 6,
												 &numFaces,
												 &numFacesInterior,
												 &numFacesBoundary);
#endif
    
    
    //template <VolumeType::enum_t _volumeType,FaceType::enum_t _faceType,typename _int_t>
    //class 

  }
  
  
};




#if 1

template <typename _derivedClass>  void space(const CRTP_MeshTopology<_derivedClass> &meshTopology_)
{
  using this_t = CRTP_MeshTopology<_derivedClass>;
  using entitykind_t = typename this_t::entitykind_t;
  
  int_t begin[entitykind_t::NumKinds+1];

  begin[0] = 0;

  int_t size = 0;
  
  static constexpr DimensionType::enum_t meshDimension = this_t::Dimension;
  for( const auto cellType : entitykind_t::All)
    {
      const auto numCellsOfType = meshTopology_.template GetNumEntities<meshDimension>(cellType);
      size += numCellsOfType *  entitykind_t::GetNumNodes(cellType);
      begin[cellType+1] = numCellsOfType*  entitykind_t::GetNumNodes(cellType);
    }
  for (unsigned int i=1;i<=entitykind_t::NumKinds;++i)
    {
      size +=begin[i];
      //      std::cout << "[" << i << "]" << begin[i] << std::endl;
    }
  
  for (unsigned int i=2;i<=entitykind_t::NumKinds;++i)
    {
      begin[i] += begin[i-1];
      //      std::cout << "[" << i << "]" << begin[i] << std::endl;
    }
  for (unsigned int i=0;i<entitykind_t::NumKinds+1;++i)
    {
      std::cout << "[" << i << "]" << begin[i] << std::endl;
    }

  int_t* c2n = (int_t*)malloc(sizeof(int_t)*size);
  int_t* c2n_m = (int_t*)calloc(entitykind_t::NumKinds,sizeof(int_t));
  int_t* c2n_n = (int_t*)calloc(entitykind_t::NumKinds,sizeof(int_t));
  int_t* c2n_ld = (int_t*)calloc(entitykind_t::NumKinds,sizeof(int_t));
  {
    int_t*pc2n = c2n;
    for( const auto cellType : entitykind_t::All)
      {
	const auto numCellsOfType = meshTopology_.template GetNumEntities<meshDimension>(cellType);
	if ( numCellsOfType > 0 )
	  {
	  
	    typename this_t::template entity_t<DimensionType::Node> cellToNodes[8];
	    for (const auto cell : meshTopology_.template GetEntities<meshDimension>(cellType) )
	      {
		const auto numNodesInCell = meshTopology_.template GetEntityToEntities<meshDimension,DimensionType::Node>(cell,
															  cellToNodes);
	      
		for (unsigned int i=0;i<numNodesInCell;++i)
		  {
		    auto nodeIndex = meshTopology_.template GetEntityIndex<DimensionType::Node>(cellToNodes[i]);
		    pc2n[i] = nodeIndex;
		    //		  out_ << " " << nodeIndex + 1;
		  }
		pc2n += numNodesInCell;
		//	      out_ << " 0" << std::endl;
		//		out_ << cell << std::endl;
	      }
	    c2n_m[cellType] = entitykind_t::GetNumNodes(cellType);
	    c2n_n[cellType] = numCellsOfType;
	    c2n_ld[cellType] = c2n_m[cellType];
	  }
      }
  }

  //  wspace_t self;
  ///  wspace_def(&self,2,entitykind_t::NumKinds,begin,c2n);
  wspace_sparse_c2n_t wspace_sparse_c2n_linear;
  wspace_sparse_c2n_def(&wspace_sparse_c2n_linear,
		 entitykind_t::NumKinds,
		 begin,
		 c2n_m,
		 c2n_n,
		 c2n,
		 c2n_ld);
  
  wspace_sparse_c2n_info(&wspace_sparse_c2n_linear,stdout);
  wspace_sparse_c2n_calculate(&wspace_sparse_c2n_linear);
  //  wspace_sparse_c2n_calculate_c2c(&wspace_sparse_c2n_linear);
#if 0  
  return out_;
#endif
};
#endif


namespace Output
{
template <typename impl>
void info(const CRTP_MeshTopology<impl>& topology,FILE*out_)
{
  //  fprintf(out_,"num nodes : %Ld\n",topology.GetNumNodes());
  //  fprintf(out_,"num cells : %Ld\n",topology.GetNumCells());
    
    using this_t = AHF::MeshTopology3D;
    using entitykind_t = typename this_t::entitykind_t;

    static constexpr DimensionType::enum_t meshDimension = this_t::Dimension;
    fprintf(out_,"dimension %d\n",meshDimension);

    for( const auto cellType : entitykind_t::All)
      {
	const auto numCellsOfType = topology.template GetNumEntities<meshDimension>(cellType);
	if ( numCellsOfType > 0 )
	  {
	    std::cout << cellType
		      << std::endl
		      << numCellsOfType
		      << std::endl;
#if 0
	    typename this_t::template entity_t<DimensionType::Node> cellToNodes[8];
	    for (const auto cell : topology.template GetEntities<meshDimension>(cellType) )
	      {		
		const auto numNodesInCell = topology.template GetEntityToEntities<meshDimension,DimensionType::Node>(cell,
															  cellToNodes);

		for (unsigned int i=0;i<numNodesInCell;++i)
		  {
		    auto nodeIndex = topology.template GetEntityIndex<DimensionType::Node>(cellToNodes[i]);
		    std::cout << " " << nodeIndex + 1;
		  }
		std::cout << " 0" << std::endl;
		//		out_ << cell << std::endl;
	      }
#endif	    
	  }	
      }


    
    

}
}
#include "wmesh_functions.h"
extern "C"
{

  wmesh_integer_t wmesh_info	(const wmesh_t*		self_,
				 FILE * out_)
  {
    return WMESH_STATUS_SUCCESS;
  };

  wmesh_integer_t wmesh_kill(wmesh_t* self_)
  {
    if (self_)
      {
	free(self_);
      }
    return 0;  
  };

  wmesh_integer_t wmesh_build_s_e2n(wmesh_t* self_,int_t dim_)
  {
    if (dim_==3)
      {
	int_t s_e2n_m[4] {2,2,2,2};
	int_t s_e2n_ld[4] {2,2,2,2};
	int_t s_e2n_n[4]{6,8,9,12};
	int_t s_e2n_ptr[4+1]{0,12,28,46,70};
	int_t v[70]{
	  // tet
	  1,2,	  
	    2,0,
	    0,1,
	    2,3,
	    0,3,
	    1,3,
	    // pyr
	    0,1,
	    1,2,
	    2,3,
	    3,0,
	    0,4,
	    1,4,
	    2,4,
	    3,4,
	    //wedge
	    0,1,
	    1,2,
	    2,0,
	    3,4,
	    4,5,
	    5,3,
	    0,3,
	    1,4,
	    2,5,
	    // hex
	    0,1,
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
	    3,7
	    };
      
	int_t * __restrict__ s_e2n_v = (int_t * __restrict__)malloc(sizeof(int_t)*s_e2n_ptr[4]);
	memcpy(s_e2n_v,v,sizeof(int_t)*s_e2n_ptr[4]);
	wint_sparsemat_def(&self_->m_s_e2n,
			   4);
	wint_sparsemat_set(&self_->m_s_e2n,
			   s_e2n_m,
			   s_e2n_n,
			   s_e2n_ld,
			   s_e2n_ptr,
			   s_e2n_v);
      }
    else if (dim_==2)
      {
	int_t s_e2n_m[2] {2,2};
	int_t s_e2n_ld[2] {2,2};
	int_t s_e2n_n[2]{3,4};
	int_t s_e2n_ptr[2+1]{0,6,14};
	int_t v[14]{ 0,1, 1,2, 2,0, 0,1, 1,2, 2,3, 3,0};
	int_t * __restrict__ s_e2n_v = (int_t * __restrict__)malloc(sizeof(int_t)*s_e2n_ptr[2]);
	memcpy(v,s_e2n_v,sizeof(int_t)*s_e2n_ptr[2]);
	wint_sparsemat_def(&self_->m_s_e2n,
			   2);
	wint_sparsemat_set(&self_->m_s_e2n,
			   s_e2n_m,
			   s_e2n_n,
			   s_e2n_ld,
			   s_e2n_ptr,
			   s_e2n_v);	
      }
    return WMESH_STATUS_SUCCESS;
  };
  
  wmesh_integer_t wmesh_analysis(wmesh_t* 		self_)
  {
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    WMESH_STATUS_CHECK(wmesh_build_s_e2n(self_,3));
    //
    // Create c2e    
    //
    {
      int_t c2e_m[4] {6,8,9,12};
      int_t c2e_ld[4] {6,8,9,12};
      int_t c2e_ptr[4+1];
      int_t c2e_n[4];
      
      for (int_t i=0;i<4;++i) c2e_n[i] = self_->m_c2n.m_data.m_n[i];
      c2e_ptr[0] = 0;
      c2e_ptr[1] = c2e_ptr[0] + c2e_n[0] * c2e_ld[0];
      c2e_ptr[2] = c2e_ptr[1] + c2e_n[1] * c2e_ld[1];
      c2e_ptr[3] = c2e_ptr[2] + c2e_n[2] * c2e_ld[2];
      c2e_ptr[4] = c2e_ptr[3] + c2e_n[3] * c2e_ld[3];
      int_t * __restrict__ c2e_v = (int_t * __restrict__)malloc(sizeof(int_t)*c2e_ptr[4]);
      wint_sparsemat_def(&self_->m_c2e,
			 4);
      wint_sparsemat_set(&self_->m_c2e,
			 c2e_m,
			 c2e_n,
			 c2e_ld,
			 c2e_ptr,
			 c2e_v);
    }

    int_t edge_idx = 0;
    for (int j=0;j<4;++j)
      {
	if (self_->m_c2n.m_data.m_n[j]>0)
	  {
	    edge_idx = 0;
	    int_t i = j;
#if 0
	    int_t s_e2n_m = 2;
	    int_t s_e2n_n = self_->m_s_e2n.m_m[j];
	    const unsigned int * __restrict__ s_e2n_v = VolumeType::EdgesToNodes[i];
	    int_t s_e2n_ld = s_e2n_m;
#endif

#define hea(_f,_i)							\
	    _f.m_m[_i],							\
	      _f.m_n[_i],					\
	      _f.m_data + _f.m_ptr[_i],				\
	      _f.m_ld[_i]
	    
	    int_t work_n;
	    int_t * work  = nullptr;
	    int_t status = wmesh_indexing_edges(self_->m_c2n.m_data.m_m[i],
						self_->m_c2n.m_data.m_n[i],
						&self_->m_c2n.m_data.m_data[self_->m_c2n.m_data.m_ptr[i]],
						self_->m_c2n.m_data.m_ld[i],
												
						self_->m_c2e.m_m[i],
						self_->m_c2e.m_n[i],
						self_->m_c2e.m_data + self_->m_c2e.m_ptr[i],
						self_->m_c2e.m_ld[i],
						
						self_->m_s_e2n.m_m[i],
						self_->m_s_e2n.m_n[i],
						self_->m_s_e2n.m_data + self_->m_s_e2n.m_ptr[i],
						self_->m_s_e2n.m_ld[i],
						
						&edge_idx,
						&work_n,
						work);
	    WMESH_STATUS_CHECK(status);
	  
	  work =  (int_t*)malloc(work_n*sizeof(int_t));
	  work[0] = self_->m_c2n.m_data.m_n[i];
	  status = wmesh_indexing_edges(self_->m_c2n.m_data.m_m[i],
					self_->m_c2n.m_data.m_n[i],
					&self_->m_c2n.m_data.m_data[self_->m_c2n.m_data.m_ptr[i]],
					self_->m_c2n.m_data.m_ld[i],
					
					self_->m_c2e.m_m[i],
					self_->m_c2e.m_n[i],
					self_->m_c2e.m_data + self_->m_c2e.m_ptr[i],
					self_->m_c2e.m_ld[i],
					
					self_->m_s_e2n.m_m[i],
					self_->m_s_e2n.m_n[i],
					self_->m_s_e2n.m_data + self_->m_s_e2n.m_ptr[i],
					self_->m_s_e2n.m_ld[i],

//					s_e2n_m,
//					s_e2n_n,
//					s_e2n_v,
//					s_e2n_ld,
					
					&edge_idx,
					&work_n,
					work);
	  
	  WMESH_STATUS_CHECK(status);
	  free(work);
	}
      }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

    std::cout << "wmesh_analysis: " << time_span.count() << " seconds.";
    std::cout << std::endl;

    return WMESH_STATUS_SUCCESS;
  }

#if 0
  wmesh_integer_t wmesh_analysis_old(wmesh_t* 		self_)
  {
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    int_t c2e_m[4] {6,8,9,12};
    int_t c2e_n[4];
    int_t c2e_ld[4] {6,8,9,12};
    int_t c2e_ptr[4+1];
    
    for (int_t i=0;i<4;++i) c2e_n[i] = self_->m_c2n.m_data.m_n[i];

    c2e_ptr[0] = 0;
    c2e_ptr[1] = c2e_ptr[0] + c2e_n[0] * c2e_ld[0];
    c2e_ptr[2] = c2e_ptr[1] + c2e_n[1] * c2e_ld[1];
    c2e_ptr[3] = c2e_ptr[2] + c2e_n[2] * c2e_ld[2];
    c2e_ptr[4] = c2e_ptr[3] + c2e_n[3] * c2e_ld[3];
    
    int_t * __restrict__ c2e_v = (int_t * __restrict__)malloc(sizeof(int_t)*c2e_ptr[4]);

    wint_sparsemat_def(&self_->m_c2e,
		       4);
    wint_sparsemat_set(&self_->m_c2e,
		       c2e_m,
		       c2e_n,
		       c2e_ld,
		       c2e_ptr,
		       c2e_v);


    int_t edge_idx = 0;
    for (int j=0;j<4;++j)
      if (self_->m_c2n.m_data.m_n[j]>0)
	{
	  edge_idx = 0;
	  int_t i = j;
	  int_t s_e2n_m = 2;
	  int_t s_e2n_n = c2e_m[j];
	  const unsigned int * __restrict__ s_e2n_v = VolumeType::EdgesToNodes[i];
	  int_t s_e2n_ld = s_e2n_m;
	  
	  int_t work_n;
	  int_t * work  = nullptr;
	  int_t status = wmesh_indexing_edges(self_->m_c2n.m_data.m_m[i],
					      self_->m_c2n.m_data.m_n[i],
					      &self_->m_c2n.m_data.m_data[self_->m_c2n.m_data.m_ptr[i]],
					      self_->m_c2n.m_data.m_ld[i],
					      
					      c2e_m[i],
					      c2e_n[i],
					      &c2e_v[c2e_ptr[i]],
					      c2e_ld[i],
					      
					      s_e2n_m,
					      s_e2n_n,
					      s_e2n_v,
					      s_e2n_ld,
					      
					      &edge_idx,
					      &work_n,
					      work);
	  
	  work =  (int_t*)malloc(work_n*sizeof(int_t));
	  work[0] = self_->m_c2n.m_data.m_n[i];
	  status = wmesh_indexing_edges(self_->m_c2n.m_data.m_m[i],
					self_->m_c2n.m_data.m_n[i],
					&self_->m_c2n.m_data.m_data[self_->m_c2n.m_data.m_ptr[i]],
					self_->m_c2n.m_data.m_ld[i],
					
					c2e_m[i],
					c2e_n[i],
					&c2e_v[c2e_ptr[i]],
					c2e_ld[i],
					
					s_e2n_m,
					s_e2n_n,
					s_e2n_v,
					s_e2n_ld,
					
					&edge_idx,
					&work_n,
					work);
	  
	  WMESH_STATUS_CHECK(status);
	  free(work);
	}



    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

    std::cout << "wmesh_analysis: " << time_span.count() << " seconds.";
    std::cout << std::endl;

    return WMESH_STATUS_SUCCESS;
  }
#endif

  wmesh_integer_t wmesh_write(const wmesh_t * 	self_,
			      const char * 	filename_)
  {
#if 0
    int i;
    for (i=0;filename_[i]!='\0';++i);
    if ( (i>=5) && !strcmp(".mesh",&filename_[i-5]))
      {
	printf("helllo\n");
	Output::Medit outputMedit(filename_);
	outputMedit << *self_->mesh;
	return WMESH_STATUS_SUCCESS;
      }
    else if ( (i>=4) && !strcmp(".vtk",&filename_[i-4]))
      {
	Output::Vtk::Writer outputVtk(filename_);
	outputVtk << *self_->mesh;
	return WMESH_STATUS_SUCCESS;
      }
    else
      {
	return WMESH_STATUS_INVALID_ARGUMENT;
      }
#endif
    return 0;
  }


  wmesh_integer_t wmesh_refine_uniform(const wmesh_t * 	self_,
				      wmesh_t ** 		refined_mesh_)
  {
    refined_mesh_[0] = (wmesh_t*)calloc(1,sizeof(wmesh_t));
    
    //
    //
    //
    



    
    std::cout << "hey" << std::endl;
    return 0;
  }


  //!
  //! @brief Get the number of entities for a specific dimension.
  //!
  wmesh_integer_t wmesh_num_entities_from_dimension	(const wmesh_t* 	self_,
							 wmesh_integer_t		dimension_,
							 wmesh_integer_p		out_num_entities_)
  {
#if 0
    if (dimension_ == 0)
      {
	out_num_entities_[0] = self_->mesh->GetTopology()->GetNumNodes();	
      }
    else if (dimension_ == 1)
      {
	out_num_entities_[0] = self_->mesh->GetTopology()->GetNumEntities<DimensionType::Edge>();	
      }
    else if (dimension_ == 2)
      {
	out_num_entities_[0] = self_->mesh->GetTopology()->GetNumEntities<DimensionType::Face>();	
      }
    else if (dimension_ == 3)
      {
	out_num_entities_[0] = self_->mesh->GetTopology()->GetNumEntities<DimensionType::Volume>();	
      }
    else
      {
	return WMESH_STATUS_INVALID_ENUM;
      }
#endif
    return WMESH_STATUS_SUCCESS;
  }

  //!
  //! @brief Get the number of entities for a specific dimension and a specific type.
  //!
  wmesh_integer_t wmesh_num_entities_from_element_type	(const wmesh_t 	*self_,
							 wmesh_integer_t 		element_,
							 wmesh_integer_p		out_num_entities_)
  {

#if 0
    out_num_entities_[0] = self_->mesh->GetTopology()->template GetNumEntities<DimensionType::Edge>(EdgeType::Edge);
    out_num_entities_[0] = self_->mesh->GetTopology()->template GetNumEntities<DimensionType::Face>(FaceType::Triangle);
    out_num_entities_[0] = self_->mesh->GetTopology()->template GetNumEntities<DimensionType::Face>(FaceType::Quadrilateral);
#endif
#if 0
    if (element_ == (wmesh_integer_t)VolumeType::Wedge)
      out_num_entities_[0] = self_->mesh->GetTopology()->template GetNumEntities<DimensionType::Volume>(VolumeType::Wedge);
    else if (element_ == (wmesh_integer_t)VolumeType::Pyramid)
      out_num_entities_[0] = self_->mesh->GetTopology()->template GetNumEntities<DimensionType::Volume>(VolumeType::Pyramid);
    else if (element_ == (wmesh_integer_t)VolumeType::Tetrahedron)
      out_num_entities_[0] = self_->mesh->GetTopology()->template GetNumEntities<DimensionType::Volume>(VolumeType::Tetrahedron);
    else if (element_ == (wmesh_integer_t)VolumeType::Hexahedron)
      out_num_entities_[0] = self_->mesh->GetTopology()->template GetNumEntities<DimensionType::Volume>(VolumeType::Hexahedron);
    else
      return WMESH_STATUS_INVALID_ENUM;
#endif
    return WMESH_STATUS_SUCCESS;
  }


  

  //!
  //! @brief Get the number of entities for a specific dimension.
  //!
  wmesh_integer_t wmesh_entity_index			(const wmesh_t* 	self_,
							 wmesh_integer_t 	entity,
							 wmesh_integer_p	out_entity_index)
  {

    return 0;
  }

  
  //!
  //! @brief Get the nodes of a specific entity of a cell.
  //!
  wmesh_integer_t wmesh_entity2nodes			(const wmesh_t* 	self_,
							 wmesh_integer_t	entity_,
							 wmesh_integer_t	localEntityDimension_,
							 wmesh_integer_t	localEntityIndex_,
							 wmesh_integer_t	entityToNodes_[])
  {

    return 0;
  }

  wmesh_integer_t wmesh_entity2entities			(const wmesh_t* 	self_,
									 wmesh_integer_t		source_dimension_,
									 wmesh_integer_t 		source_entity_,
									 wmesh_integer_t		target_dimension_,
									 wmesh_integer_t 		target_entities_[])
  {

    return 0;
  }

};
#endif
