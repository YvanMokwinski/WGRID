#pragma once
#include "wmesh-enums.hpp"

struct wmesh_enum_element_t
{
public:
  typedef enum : wmesh_int_t    
  {
    NODE = WMESH_ELEMENT_NODE,
      EDGE = WMESH_ELEMENT_EDGE,
      TRIANGLE = WMESH_ELEMENT_TRIANGLE,
      QUADRILATERAL = WMESH_ELEMENT_QUADRILATERAL,
      TETRAHEDRON =WMESH_ELEMENT_TETRAHEDRON,
      PYRAMID=WMESH_ELEMENT_PYRAMID,
      WEDGE=WMESH_ELEMENT_WEDGE,
      HEXAHEDRON=WMESH_ELEMENT_HEXAHEDRON
      } kind_t;

private: kind_t m_kind{};
  
public: static constexpr kind_t all[8] = {    NODE,
					      EDGE,
					      TRIANGLE,
					      QUADRILATERAL,
					      TETRAHEDRON,
					      PYRAMID,
					      WEDGE,
					      HEXAHEDRON};
  
public: inline static bool is_valid(kind_t kind_) 
  {
    return
      (NODE == kind_) ||
      (EDGE == kind_) ||
      (TRIANGLE == kind_) ||
      (QUADRILATERAL == kind_) ||
      (TETRAHEDRON == kind_) ||
      (PYRAMID == kind_) ||
      (WEDGE == kind_) ||
      (HEXAHEDRON == kind_);
  };
  
public: inline static bool is_invalid(kind_t kind_) 
  {
    return !is_valid(kind_);
  };
  
public: inline wmesh_enum_element_t(kind_t kind_) noexcept : m_kind(kind_){std::cout << "constructor" << std::endl; };
public: inline wmesh_enum_element_t() noexcept = delete;
public: inline operator kind_t() const noexcept{ std::cout << "cast" << std::endl; return m_kind; };
  
};

constexpr wmesh_enum_element_t::kind_t wmesh_enum_element_t::all[8];




struct wmesh_enum_nodes_t
{
    static constexpr wmesh_int_t num_kinds = 2;
public:
  typedef enum : wmesh_int_t    
  {
    LAGRANGE = WMESH_NODES_FAMILY_LAGRANGE,
      GAUSSLOBATTO = WMESH_NODES_FAMILY_GAUSSLOBATTO
      } kind_t;
  
private: kind_t m_kind{};  
public: static constexpr kind_t all[num_kinds] = {    LAGRANGE,
					      GAUSSLOBATTO };
  
public: inline static bool is_valid(kind_t kind_) 
  {
    return
      (LAGRANGE == kind_) ||
      (GAUSSLOBATTO == kind_);
  };
  
public: inline static bool is_invalid(kind_t kind_) 
  {
    return !is_valid(kind_);
  };
  
public: inline wmesh_enum_nodes_t(kind_t kind_) noexcept : m_kind(kind_){};
public: inline wmesh_enum_nodes_t() noexcept = delete;
public: inline operator kind_t() const noexcept{ return m_kind; };
  
};

constexpr wmesh_enum_nodes_t::kind_t wmesh_enum_nodes_t::all[wmesh_enum_nodes_t::num_kinds];


struct wmesh_enum_shape_t
{
public:
  static constexpr wmesh_int_t num_kinds = 3;
  typedef enum : wmesh_int_t    
  {
    LAGRANGE = WMESH_SHAPE_FAMILY_LAGRANGE,
      LEGENDRE = WMESH_SHAPE_FAMILY_LEGENDRE,
      ORTHOGONAL = WMESH_SHAPE_FAMILY_ORTHOGONAL
      } kind_t;
  
private: kind_t m_kind{};  
public: static constexpr kind_t all[num_kinds] = {      LAGRANGE,
							LEGENDRE,
							ORTHOGONAL};
  
public: inline static bool is_valid(kind_t kind_) 
  {
    return
      (LAGRANGE == kind_) ||
      (LEGENDRE == kind_) ||
      (ORTHOGONAL == kind_);
  };
  
public: inline static bool is_invalid(kind_t kind_) 
  {
    return !is_valid(kind_);
  };
  
public: inline wmesh_enum_shape_t(kind_t kind_) noexcept : m_kind(kind_){};
public: inline wmesh_enum_shape_t() noexcept = delete;
public: inline operator kind_t() const noexcept{ return m_kind; };
  
};

constexpr wmesh_enum_shape_t::kind_t wmesh_enum_shape_t::all[wmesh_enum_shape_t::num_kinds];







#if 0
struct wmesh_enum_nodes_family
{
  typedef enum : wmesh_int_t    
  {
    LAGRANGE    = WMESH_NODES_FAMILY_LAGRANGE,
      GAUSSLOBATTO= WMESH_NODES_FAMILY_GAUSSLOBATTO } kind_t;
  wmesh_int_t num_kinds = 2;  
};
typedef enum : wmesh_int_t    
  {
    LAGRANGE    = WMESH_SHAPE_FAMILY_LAGRANGE,
      LEGENDRE  = WMESH_SHAPE_FAMILY_LEGENDRE,
      ORTHOGONAL = WMESH_SHAPE_FAMILY_ORTHOGONAL} kind_t;


struct wmesh_enum_shape_family
{
  typedef enum : wmesh_int_t    
  {
    LAGRANGE    = WMESH_SHAPE_FAMILY_LAGRANGE,
      LEGENDRE  = WMESH_SHAPE_FAMILY_LEGENDRE,
      ORTHOGONAL = WMESH_SHAPE_FAMILY_ORTHOGONAL} kind_t;
  wmesh_int_t num_kinds = 3;
};
#endif



#if 0
#define WMESH_NODES_FAMILY_LAGRANGE 		0
#define WMESH_NODES_FAMILY_GAUSSLOBATTO 	1
#define WMESH_NODES_FAMILY_ALL 			2

#define WMESH_SHAPE_FAMILY_LAGRANGE   		0
#define WMESH_SHAPE_FAMILY_LEGENDRE 		1
#define WMESH_SHAPE_FAMILY_ORTHOGONAL 		2
#define WMESH_SHAPE_FAMILY_ALL 			3
#endif
#if 0
struct wmesh_enum_element_t
{
public:
  typedef enum : wmesh_int_t    
  {
    NODE = WMESH_ELEMENT_NODE,
      EDGE = WMESH_ELEMENT_EDGE,
      TRIANGLE = WMESH_ELEMENT_TRIANGLE,
      QUADRILATERAL = WMESH_ELEMENT_QUADRILATERAL,
      TETRAHEDRON =WMESH_ELEMENT_TETRHAEDRON,
      PYRAMID=WMESH_ELEMENT_PYRAMID,
      WEDGE=WMESH_ELEMENT_WEDGE,
      HEXAHEDRON=WMESH_ELEMENT_HEXAHEDRON
      } kind_t;

private: kind_t m_kind{};
  
public: static constexpr kind_t all[8] = {    NODE,
					      EDGE,
					      TRIANGLE,
					      QUADRILATERAL,
					      TETRAHEDRON,
					      PYRAMID,
					      WEDGE,
					      HEXAHEDRON};
  
public: inline static bool is_valid(kind_t kind_) 
  {
    return
      (NODE == kind_) ||
      (EDGE == kind_) ||
      (TRIANGLE == kind_) ||
      (QUADRILATERAL == kind_) ||
      (TETRAHEDRON == kind_) ||
      (PYRAMID == kind_) ||
      (WEDGE == kind_) ||
      (HEXAHEDRON == kind_);
    return kind_ != block && kind_ != interleave;
  };
  
public: inline static bool is_invalid(kind_t kind_) 
  {
    return !is_valid(kind_);
  };
  
public: inline wmesh_enum_element_t(kind_t kind_) noexcept;
public: inline wmesh_enum_element_t() noexcept;
public: inline operator kind_t() const noexcept;
  
};
#endif
