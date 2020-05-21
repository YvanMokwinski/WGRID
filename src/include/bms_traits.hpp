#pragma once



template<wmesh_int_t TOPODIM>
struct bms_traits_topodim;

template<>
struct bms_traits_topodim<WMESH_TOPODIM_VOLUME>
{
  static constexpr wmesh_int_t s_ntypes = 4;
};

template<>
struct bms_traits_topodim<WMESH_TOPODIM_FACE>
{
  static constexpr wmesh_int_t s_ntypes = 2;
};

template<>
struct bms_traits_topodim<WMESH_TOPODIM_EDGE>
{
  static constexpr wmesh_int_t s_ntypes = 1;
};

template<>
struct bms_traits_topodim<WMESH_TOPODIM_NODE>
{
  static constexpr wmesh_int_t s_ntypes = 1;
};





template<wmesh_int_t ELEMENT>
struct bms_traits_element;

template<>
struct bms_traits_element<WMESH_ELEMENT_NODE>
{
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_NODE;
  static constexpr wmesh_int_t s_num_nodes = 1;
  static constexpr wmesh_int_t s_num_edges = 0;
  static constexpr wmesh_int_t s_num_triangles = 0;
  static constexpr wmesh_int_t s_num_quadrilaterals = 0;
};

template<>
struct bms_traits_element<WMESH_ELEMENT_EDGE>
{
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_EDGE;
  static constexpr wmesh_int_t s_num_nodes = 2;
  static constexpr wmesh_int_t s_num_edges = 1;
  static constexpr wmesh_int_t s_num_triangles = 0;
  static constexpr wmesh_int_t s_num_quadrilaterals = 0;
};

template<>
struct bms_traits_element<WMESH_ELEMENT_TRIANGLE>
{
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_FACE;
  static constexpr wmesh_int_t s_num_nodes = 3;
  static constexpr wmesh_int_t s_num_edges = 3;
  static constexpr wmesh_int_t s_num_triangles = 1;
  static constexpr wmesh_int_t s_num_quadrilaterals = 0;
};

template<>
struct bms_traits_element<WMESH_ELEMENT_QUADRILATERAL>
{
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_FACE;
  static constexpr wmesh_int_t s_num_nodes = 4;
  static constexpr wmesh_int_t s_num_edges = 4;
  static constexpr wmesh_int_t s_num_triangles = 0;
  static constexpr wmesh_int_t s_num_quadrilaterals = 1;

};

template<>
struct bms_traits_element<WMESH_ELEMENT_TETRAHEDRON>
{
  static constexpr wmesh_int_t s_num_nodes = 4;
  static constexpr wmesh_int_t s_num_edges = 6;
  static constexpr wmesh_int_t s_num_triangles = 4;
  static constexpr wmesh_int_t s_num_quadrilaterals = 0;
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_VOLUME;
};

template<>
struct bms_traits_element<WMESH_ELEMENT_PYRAMID>
{
  static constexpr wmesh_int_t s_num_nodes = 5;
  static constexpr wmesh_int_t s_num_edges = 8;
  static constexpr wmesh_int_t s_num_triangles = 4;
  static constexpr wmesh_int_t s_num_quadrilaterals = 1;
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_VOLUME;
};

template<>
struct bms_traits_element<WMESH_ELEMENT_WEDGE>
{
  static constexpr wmesh_int_t s_num_nodes = 6;
  static constexpr wmesh_int_t s_num_edges = 9;
  static constexpr wmesh_int_t s_num_triangles = 2;
  static constexpr wmesh_int_t s_num_quadrilaterals = 3;
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_VOLUME;
};

template<>
struct bms_traits_element<WMESH_ELEMENT_HEXAHEDRON>
{
  static constexpr wmesh_int_t s_num_nodes = 8;
  static constexpr wmesh_int_t s_num_edges = 12;
  static constexpr wmesh_int_t s_num_triangles = 0;
  static constexpr wmesh_int_t s_num_quadrilaterals = 6;
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_VOLUME;
};
