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
};

template<>
struct bms_traits_element<WMESH_ELEMENT_EDGE>
{
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_EDGE;
};

template<>
struct bms_traits_element<WMESH_ELEMENT_TRIANGLE>
{
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_FACE;
};

template<>
struct bms_traits_element<WMESH_ELEMENT_QUADRILATERAL>
{
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_FACE;
};

template<>
struct bms_traits_element<WMESH_ELEMENT_TETRAHEDRON>
{
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_VOLUME;
};

template<>
struct bms_traits_element<WMESH_ELEMENT_PYRAMID>
{
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_VOLUME;
};

template<>
struct bms_traits_element<WMESH_ELEMENT_WEDGE>
{
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_VOLUME;
};

template<>
struct bms_traits_element<WMESH_ELEMENT_HEXAHEDRON>
{
  static constexpr wmesh_int_t s_topodim = WMESH_TOPODIM_VOLUME;
};
