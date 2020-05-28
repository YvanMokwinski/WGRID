
#pragma once

template<wmesh_int_t 	degree_,
	 wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_template_shape_orthogonal_splz
{  
  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_);
};


template<typename T>
struct bms_template_shape_orthogonal_splz<0,WMESH_ELEMENT_EDGE,T>
{

  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T zero  = static_cast<T>(0);
    static constexpr T one  = static_cast<T>(1);
    static constexpr T two  = static_cast<T>(2);
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = zero;
      }
    else
      {
	b_[b_inc_*0] = one / wmesh_math<T>::xsqrt(two);
      }
  }
  
};

template<typename T>
struct bms_template_shape_orthogonal_splz<0,WMESH_ELEMENT_TRIANGLE,T>
{

  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T zero  = static_cast<T>(0);
    static constexpr T two  = static_cast<T>(2);
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = zero;
      }
    else
      {
	b_[b_inc_*0] = wmesh_math<T>::xsqrt(two);
      }
  }
  
};

template<typename T>
struct bms_template_shape_orthogonal_splz<0,WMESH_ELEMENT_QUADRILATERAL,T>
{

  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T zero  = static_cast<T>(0);
    static constexpr T one  = static_cast<T>(1);
    static constexpr T two  = static_cast<T>(2);
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = zero;
      }
    else
      {
	b_[b_inc_*0] = one / two;
      }
  }
  
};

template<typename T>
struct bms_template_shape_orthogonal_splz<0,WMESH_ELEMENT_TETRAHEDRON,T>
{

  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T zero  = static_cast<T>(0);
    static constexpr T six  = static_cast<T>(6);
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = zero;
      }
    else
      {
	b_[b_inc_*0] = wmesh_math<T>::xsqrt(six);
      }
  }
  
};


template<typename T>
struct bms_template_shape_orthogonal_splz<0,WMESH_ELEMENT_PYRAMID,T>
{

  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T zero    = static_cast<T>(0);
    static constexpr T two    = static_cast<T>(2);
    static constexpr T three  = static_cast<T>(3);
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = zero;
      }
    else
      {
	b_[b_inc_*0] = wmesh_math<T>::xsqrt(three)/two;
      }
  }
  
};

template <typename T>
struct bms_template_shape_orthogonal_splz<0,WMESH_ELEMENT_WEDGE,T>
{

  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T zero  = static_cast<T>(0);
    static constexpr T two  = static_cast<T>(2);
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = zero;
      }
    else
      {
	b_[b_inc_*0] = wmesh_math<T>::xsqrt(two);
      }
  }
  
  
};

template <typename T>
struct bms_template_shape_orthogonal_splz<0,WMESH_ELEMENT_HEXAHEDRON,T>
{


  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T zero  = static_cast<T>(0);
    static constexpr T one  = static_cast<T>(1);
    static constexpr T eight  = static_cast<T>(8);
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = zero;
      }
    else
      {
	b_[b_inc_*0] = one / wmesh_math<T>::xsqrt(eight);
      }
  }
  
  
};


