
#pragma once

template<wmesh_int_t 	degree_,
	 wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_template_shape_legendre_splz
{  
  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_);
};


template<wmesh_int_t ELEMENT_,typename T>
struct bms_template_shape_legendre_splz<0,ELEMENT_,T>
{

  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T one  = static_cast<T>(1);
    static constexpr T zero = static_cast<T>(0);
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = zero;
      }
    else
      {
	b_[b_inc_*0] = one;
      }
  }
  
};

template<typename T>
struct bms_template_shape_legendre_splz<1,WMESH_ELEMENT_EDGE,T>
{
  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T zero = static_cast<T>(0);
    static constexpr T one = static_cast<T>(1);
    const T
      r = lcoo_[lcoo_inc_* 0];
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = one;
      }
    else
      {
	b_[b_inc_*0] = one;
	b_[b_inc_*1] = r;
      }
  }
};

template<typename T>
struct bms_template_shape_legendre_splz<1,WMESH_ELEMENT_TRIANGLE,T>
{

  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T one = static_cast<T>(1);
    static constexpr T zero = static_cast<T>(0);
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = -one;
	b_[b_inc_*1] = one;
	b_[b_inc_*2] = zero;	
      }
    else if (diff_[1] > 0)
      {
	b_[b_inc_*0] = -one;
	b_[b_inc_*1] = zero;
	b_[b_inc_*2] = one;	
      }
    else
      {
	b_[b_inc_*0] = one - (lcoo_[lcoo_inc_ * 0] + lcoo_[lcoo_inc_ * 1]);
	b_[b_inc_*1] = lcoo_[lcoo_inc_ * 0];
	b_[b_inc_*2] = lcoo_[lcoo_inc_ * 1];
      }
  }
};

template<typename T>
struct bms_template_shape_legendre_splz<1,WMESH_ELEMENT_QUADRILATERAL,T>
{

  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T one = static_cast<T>(1);
    static constexpr T four = static_cast<T>(4);
    const T
      r = lcoo_[lcoo_inc_* 0],
      s = lcoo_[lcoo_inc_* 1];
    
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = -(one - s) / four;
	b_[b_inc_*1] = (one - s) / four;
	b_[b_inc_*2] = (one + s) / four;
	b_[b_inc_*3] = -(one + s) / four;
      }
    else if (diff_[1] > 0)
      {
	b_[b_inc_*0] = -(one - r) / four;
	b_[b_inc_*1] = -(one + r) / four;
	b_[b_inc_*2] = (one + r) / four;
	b_[b_inc_*3] = (one - r) / four;
      }
    else
      {
	b_[b_inc_*0] = (one - r)*(one - s) / four;
	b_[b_inc_*1] = (one + r)*(one - s) / four;
	b_[b_inc_*2] = (one + r)*(one + s) / four;
	b_[b_inc_*3] = (one - r)*(one + s) / four;
      }
  }


};

template<typename T>
struct bms_template_shape_legendre_splz<1,WMESH_ELEMENT_TETRAHEDRON,T>
{

  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T one = static_cast<T>(1);
    static constexpr T zero = static_cast<T>(0);
    const T
      r = lcoo_[lcoo_inc_*0],
      s = lcoo_[lcoo_inc_*1],
      t = lcoo_[lcoo_inc_*2];
    
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = -one;
	b_[b_inc_*1] = one;
	b_[b_inc_*2] = zero; 
	b_[b_inc_*3] = zero;
      }
    else if (diff_[1] > 0)
      {
	b_[b_inc_*0] = -one;
	b_[b_inc_*1] = zero;
	b_[b_inc_*2] = one;
	b_[b_inc_*3] = zero;
      }
    else if (diff_[2] > 0)
      {
	b_[b_inc_*0] = -one;
	b_[b_inc_*1] = zero;
	b_[b_inc_*2] = zero;
	b_[b_inc_*3] = one;
      }
    else
      {
	b_[b_inc_*0] = one - (r + s + t);
	b_[b_inc_*1] = r;
	b_[b_inc_*2] = s;
	b_[b_inc_*3] = t;
      }
  }



};

template<typename T>
struct bms_template_shape_legendre_splz<1,WMESH_ELEMENT_WEDGE,T>
{

  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T one = static_cast<T>(1);
    static constexpr T zero = static_cast<T>(0);
    const T
      r = lcoo_[lcoo_inc_*0],
      s = lcoo_[lcoo_inc_*1],
      t = lcoo_[lcoo_inc_*2];
    
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = -(one - t);
	b_[b_inc_*1] = (one - t);
	b_[b_inc_*2] = zero;
	b_[b_inc_*3] = -t;
	b_[b_inc_*4] = t;
	b_[b_inc_*5] = zero;
      }
    else if (diff_[1] > 0)
      {
	b_[b_inc_*0] = -(one - t);
	b_[b_inc_*1] = zero;
	b_[b_inc_*2] = (one - t);
	b_[b_inc_*3] = -t;
	b_[b_inc_*4] = zero;
	b_[b_inc_*5] = t;
      }
    else if (diff_[2] > 0)
      {
	b_[b_inc_*0] = -( one - (r + s) );
	b_[b_inc_*1] = -r;
	b_[b_inc_*2] = -s;
	b_[b_inc_*3] = ( one - (r + s) );
	b_[b_inc_*4] = r;
	b_[b_inc_*5] = s;
      }
    else
      {
	b_[b_inc_*0] = ( one - (r + s) ) * (one - t);
	b_[b_inc_*1] = r * (one - t);
	b_[b_inc_*2] = s * (one - t);
	b_[b_inc_*3] = ( one - (r + s) ) * t;
	b_[b_inc_*4] = r * t;
	b_[b_inc_*5] = s * t;
      }
  }


};

template<typename T>
struct bms_template_shape_legendre_splz<1,WMESH_ELEMENT_HEXAHEDRON,T>
{

  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T one = static_cast<T>(1);
    static constexpr T eight = static_cast<T>(8);
    const T
      r = lcoo_[lcoo_inc_*0],
      s = lcoo_[lcoo_inc_*1],
      t = lcoo_[lcoo_inc_*2];
    
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = -(one - s)*(one - t) / eight;
	b_[b_inc_*1] = (one - s)*(one - t) / eight;
	b_[b_inc_*2] = (one + s)*(one - t) / eight;
	b_[b_inc_*3] = -(one + s)*(one - t) / eight;
	b_[b_inc_*4] = -(one - s)*(one + t) / eight;
	b_[b_inc_*5] = (one - s)*(one + t) / eight;
	b_[b_inc_*6] = (one + s)*(one + t) / eight;
	b_[b_inc_*7] = -(one + s)*(one + t) / eight;
      }
    else if (diff_[1] > 0)
      {
	b_[b_inc_*0] = -(one - r)*(one - t) / eight;
	b_[b_inc_*1] = -(one + r)*(one - t) / eight;
	b_[b_inc_*2] =  (one + r)*(one - t) / eight;
	b_[b_inc_*3] =  (one - r)*(one - t) / eight;
	b_[b_inc_*4] = -(one - r)*(one + t) / eight;
	b_[b_inc_*5] = -(one + r)*(one + t) / eight;
	b_[b_inc_*6] =  (one + r)*(one + t) / eight;
	b_[b_inc_*7] =  (one - r)*(one + t) / eight;
      }
    else if (diff_[2] > 0)
      {
	b_[b_inc_*0] = -(one - r)*(one - s) / eight;
	b_[b_inc_*1] = -(one + r)*(one - s) / eight;
	b_[b_inc_*2] = -(one + r)*(one + s) / eight;
	b_[b_inc_*3] = -(one - r)*(one + s) / eight;
	b_[b_inc_*4] =  (one - r)*(one - s) / eight;
	b_[b_inc_*5] =  (one + r)*(one - s) / eight;
	b_[b_inc_*6] =  (one + r)*(one + s) / eight;
	b_[b_inc_*7] =  (one - r)*(one + s) / eight;
      }
    else
      {
	b_[b_inc_*0] = (one - r)*(one - s)*(one - t) / eight;
	b_[b_inc_*1] = (one + r)*(one - s)*(one - t) / eight;
	b_[b_inc_*2] = (one + r)*(one + s)*(one - t) / eight;
	b_[b_inc_*3] = (one - r)*(one + s)*(one - t) / eight;
	b_[b_inc_*4] = (one - r)*(one - s)*(one + t) / eight;
	b_[b_inc_*5] = (one + r)*(one - s)*(one + t) / eight;
	b_[b_inc_*6] = (one + r)*(one + s)*(one + t) / eight;
	b_[b_inc_*7] = (one - r)*(one + s)*(one + t) / eight;
      }
  }

};

template<typename T>
struct bms_template_shape_legendre_splz<1,WMESH_ELEMENT_PYRAMID,T>
{


  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {

    static constexpr T four = static_cast<T>(4);
    static constexpr T one = static_cast<T>(1);
    static constexpr T zero = static_cast<T>(0);
    const T
      r = lcoo_[lcoo_inc_*0],
      s = lcoo_[lcoo_inc_*1],
      t = lcoo_[lcoo_inc_*2];
    
    if (diff_[0] > 0)
      {
	if (t < 1.0)
	  {
	    b_[b_inc_*0] = ( -one + s / (one-t)) / four;
	    b_[b_inc_*1] = (  one  - s / (one-t)) / four;
	    b_[b_inc_*2] = (  one  + s / (one-t)) / four;
	    b_[b_inc_*3] = ( -one - s / (one-t)) / four;
	    b_[b_inc_*4] = zero;
	  }
	else
	  {
	    b_[b_inc_*0] = -one / four;
	    b_[b_inc_*1] = one / four;
	    b_[b_inc_*2] = one / four;
	    b_[b_inc_*3] = -one / four;
	    b_[b_inc_*4] = zero;
	  }
      }
    else if (diff_[1] > 0)
      {
	if (t < 1.0)
	  {
	    b_[b_inc_*0] = ( -one + r / (one-t)) / four;
	    b_[b_inc_*1] = ( -one  - r / (one-t)) / four;
	    b_[b_inc_*2] = (  one  + r / (one-t)) / four;
	    b_[b_inc_*3] = (  one  - r / (one-t)) / four;
	    b_[b_inc_*4] = zero;
	  }
	else
	  {
	    b_[b_inc_*0] = -one / four;
	    b_[b_inc_*1] = -one / four;
	    b_[b_inc_*2] = one / four;
	    b_[b_inc_*3] = one / four;
	    b_[b_inc_*4] = zero;
	  }
      }
    else if (diff_[2] > 0)
      {
	if (t < 1.0)
	  {
	    b_[b_inc_*0] = ( -one + r*s / (one-t) / (one-t)) / four;
	    b_[b_inc_*1] = ( -one - r*s / (one-t) / (one-t)) / four;
	    b_[b_inc_*2] = ( -one + r*s / (one-t) / (one-t)) / four;
	    b_[b_inc_*3] = ( -one - r*s / (one-t) / (one-t)) / four;
	    b_[b_inc_*4] = one;
	  }
	else
	  {
	    b_[b_inc_*0] = -one / four;
	    b_[b_inc_*1] = -one / four;
	    b_[b_inc_*2] = -one / four;
	    b_[b_inc_*3] = -one / four;
	    b_[b_inc_*4] = one;
	  }
      }
    else
      {
	if (t < 1.0)
	  {
	    b_[b_inc_*0] = ( one - r - s - t + r*s / (one-t)) / four;
	    b_[b_inc_*1] = ( one + r - s - t - r*s / (one-t)) / four;
	    b_[b_inc_*2] = ( one + r + s - t + r*s / (one-t)) / four;
	    b_[b_inc_*3] = ( one - r + s - t - r*s / (one-t)) / four;
	    b_[b_inc_*4] = t;
	  }
	else
	  {
	    b_[b_inc_*0] = zero;
	    b_[b_inc_*1] = zero;
	    b_[b_inc_*2] = zero;
	    b_[b_inc_*3] = zero;
	    b_[b_inc_*4] = one;
	  }
      }
  }

  
};
