#pragma once

template<wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_template_shape_monomial_functor
{
  wmesh_int_t 	m_degree;
  void basis(const_wmesh_int_p 	diff_,
	     const T * 		lcoo_,
	     wmesh_int_t 	lcoo_inc_,
	     T * 		b_,
	     wmesh_int_t 	b_inc_)
  {    
  };
  
  bms_template_shape_monomial_functor(wmesh_int_t degree_) : m_degree(degree_)
  {    
  };  

};





template<typename 	T>
struct bms_template_shape_monomial_functor<WMESH_ELEMENT_EDGE,T>
{
  wmesh_int_t 	m_degree;
  void basis(const_wmesh_int_p 	diff_,
	     const T * 		lcoo_,
	     wmesh_int_t 	lcoo_inc_,
	     T * 		b_,
	     wmesh_int_t 	b_inc_)
  {
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = static_cast<T>(1);
	if (m_degree>0)
	  {
	    for (wmesh_int_t i=1;i<=m_degree;++i)
	      {
		b_[b_inc_*i] = b_[b_inc_*(i-1)] * lcoo_[lcoo_inc_*0];
	      }
	  }
      }
    else
      {
	b_[b_inc_*0] = static_cast<T>(0);
	if (m_degree>0)
	  {
	    b_[b_inc_*1] = static_cast<T>(1);
	    for (wmesh_int_t i=2;i<=m_degree;++i)
	      {
		b_[b_inc_*i] = b_[b_inc_*(i-1)] * lcoo_[lcoo_inc_*0] * static_cast<T>(i);
	      }
	  }
      }
  };
  
  bms_template_shape_monomial_functor(wmesh_int_t degree_) : m_degree(degree_)
  {    
  };  

};


template<typename 	T>
struct bms_template_shape_monomial_functor<WMESH_ELEMENT_TRIANGLE,T>
{
  wmesh_int_t 	m_degree;
  void basis(const_wmesh_int_p 	diff_,
	     const T * 		lcoo_,
	     wmesh_int_t 	lcoo_inc_,
	     T * 		b_,
	     wmesh_int_t 	b_inc_)
  {
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = static_cast<T>(1);
	if (m_degree>0)
	  {
	    for (wmesh_int_t i=1;i<=m_degree;++i)
	      {
		b_[b_inc_*i] = b_[b_inc_*(i-1)] * lcoo_[lcoo_inc_*0];
	      }
	  }
      }
    else if (diff_[0] > 0)
      {
	b_[b_inc_*0] = static_cast<T>(1);
	if (m_degree>0)
	  {
	    for (wmesh_int_t i=1;i<=m_degree;++i)
	      {
		b_[b_inc_*i] = b_[b_inc_*(i-1)] * lcoo_[lcoo_inc_*0];
	      }
	  }
      }
    else
      {
	static constexpr T one  = static_cast<T>(1);
	b_[b_inc_*0] = one;
	if (this->m_degree > 0)
	  {
	    wmesh_int_t shift = 1;

	    const T
	      r = lcoo_[lcoo_inc_* 0],
	      s = lcoo_[lcoo_inc_* 1];
		
	    b_[b_inc_*1] = r;
	    b_[b_inc_*2] = s;

	    if (this->m_degree>1)
	      {
		for (wmesh_int_t i=2;i<=this->m_degree;++i)
		  {		    
		    {
		      wmesh_int_t j = 0;
		      b_[b_inc_*i] = b_[b_inc_*(shift + j)] * r;
		    }
		    
		    for (wmesh_int_t j=1;j<i;++i)
		      {
			b_[b_inc_*(shift + (i-1) )] = b_[b_inc_*(shift + (i-1))] * b_[b_inc_*(shift + i)];
			// b_[b_inc_*(shift + (i-1) )] = b_[b_inc_*(shift + (i-1))] * b_[b_inc_*(shift + i)];
		      }
			   
		    {
		      wmesh_int_t j = i;
		      b_[b_inc_*i] = b_[b_inc_*(shift + j)] * s;
		    }
		    shift+=i;
		  }		
	      }
	  }
	
      }
    
  };
  
  bms_template_shape_monomial_functor(wmesh_int_t degree_) : m_degree(degree_)
  {    
  };  

};





template<typename 	T>
struct bms_template_shape_monomial_functor<WMESH_ELEMENT_QUADRILATERAL,T>
{
  wmesh_int_t 	m_degree;
  void basis(const_wmesh_int_p 	diff_,
	     const T * 		lcoo_,
	     wmesh_int_t 	lcoo_inc_,
	     T * 		b_,
	     wmesh_int_t 	b_inc_)
  {
    
  };
  
  bms_template_shape_monomial_functor(wmesh_int_t degree_) : m_degree(degree_)
  {    
  };  

};


template<typename 	T>
struct bms_template_shape_monomial_functor<WMESH_ELEMENT_TETRAHEDRON,T>
{
  wmesh_int_t 	m_degree;
  void basis(const_wmesh_int_p 	diff_,
	     const T * 		lcoo_,
	     wmesh_int_t 	lcoo_inc_,
	     T * 		b_,
	     wmesh_int_t 	b_inc_)
  {
    
  };
  
  bms_template_shape_monomial_functor(wmesh_int_t degree_) : m_degree(degree_)
  {    
  };  

};


template<typename 	T>
struct bms_template_shape_monomial_functor<WMESH_ELEMENT_PYRAMID,T>
{
  wmesh_int_t 	m_degree;
  void basis(const_wmesh_int_p 	diff_,
	     const T * 		lcoo_,
	     wmesh_int_t 	lcoo_inc_,
	     T * 		b_,
	     wmesh_int_t 	b_inc_)
  {
    
  };
  
  bms_template_shape_monomial_functor(wmesh_int_t degree_) : m_degree(degree_)
  {    
  };  

};


template<typename 	T>
struct bms_template_shape_monomial_functor<WMESH_ELEMENT_WEDGE,T>
{
  wmesh_int_t 	m_degree;
  void basis(const_wmesh_int_p 	diff_,
	     const T * 		lcoo_,
	     wmesh_int_t 	lcoo_inc_,
	     T * 		b_,
	     wmesh_int_t 	b_inc_)
  {
    
  };
  
  bms_template_shape_monomial_functor(wmesh_int_t degree_) : m_degree(degree_)
  {    
  };  

};


template<typename 	T>
struct bms_template_shape_monomial_functor<WMESH_ELEMENT_HEXAHEDRON,T>
{
  wmesh_int_t 	m_degree;
  void basis(const_wmesh_int_p 	diff_,
	     const T * 		lcoo_,
	     wmesh_int_t 	lcoo_inc_,
	     T * 		b_,
	     wmesh_int_t 	b_inc_)
  {
    
  };
  
  bms_template_shape_monomial_functor(wmesh_int_t degree_) : m_degree(degree_)
  {    
  };  

};



template<wmesh_int_t 	degree_,
	 wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_template_shape_monomial_splz
{  
  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_);
};


template<wmesh_int_t ELEMENT_,typename T>
struct bms_template_shape_monomial_splz<0,ELEMENT_,T>
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
struct bms_template_shape_monomial_splz<1,WMESH_ELEMENT_EDGE,T>
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
      r = lcoo_[lcoo_inc_* 0];
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = one;
	b_[b_inc_*1] = r;
      }
    else
      {
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = one;
      }    
  }
};


template<typename T>
struct bms_template_shape_monomial_splz<1,WMESH_ELEMENT_TRIANGLE,T>
{

  static inline void basis(const_wmesh_int_p 	diff_,
			   const T * 		lcoo_,
			   wmesh_int_t 		lcoo_inc_,
			   T * 			b_,
			   wmesh_int_t 		b_inc_)
  {
    static constexpr T one = static_cast<T>(1);
    static constexpr T zero = static_cast<T>(0);
    const T r = lcoo_[lcoo_inc_ * 0];
    const T s = lcoo_[lcoo_inc_ * 1];    
    if (diff_[0] == 1)
      {
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = one;
	b_[b_inc_*2] = zero;	
      }
    else if (diff_[1] == 1)
      {
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = zero;
	b_[b_inc_*2] = one;	
      }
    else
      {
	b_[b_inc_*0] = one;
	b_[b_inc_*1] = r;
	b_[b_inc_*2] = s;
      }
  }
};

template<typename T>
struct bms_template_shape_monomial_splz<1,WMESH_ELEMENT_QUADRILATERAL,T>
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
      r = lcoo_[lcoo_inc_* 0],
      s = lcoo_[lcoo_inc_* 1];
    
    if (diff_[0] > 0)
      {
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = one;
	b_[b_inc_*2] = zero;
	b_[b_inc_*3] = s;
      }
    else if (diff_[1] > 0)
      {
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = zero;
	b_[b_inc_*2] = one;
	b_[b_inc_*3] = r;
      }
    else
      {
	b_[b_inc_*0] = one;
	b_[b_inc_*1] = r;
	b_[b_inc_*2] = s;
	b_[b_inc_*3] = r*s;
      }
  }


};

template<typename T>
struct bms_template_shape_monomial_splz<1,WMESH_ELEMENT_TETRAHEDRON,T>
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
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = one;
	b_[b_inc_*2] = zero;
	b_[b_inc_*3] = zero;
      }
    else if (diff_[1] > 0)
      {
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = zero;
	b_[b_inc_*2] = one;
	b_[b_inc_*3] = zero;
      }
    else if (diff_[2] > 0)
      {
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = zero;
	b_[b_inc_*2] = zero;
	b_[b_inc_*3] = one;
      }
    else
      {
	b_[b_inc_*0] = one;
	b_[b_inc_*1] = r;
	b_[b_inc_*2] = s;
	b_[b_inc_*3] = t;
      }
  }



};

template<typename T>
struct bms_template_shape_monomial_splz<1,WMESH_ELEMENT_WEDGE,T>
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
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = one;
	b_[b_inc_*2] = zero;
	b_[b_inc_*3] = zero;
	b_[b_inc_*4] = t;
	b_[b_inc_*5] = zero;
      }
    else if (diff_[1] > 0)
      {
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = zero;
	b_[b_inc_*2] = one;
	b_[b_inc_*3] = zero;
	b_[b_inc_*4] = zero;
	b_[b_inc_*5] = t;
      }
    else if (diff_[2] > 0)
      {
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = zero;
	b_[b_inc_*2] = zero;
	b_[b_inc_*3] = one;
	b_[b_inc_*4] = r;
	b_[b_inc_*5] = s;
      }
    else
      {
	b_[b_inc_*0] = one;
	b_[b_inc_*1] = r;
	b_[b_inc_*2] = s;
	b_[b_inc_*3] = t;
	b_[b_inc_*4] = r*t;
	b_[b_inc_*5] = s*t;
      }
  }


};

template<typename T>
struct bms_template_shape_monomial_splz<1,WMESH_ELEMENT_HEXAHEDRON,T>
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
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = one;
	b_[b_inc_*2] = zero;
	b_[b_inc_*3] = zero;
	b_[b_inc_*4] = s;
	b_[b_inc_*5] = t;
	b_[b_inc_*6] = zero;
	b_[b_inc_*7] = s*t;
      }
    else if (diff_[1] > 0)
      {
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = zero;
	b_[b_inc_*2] = one;
	b_[b_inc_*3] = zero;
	b_[b_inc_*4] = r;
	b_[b_inc_*5] = zero;
	b_[b_inc_*6] = t;
	b_[b_inc_*7] = r*t;
      }
    else if (diff_[2] == 1)
      {
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = zero;
	b_[b_inc_*2] = zero;
	b_[b_inc_*3] = one;
	b_[b_inc_*4] = zero;
	b_[b_inc_*5] = r;
	b_[b_inc_*6] = s;
	b_[b_inc_*7] = r*s;
      }
    else
      {
	b_[b_inc_*0] = one;
	b_[b_inc_*1] = r;
	b_[b_inc_*2] = s;
	b_[b_inc_*3] = t;
	b_[b_inc_*4] = r*s;
	b_[b_inc_*5] = r*t;
	b_[b_inc_*6] = s*t;
	b_[b_inc_*7] = r*s*t;
      }
  }

};

template<typename T>
struct bms_template_shape_monomial_splz<1,WMESH_ELEMENT_PYRAMID,T>
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
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = one;
	b_[b_inc_*2] = zero;
	b_[b_inc_*3] = zero;
	b_[b_inc_*4] = (t<one) ? s / (one-t) : one;
      }
    else if (diff_[1] > 0)
      {
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = zero;
	b_[b_inc_*2] = one;
	b_[b_inc_*3] = zero;
	b_[b_inc_*4] = (t<one) ? r / (one-t) : one;
      }
    else if (diff_[2] > 0)
      {
	b_[b_inc_*0] = zero;
	b_[b_inc_*1] = zero;
	b_[b_inc_*2] = zero;
	b_[b_inc_*3] = one;
	b_[b_inc_*4] = (t<one) ? (r*s / (one-t) / (one-t)) : one;
      }
    else
      {
	b_[b_inc_*0] = one;
	b_[b_inc_*1] = r;
	b_[b_inc_*2] = s;
	b_[b_inc_*3] = t;
	b_[b_inc_*4] = (t<one) ? r*s / (one-t) : one;
      }
  }

  
};
