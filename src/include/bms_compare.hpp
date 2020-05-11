#pragma once
#include "wmesh-types.h"

template<wmesh_int_t XTYPE>
struct bms_compare
{
public:
  
  static inline bool 		comp_dna(const_wmesh_int_p this_dna_,
					 const_wmesh_int_p that_dna_) noexcept;
  
  static inline void 		hash_dna(const_wmesh_int_p	x2n_,
					 wmesh_int_p 	dna_) noexcept;
  
  static inline wmesh_int_t 	hash(const_wmesh_int_p 	dna_,
				     wmesh_int_t	bound_) noexcept;
};

template<>
struct bms_compare<WMESH_ELEMENT_QUADRILATERAL>
{
public:
  
  static inline bool comp_dna(const_wmesh_int_p this_dna_,
			      const_wmesh_int_p that_dna_) noexcept
  {
    return ((this_dna_[0] == that_dna_[0]) &&(this_dna_[1] == that_dna_[1]));
  }
  
  static inline void hash_dna(const_wmesh_int_p	this_x2n_,
			      wmesh_int_p 	dna_) noexcept
  {
    wmesh_int_t a = this_x2n_[0];
    int ia = 0;
    for (int i=1;i<4;++i)
      {
	if (a > this_x2n_[i])
	  {
	    a = this_x2n_[i];
	    ia = i;
	  }
      }
    dna_[0] = a;
    dna_[1] = this_x2n_[(ia+2)%4];
  };  
  
  static inline wmesh_int_t 	hash	(const_wmesh_int_p 	dna_,
					 wmesh_int_t 		bound_) noexcept
  {
    const wmesh_int_t h = (31 * dna_[0] + 57*dna_[1]) % bound_;
    return (h<0)
      ? -h
      : h;
  };  
};

template<>
struct bms_compare<WMESH_ELEMENT_TRIANGLE>
{
public:
  static inline bool comp_dna(const_wmesh_int_p this_dna_,
			      const_wmesh_int_p that_dna_) noexcept
  {
    return ((this_dna_[0] == that_dna_[0]) && (this_dna_[1] == that_dna_[1]) && (this_dna_[2] == that_dna_[2]));
  }
  
  static inline void hash_dna(const_wmesh_int_p	c2n_,
			      wmesh_int_p 	dna_) noexcept
  {
    wmesh_int_t a 	= (c2n_[0] < c2n_[1]) ? c2n_[0] : c2n_[1];
    wmesh_int_t b 	= (c2n_[0] < c2n_[1]) ? c2n_[1] : c2n_[0];    
    a 			= (a < c2n_[2]) ? a : c2n_[2];
    b 			= (b < c2n_[2]) ? c2n_[2] : b;
    wmesh_int_t sum 	= c2n_[0]+c2n_[1]+c2n_[2];    

    dna_[0] = a;
    dna_[1] = b;
    dna_[2] = sum;
    
  };  
  
  static inline wmesh_int_t 	hash	(const_wmesh_int_p dna_,
					 wmesh_int_t bound_) noexcept
  {
    const wmesh_int_t h = (31 * dna_[0] + 57*dna_[1] + 79*dna_[2]) % bound_;
    return (h<0)
      ? -h
      : h;
  };  
};


template<>
struct bms_compare<WMESH_ELEMENT_EDGE>
{
public:
  
  static inline bool comp_dna(const_wmesh_int_p this_dna_,
			      const_wmesh_int_p that_dna_) noexcept
  {
    return ( (this_dna_[0] == that_dna_[0]) &&(this_dna_[1] == that_dna_[1]) );
  }
  
  static inline void 		hash_dna(const_wmesh_int_p	x2n_,
					 wmesh_int_p 	dna_) noexcept
  {
    bool ordered = (x2n_[0] < x2n_[1]);
    dna_[0] 	= ordered ? x2n_[0] : x2n_[1];
    dna_[1] 	= ordered ? x2n_[1] : x2n_[0];    
  };  
  
  static inline wmesh_int_t 	hash(const_wmesh_int_p 	dna_,
				     wmesh_int_t	bound_) noexcept
  {    
    const wmesh_int_t h = (31 * dna_[0] + 57*dna_[1]) % bound_;
    return (h<0)
      ? -h
      : h;
  };  
};

