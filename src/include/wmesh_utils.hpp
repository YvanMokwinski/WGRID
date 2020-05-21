#pragma once
#include "wmesh.h"

template <wmesh_int_t ELEMENT> inline wmesh_int_t wfe_ndofs_template(wmesh_int_t d_) {return -1;}

template <> inline wmesh_int_t wfe_ndofs_template<WMESH_ELEMENT_EDGE>(wmesh_int_t d_) { return d_+1; }

template <> inline wmesh_int_t wfe_ndofs_template<WMESH_ELEMENT_TRIANGLE>(wmesh_int_t d_) { return ((d_+1)*(d_+2)) / 2; }
template <> inline wmesh_int_t wfe_ndofs_template<WMESH_ELEMENT_QUADRILATERAL>(wmesh_int_t d_) { return ((d_+1)*(d_+1)); }

template <> inline wmesh_int_t wfe_ndofs_template<WMESH_ELEMENT_TETRAHEDRON>(wmesh_int_t d_) { return ((d_+1)*(d_+2)*(d_+3)) / 6; }
template <> inline wmesh_int_t wfe_ndofs_template<WMESH_ELEMENT_PYRAMID>(wmesh_int_t d_) { return ((d_+2)*(d_+1)*(2*d_+3)) / 6; }
template <> inline wmesh_int_t wfe_ndofs_template<WMESH_ELEMENT_WEDGE>(wmesh_int_t d_) { return ((d_+1)*(d_+1)*(d_+2)) / 2; }
template <> inline wmesh_int_t wfe_ndofs_template<WMESH_ELEMENT_HEXAHEDRON>(wmesh_int_t d_) { return (d_+1)*(d_+1)*(d_+1); }




template<typename T>
wmesh_status_t wmesh_raw_bounding_box(wmesh_int_t 		coo_m_,
				      wmesh_int_t 		coo_n_,
				      const T *__restrict__ 	coo_v_,
				      wmesh_int_t 		coo_ld_,
				      T *__restrict__ 		box_);


extern "C"
{

  wmesh_status_t wfe_ndofs(wmesh_int_t element_,
			   wmesh_int_t d_,
			   wmesh_int_p ndofs_);


  const char * file_extension(const char * filename_);
  
  

  
  
unsigned long long int hilbert_coordinate(double*__restrict__ 	crd,
					  const double*__restrict__ box,
					  int 	itr);


}


#if 0
#include <limits>

template <typename _int_t,int _nbits> struct encoding
{
private: static constexpr const bool s_is_signed = std::is_signed<_int_t>::value;
private: static constexpr const _int_t s_one = _int_t(1);
private: static constexpr const _int_t s_two = 2;  
public: static constexpr int s_totalBits = sizeof(_int_t)*8 - (s_is_signed ? 1 : 0);
public: static constexpr int s_availableBits = s_totalBits - _nbits;  
public: static constexpr const _int_t s_lowLimit = (s_one << _nbits);  
public: static constexpr const _int_t s_upLimit = std::numeric_limits<_int_t>::max()/s_lowLimit+1;
public: static constexpr const _int_t s_signedValue = s_is_signed ? (s_one << s_totalBits) : 0;  
private: static constexpr _int_t s_c20 = ( ((s_one << (_nbits+(s_is_signed ? 1 : 0))) - s_one) << s_availableBits );
private: static constexpr _int_t s_c2 = (_int_t)~s_c20;
public: static inline void info() noexcept
  {
    std::cout << "nbits     " << s_totalBits <<std::endl;
    std::cout << "Low nbits " << _nbits <<std::endl;
    std::cout << "Up  nbits " << s_availableBits <<std::endl;
    std::cout << "Low limit " << s_lowLimit <<std::endl;
    std::cout << "Up  limit " << s_upLimit <<std::endl;
    
    _int_t s = s_signedValue;
    std::cout << "SignedValue " << s << std::endl;
    print(sizeof(_int_t),&s);
    std::cout << "~s_c2" << std::endl;
    s = ~s_c2;
    print(sizeof(_int_t),&s);
    std::cout << "s_c2" << std::endl;
    s = s_c2;
    print(sizeof(_int_t),&s);
    // std::cout << "s_c2 " << s_c <<std::endl;
  };

private: template <typename R,typename T = _int_t> using enable_if_signed_t = typename std::enable_if<std::is_signed<T>::value,R >::type;
private: template <typename R,typename T = _int_t> using enable_if_unsigned_t = typename std::enable_if<!std::is_signed<T>::value,R >::type;

private: template <typename T = _int_t> using signed__int_t = enable_if_signed_t<T>;
private: template <typename T = _int_t> using unsigned__int_t = enable_if_unsigned_t<T>;
  
  //!
  //!
  //!
public: template <typename T = _int_t>
static inline constexpr typename std::enable_if<std::is_signed<T>::value,bool >::type
is_positive(const T e_) noexcept
  {
    return (e_ & s_signedValue) == 0;
  };
  
  
public: template <typename T = _int_t>
static inline constexpr signed__int_t<T> low_value(const T e_) noexcept
  {
    return  ( ( is_positive(e_) ? e_ : s_signedValue ^ e_ ) >> s_availableBits );
  };
  
public: template <typename T = _int_t>
static inline constexpr unsigned__int_t<T> low_value(const T e_) noexcept
  {
    return e_ >> s_availableBits;
  };

  
public: template <typename T = _int_t>
static inline constexpr signed__int_t<T> up_value(const T e_) noexcept
  {
    return (is_positive(e_) ?  e_ & s_c2 : -(e_ & s_c2));
  };
  
public: template <typename T = _int_t>
static inline constexpr unsigned__int_t<T> up_value(const T e_) noexcept
  {
    return e_ & s_c2;
  };
  
  
public: template <typename T = _int_t>
static inline constexpr signed__int_t<T> encod(const T up_,
					       const T low_) noexcept
  {
    return
      is_positive(up_)
      ? (up_  | (low_ << s_availableBits) )
      : (-up_ | (low_ << s_availableBits) ) | s_signedValue;      
  };
  
public: template <typename T = _int_t>
static inline constexpr unsigned__int_t<T> encod(const T up_,
						 const T low_) noexcept
  {
    return up_  | (low_ << s_availableBits);
  };
  
  // assumes little endian
public: static void print(size_t const size, void const * const ptr)
  {
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;
    printf("##\n");
    for (i=size-1;i>=0;i--)
      {
	for (j=7;j>=0;j--)
	  {
	    byte = b[i] & (1<<j);
	    byte >>= j;
	    printf("%u", byte);
	  }
      }
    puts("");
  };

  
};
#endif
