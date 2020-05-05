#pragma once

#include <limits>

template <typename _int_t,int _nbits> struct GenericEncoding
{

private: static constexpr const bool s_is_signed = std::is_signed<_int_t>::value;
private: static constexpr const _int_t s_one = _int_t(1);
private: static constexpr const _int_t s_two = 2;  
  //!
  //!
  //! 
public: static constexpr int s_totalBits = sizeof(_int_t)*8 - (s_is_signed ? 1 : 0);
public: static constexpr int s_availableBits = s_totalBits - _nbits;  
  
public: static constexpr const _int_t s_lowLimit = (s_one << _nbits);  
public: static constexpr const _int_t s_upLimit = std::numeric_limits<_int_t>::max()/s_lowLimit+1;
public: static constexpr const _int_t s_signedValue = s_is_signed ? (s_one << s_totalBits) : 0;
  
private: static constexpr _int_t s_c20 = ( ((s_one << (_nbits+(s_is_signed ? 1 : 0))) - s_one) << s_availableBits );
private: static constexpr _int_t s_c2 = (_int_t)~s_c20;
  
public: static inline void Info() noexcept
  {
    std::cout << "nbits     " << s_totalBits <<std::endl;
    std::cout << "Low nbits " << _nbits <<std::endl;
    std::cout << "Up  nbits " << s_availableBits <<std::endl;
    std::cout << "Low limit " << s_lowLimit <<std::endl;
    std::cout << "Up  limit " << s_upLimit <<std::endl;
    
    _int_t s = s_signedValue;
    std::cout << "SignedValue " << s << std::endl;
    PrintBits(sizeof(_int_t),&s);
    std::cout << "~s_c2" << std::endl;
    s = ~s_c2;
    PrintBits(sizeof(_int_t),&s);
    std::cout << "s_c2" << std::endl;
    s = s_c2;
    PrintBits(sizeof(_int_t),&s);
    // std::cout << "s_c2 " << s_c <<std::endl;
  };

private: template <typename R,typename T = _int_t> using enable_if_signed_t = typename std::enable_if<std::is_signed<T>::value,R >::type;
private: template <typename R,typename T = _int_t> using enable_if_unsigned_t = typename std::enable_if<!std::is_signed<T>::value,R >::type;

private: template <typename T = _int_t> using signed__int_t = enable_if_signed_t<T>;
private: template <typename T = _int_t> using unsigned__int_t = enable_if_unsigned_t<T>;
  
  //!
  //!
  //!
#if 0
public: template <typename T = _int_t>
static inline constexpr typename std::enable_if<std::is_signed<T>::value,bool >::type
IsPositive(const T e_) noexcept
  {
    return (e_ & s_signedValue) == 0;
  };
#endif  
  
public: template <typename T = _int_t>
static inline constexpr signed__int_t<T> Low(const T e_) noexcept
  {
    //    return  ( ( IsPositive(e_) ? e_ : s_signedValue ^ e_ ) >> s_availableBits );
    return e_  >> s_availableBits;
  };
  
public: template <typename T = _int_t>
static inline constexpr unsigned__int_t<T> Low(const T e_) noexcept
  {
    return e_ >> s_availableBits;
  };

  
public: template <typename T = _int_t>
static inline constexpr signed__int_t<T> Up(const T e_) noexcept
  {
    //    return (IsPositive(e_) ?  e_ & s_c2 : -(e_ & s_c2));
        return e_ & s_c2;
  };
  
public: template <typename T = _int_t>
static inline constexpr unsigned__int_t<T> Up(const T e_) noexcept
  {
    return e_ & s_c2;
  };
  
  
public: template <typename T = _int_t>
static inline constexpr signed__int_t<T> Encod(const T up_,
					      const T low_) noexcept
  {
#if 0
    return
      IsPositive(up_)
      ? (up_  | (low_ << s_availableBits) )
      : (-up_ | (low_ << s_availableBits) ) | s_signedValue;
#endif
    return (up_  | (low_ << s_availableBits) );
  };
  
public: template <typename T = _int_t>
static inline constexpr unsigned__int_t<T> Encod(const T up_,
						const T low_) noexcept
  {
    return up_  | (low_ << s_availableBits);
  };
  
  
  // assumes little endian
public: static void PrintBits(size_t const size, void const * const ptr)
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

  
public: template <typename T = _int_t>
static inline typename std::enable_if<std::is_signed<T>::value>::type Test()
  {
    _int_t up, low;
    std::cout << "testing from -100 to 100"  << std::endl;
    for (up=-100;up<=100;++up)
      {
	for (low = 0;low < 4;++low)
	  {
	    auto encodedValue = Encod(up, low);
	    if (Low(encodedValue) != low)
	      {
		std::cout << "low is wrong, expected " << low << " but got " <<  Low(encodedValue) << std::endl;
		exit(1);
	      }
	    if (Up(encodedValue) != up)
	      {
		std::cout << "up is wrong, expected " << up << " but got " <<  Up(encodedValue) << std::endl;
		exit(1);
	      }
	    
	  }
      }
    std::cout << "passed"
	      << std::endl
	      << std::endl;

    std::cout << "testing ["<< s_upLimit - 3 << "," << s_upLimit << "[" << std::endl;    
    std::cout << "testing from "<< s_upLimit - 3 << " to " << s_upLimit << std::endl;    
    for (up=s_upLimit-3;up < s_upLimit;++up)
      {
	for (low = 0;low < 4;++low)
	  {
	    auto encodedValue = Encod(up, low);
	    if (Low(encodedValue) != low)
	      {
		std::cout << "low is wrong, expected " << low << " but got " <<  Low(encodedValue) << std::endl;
		exit(1);
	      }
	    if (Up(encodedValue) != up)
	      {
		std::cout << "up is wrong, expected " << up << " but got " <<  Up(encodedValue) << std::endl;
		exit(1);
	      }
	    
	  }
      }
    std::cout << "passed"
	      << std::endl
	      << std::endl;

    std::cout << "testing ["<< -(s_upLimit - 1) << "," << -(s_upLimit-10) << "[" << std::endl;    
    for (up=-(s_upLimit-1);up < -(s_upLimit-10);++up)
      {	
	for (low = 0;low < 4;++low)
	  {
	    auto encodedValue = Encod(up,low);
	    if (Low(encodedValue) != low)
	      {
		std::cout << "low is wrong, expected " << low << " but got " <<  Low(encodedValue) << std::endl;
		exit(1);
	      }
	    if (Up(encodedValue) != up)
	      {
		std::cout << "up is wrong, expected " << up << " but got " <<  Up(encodedValue) << std::endl;
		exit(1);
	      }	    
	  }
      }

    std::cout << "passed"
	      << std::endl
	      << std::endl;

  };

  
public: template <typename T = _int_t>
static inline typename std::enable_if<!std::is_signed<T>::value >::type Test()
  {    
    T up, low;
    std::cout << "testing from 0 to 100"  << std::endl;
    for (up=0;up<=100;++up)
      {
	for (low = 0;low < 4;++low)
	  {
	    auto encodedValue = Encod(up, low);
	    if (Low(encodedValue) != low)
	      {
		std::cout << "low is wrong, expected " << low << " but got " <<  Low(encodedValue) << std::endl;
		exit(1);
	      }
	    if (Up(encodedValue) != up)
	      {
		std::cout << "up is wrong, expected " << up << " but got " <<  Up(encodedValue) << std::endl;
		exit(1);
	      }
	    
	  }
      }
    std::cout << "passed"
	      << std::endl
	      << std::endl;

    std::cout << "testing ["<< s_upLimit - 3 << "," << s_upLimit << "[" << std::endl;    
    std::cout << "testing from "<< s_upLimit - 3 << " to " << s_upLimit << std::endl;    
    for (up=s_upLimit-3;up < s_upLimit;++up)
      {
	for (low = 0;low < 4;++low)
	  {
	    auto encodedValue = Encod(up, low);
	    if (Low(encodedValue) != low)
	      {
		std::cout << "low is wrong, expected " << low << " but got " <<  Low(encodedValue) << std::endl;
		exit(1);
	      }
	    if (Up(encodedValue) != up)
	      {
		std::cout << "up is wrong, expected " << up << " but got " <<  Up(encodedValue) << std::endl;
		exit(1);
	      }
	    
	  }
      }
    std::cout << "passed"
	      << std::endl
	      << std::endl;


  };
  
};
