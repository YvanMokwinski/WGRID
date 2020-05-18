#pragma once
#include <ostream>
#include "storage_t.h"

namespace WCOMMON
{
  template <typename integer_t>
  struct storage_t
  {
  public: typedef enum ekind : integer_t    
    {
      unknown = WCOMMON_STORAGE_UNKNOWN,
	block = WCOMMON_STORAGE_BLOCK,
	interleave = WCOMMON_STORAGE_INTERLEAVE
	} kind_t;
    
  public: static constexpr kind_t all[3] = {unknown,
					    block,
					    interleave};
  public: inline static bool is_invalid(kind_t kind_) 
    {
      return kind_ != block && kind_ != interleave;
    };
    
  public: inline static bool is_invalid(integer_t kind_) 
    {
      return kind_ != block && kind_ != interleave;
    };
  public: inline storage_t(kind_t kind_) noexcept;
  public: inline storage_t() noexcept;
  public: inline operator kind_t() const noexcept;
  private: kind_t m_kind{};
  };
  
  template <typename integer_t>
  inline storage_t<integer_t>::storage_t(storage_t::kind_t kind_) noexcept : m_kind(kind_) {};
  template <typename integer_t>
  inline storage_t<integer_t>::storage_t() noexcept {};
  template <typename integer_t>
  inline storage_t<integer_t>::operator storage_t::kind_t() const noexcept { return this->m_kind; };
  
#if 0
  std::ostream& operator<<(std::ostream&out_,
			   const storage_t::kind_t&kind_)
  {
#ifdef define_case
#error macro 'define_case' is already defined.
#endif
#define define_case(_kind) case storage_t::_kind: out_ << #_kind; break
    switch(kind_)
      {
	define_case(unknown);
	define_case(block);
	define_case(interleave);
      }
#undef define_case
    return out_;
  };
#endif
  
};

namespace WFE
{
  using storage_t = WCOMMON::storage_t<wfe_int_t>;
}
