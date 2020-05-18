#pragma once
#include <math.h>

//! @brief Set of mathematical static functions. 
//! @tparam _float_type The real type to use.
template <typename T> struct wmesh_math
{
  using real_t = T;
  //! @brief Compute the squared root.
  //! @param x_ The value from which we compute the squared root.
  //! @return The squared root.
  static constexpr real_t sqrt(const real_t&x_) noexcept;
  //! @brief Compute the absolute value.
  //! @param x_ The value from which we compute the absolue value.
  //! @return The absolute value.
  static constexpr real_t abs(const real_t&x_) noexcept;
  //! @brief Compute the power x^y.
  //! @param x_ The base of the power.
  //! @param y_ The exponent of the power.
  //! @return The power x^y.
  static constexpr real_t pow(const real_t&x_,const real_t&y_) noexcept;
  //! @brief Compute the maximum between two real.
  //! @param a_ The first operand.
  //! @param b_ The second operand.
  //! @return MAX(a,b)
  static constexpr real_t max(const real_t&a_,const real_t&b_) noexcept;
  //! @brief Compute the maximum between two real.
  //! @param a_ The first operand.
  //! @param b_ The second operand.
  //! @return MAX(a,b)
  static constexpr real_t min(const real_t&a_,const real_t&b_) noexcept;

  //! @brief Set the value to zero if the absolute value is less than a tolerance.
  //! @param The value.
  //! @param The tolerance.
  //! @return The filtered value.
  static constexpr real_t filter(const real_t&value_,const real_t&tolerance_) noexcept;

};

//! @copydoc wmesh_math
template <> struct wmesh_math<double>
{
  static constexpr unsigned int num_digits = 15;
  static constexpr double zero = 0.0;
  static constexpr double epsilon = 2.22044604925031308e-16;
  //! @copydoc wmesh_math::pow(const real_t&x_,const real_t&y_)
  static constexpr double pow(const double&x_,const double&y_) noexcept
  {
    return pow(x_,y_);
  };
  //! @copydoc wmesh_math::sqrt(const real_t&x_)
  static double sqrt(const double&x_) noexcept
  {
    return  sqrt(x_);
  };
  //! @copydoc wmesh_math::abs(const real_t&x_)
  static constexpr double abs(const double&x_) noexcept
  {
    return  (x_<0.0) ? -x_ : x_;
  };
  //! @copydoc wmesh_math::max(const real_t&,const real_t&)
  static constexpr double max(const double&a_,const double&b_) noexcept
  {
    return (a_>=b_) ? a_ : b_;
  };
  //! @copydoc wmesh_math::min(const real_t&,const real_t&)
  static constexpr double min(const double&a_,const double&b_) noexcept
  {
    return (a_<=b_) ? a_ : b_;
  };

  //! @copydoc wmesh_math::filter(const real_t&,const real_t&)
  static constexpr double filter(const double&value_,const double&tolerance_) noexcept
  {
    return (value_ < -tolerance_) ? value_ : ( (value_ > tolerance_) ? value_ : 0.0  );
  };

};

//! @copydoc wmesh_math
template <> struct wmesh_math<float>
{
  static constexpr unsigned int num_digits = 7;
  static constexpr float zero = 0.0f;
  static constexpr float epsilon = 1.19209290e-07f;
  //! @copydoc wmesh_math::pow(const real_t&x_,const real_t&y_)
  static constexpr float pow(const float&x_,const float&y_) noexcept
  {
    return powf(x_,y_);
  };
  //! @copydoc wmesh_math::sqrt(const real_t&x_)
  static constexpr float sqrt(const float&x_) noexcept
  {
    return sqrtf(x_);
  };
  //! @copydoc wmesh_math::abs(const real_t&x_)
  static constexpr float abs(const float&x_) noexcept
  {
    return  (x_<0.0) ? -x_ : x_;
  };

  //! @copydoc wmesh_math::max(const real_t&,const real_t&)
  static constexpr float max(const float&a_,const float&b_) noexcept
  {
    return (a_>=b_) ? a_ : b_;
  };
  //! @copydoc wmesh_math::min(const real_t&,const real_t&)
  static constexpr float min(const float&a_,const float&b_) noexcept
  {
    return (a_<=b_) ? a_ : b_;
  };

  //! @copydoc wmesh_math::filter(const real_t&,const real_t&)
  static constexpr float filter(const float &value_,const float&tolerance_) noexcept
  {
    return (value_ < -tolerance_) ? value_ : ( (value_ > tolerance_) ? value_ : 0.0f  );
  };

};

//! @copydoc wmesh_math
template <> struct wmesh_math<long double>
{
  static constexpr unsigned int num_digits = 17;
  static constexpr long double zero = 0.0L;
  static constexpr long double epsilon = 1.08420217248550443401e-19L;
  //! @copydoc wmesh_math::pow(const real_t&x_,const real_t&y_)
  static constexpr long double pow(const long double&x_,const long double&y_) noexcept
  {
    return powl(x_,y_);
  };
  //! @copydoc wmesh_math::sqrt(const real_t&x_)
  static constexpr long double sqrt(const long double&x_) noexcept
  {
    return sqrtl(x_);
  };
  //! @copydoc wmesh_math::abs(const real_t&x_)
  static constexpr long double abs(const long double&x_) noexcept
  {
    return  (x_<0.0) ? -x_ : x_;
  };

  //! @copydoc wmesh_math::max(const real_t&,const real_t&)
  static constexpr long double max(const long double&a_,const long double&b_) noexcept
  {
    return (a_>=b_) ? a_ : b_;
  };
  //! @copydoc wmesh_math::min(const real_t&,const real_t&)
  static constexpr long double min(const long double&a_,const long double&b_) noexcept
  {
    return (a_<=b_) ? a_ : b_;
  };

  //! @copydoc wmesh_math::filter(const real_t&,const real_t&)
  static constexpr long double filter(const long double &value_,const long double&tolerance_) noexcept
  {
    return (value_ < -tolerance_) ? value_ : ( (value_ > tolerance_) ? value_ : 0.0L  );
  };

};
