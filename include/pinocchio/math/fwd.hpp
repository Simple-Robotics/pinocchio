//
// Copyright (c) 2016-2018 CNRS
// Copyright (c) 2018-2025 INRIA
//

#ifndef __pinocchio_math_fwd_hpp__
#define __pinocchio_math_fwd_hpp__

#include "pinocchio/fwd.hpp"

#include <math.h>
#include <boost/math/constants/constants.hpp>
#include <type_traits>

namespace pinocchio
{

  template<typename T>
  struct EigenMatrixExpression;

  template<typename T>
  struct is_floating_point : ::std::is_floating_point<T>
  {
  };

  ///
  /// \brief Returns the value of PI according to the template parameters Scalar
  ///
  /// \tparam Scalar The scalar type of the return pi value
  ///
  template<typename Scalar>
  const Scalar PI()
  {
    return boost::math::constants::pi<Scalar>();
  }

  ///  \brief Foward declaration of TaylorSeriesExpansion.
  template<typename Scalar>
  struct TaylorSeriesExpansion;

  namespace math
  {

#define PINOCCHIO_OVERLOAD_MATH_UNARY_OPERATOR(name)                                               \
  template<typename Scalar>                                                                        \
  Scalar name(const Scalar & value)                                                                \
  {                                                                                                \
    using std::name;                                                                               \
    return name(value);                                                                            \
  }

#define PINOCCHIO_OVERLOAD_MATH_BINARY_OPERATOR(name)                                              \
  namespace internal                                                                               \
  {                                                                                                \
    template<typename T1, typename T2>                                                             \
    struct return_type_##name                                                                      \
    {                                                                                              \
      typedef T1 type;                                                                             \
    };                                                                                             \
    template<typename T1, typename T2>                                                             \
    struct call_##name                                                                             \
    {                                                                                              \
      static inline typename return_type_##name<T1, T2>::type run(const T1 & a, const T2 & b)      \
      {                                                                                            \
        using std::name;                                                                           \
        return name(a, b);                                                                         \
      }                                                                                            \
    };                                                                                             \
  }                                                                                                \
  template<typename T1, typename T2>                                                               \
  inline typename internal::return_type_##name<T1, T2>::type name(const T1 & a, const T2 & b)      \
  {                                                                                                \
    return internal::call_##name<T1, T2>::run(a, b);                                               \
  }

    PINOCCHIO_OVERLOAD_MATH_UNARY_OPERATOR(fabs)
    PINOCCHIO_OVERLOAD_MATH_UNARY_OPERATOR(sqrt)
    PINOCCHIO_OVERLOAD_MATH_UNARY_OPERATOR(atan)
    PINOCCHIO_OVERLOAD_MATH_UNARY_OPERATOR(acos)
    PINOCCHIO_OVERLOAD_MATH_UNARY_OPERATOR(asin)
    PINOCCHIO_OVERLOAD_MATH_UNARY_OPERATOR(cos)
    PINOCCHIO_OVERLOAD_MATH_UNARY_OPERATOR(sin)
    PINOCCHIO_OVERLOAD_MATH_UNARY_OPERATOR(log)
    PINOCCHIO_OVERLOAD_MATH_UNARY_OPERATOR(exp)

    PINOCCHIO_OVERLOAD_MATH_BINARY_OPERATOR(pow)
    PINOCCHIO_OVERLOAD_MATH_BINARY_OPERATOR(min)
    PINOCCHIO_OVERLOAD_MATH_BINARY_OPERATOR(max)
    PINOCCHIO_OVERLOAD_MATH_BINARY_OPERATOR(atan2)

    template<typename T>
    inline T square(const T & value)
    {
      return value * value;
    }
  } // namespace math

  /// \brief Shortcut for calling ::Eigen::NumTraits<T>::dummy_precision()
  template<typename T>
  inline EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR T dummy_precision()
  {
    return ::Eigen::NumTraits<T>::dummy_precision();
  }

} // namespace pinocchio

#endif // #ifndef __pinocchio_math_fwd_hpp__
