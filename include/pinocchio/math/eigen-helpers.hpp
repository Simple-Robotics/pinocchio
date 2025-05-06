//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_math_eigen_helpers_hpp__
#define __pinocchio_math_eigen_helpers_hpp__

#include "pinocchio/math/fwd.hpp"

namespace pinocchio
{

#define PINOCCHIO_EIGEN_HELPER_FUNC(method)                                                        \
  template<typename Matrix>                                                                        \
  void method(const Eigen::MatrixBase<Matrix> & mat)                                               \
  {                                                                                                \
    mat.const_cast_derived().method();                                                             \
  }

  PINOCCHIO_EIGEN_HELPER_FUNC(setZero);
  PINOCCHIO_EIGEN_HELPER_FUNC(setOnes);
  PINOCCHIO_EIGEN_HELPER_FUNC(setIdentity);

#undef PINOCCHIO_EIGEN_HELPER_FUNC

} // namespace pinocchio

#endif // ifndef __pinocchio_math_eigen_helpers_hpp__
