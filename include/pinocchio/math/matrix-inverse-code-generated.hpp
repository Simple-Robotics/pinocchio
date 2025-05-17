//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_math_matrix_inverse_code_generated_hpp__
#define __pinocchio_math_matrix_inverse_code_generated_hpp__

#include "pinocchio/math/matrix.hpp"

namespace pinocchio
{
  namespace internal
  {
    template<int RowsAtCompileTime, int ColsAtCompileTime = RowsAtCompileTime>
    struct MatrixInversionCodeGeneratedImpl
    {
      template<typename M1, typename M2>
      static EIGEN_STRONG_INLINE void run(
        const Eigen::MatrixBase<M1> & /*matrix*/, const Eigen::MatrixBase<M2> & /*matrix_inverse*/)
      {
        // static_assert(false, "Not implemented.");
        assert(false && "Not implemented.");
      }
    };
  } // namespace internal
} // namespace pinocchio

#include "pinocchio/math/details/matrix-inverse-1x1.hpp"
#include "pinocchio/math/details/matrix-inverse-2x2.hpp"
#include "pinocchio/math/details/matrix-inverse-3x3.hpp"
#include "pinocchio/math/details/matrix-inverse-4x4.hpp"
#include "pinocchio/math/details/matrix-inverse-5x5.hpp"
#include "pinocchio/math/details/matrix-inverse-6x6.hpp"
#include "pinocchio/math/details/matrix-inverse-7x7.hpp"
#include "pinocchio/math/details/matrix-inverse-8x8.hpp"
#include "pinocchio/math/details/matrix-inverse-9x9.hpp"
#include "pinocchio/math/details/matrix-inverse-10x10.hpp"
#include "pinocchio/math/details/matrix-inverse-11x11.hpp"
#include "pinocchio/math/details/matrix-inverse-12x12.hpp"

#endif // ifndef __pinocchio_math_matrix_inverse_code_generated_hpp__
