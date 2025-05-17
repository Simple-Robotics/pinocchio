//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_math_details_matrix_inversion_3x3_hpp__
#define __pinocchio_math_details_matrix_inversion_3x3_hpp__

#include "pinocchio/math/matrix.hpp"

namespace pinocchio
{
  namespace internal
  {
    template<>
    struct MatrixInversionCodeGeneratedImpl<3>
    {
      template<typename M1, typename M2>
      static EIGEN_STRONG_INLINE void
      run(const Eigen::MatrixBase<M1> & matrix, const Eigen::MatrixBase<M2> & matrix_inverse_)
      {
        typedef typename M1::Scalar Scalar;

        assert(is_symmetric(matrix, math::sqrt(dummy_precision<Scalar>())));

        const auto & input_vec = matrix.reshaped();
        auto & matrix_inverse = matrix_inverse_.const_cast_derived();
        auto output_vec = matrix_inverse.reshaped();
      }
    };
  } // namespace internal
} // namespace pinocchio

#endif // ifndef __pinocchio_math_details_matrix_inversion_3x3_hpp__
