//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_math_details_matrix_inversion_2x2_hpp__
#define __pinocchio_math_details_matrix_inversion_2x2_hpp__

#include "pinocchio/math/matrix.hpp"

namespace pinocchio
{
  namespace internal
  {
    template<>
    struct MatrixInversionCodeGeneratedImpl<2>
    {
      template<typename M1, typename M2>
      static void
      run(const Eigen::MatrixBase<M1> & matrix, const Eigen::MatrixBase<M2> & matrix_inverse_)
      {
        typedef typename M1::Scalar Scalar;

        assert(is_symmetric(matrix));

        const auto & input_vec = matrix.reshaped();
        auto & matrix_inverse = matrix_inverse_.const_cast_derived();
        auto output_vec = matrix_inverse.reshaped();

        Scalar a0, a1, a2, a3;
        a0 = input_vec[0];
        a1 = input_vec[3];
        a2 = input_vec[2];
        a2 = (a2 / a1);
        a3 = math::square(a2);
        a3 = (a1 * a3);
        a0 = (a0 - a3);
        a3 = (1. / a0);
        output_vec[0] = a3;
        a0 = (a2 * a3);
        a3 = (-a0);
        output_vec[1] = a3;
        output_vec[2] = a3;
        a1 = (1. / a1);
        a2 = (a2 * a0);
        a1 = (a1 + a2);
        output_vec[3] = a1;
      }
    };
  } // namespace internal
} // namespace pinocchio

#endif // ifndef __pinocchio_math_details_matrix_inversion_2x2_hpp__
