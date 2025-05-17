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
      static void
      run(const Eigen::MatrixBase<M1> & matrix, const Eigen::MatrixBase<M2> & matrix_inverse_)
      {
        typedef typename M1::Scalar Scalar;

        assert(is_symmetric(matrix));

        const auto & input_vec = matrix.reshaped();
        auto & matrix_inverse = matrix_inverse_.const_cast_derived();
        auto output_vec = matrix_inverse.reshaped();

        Scalar a0, a1, a2, a3, a4, a5, a6, a7, a8;
        a0 = input_vec[0];
        a1 = input_vec[4];
        a2 = input_vec[8];
        a3 = input_vec[7];
        a3 = (a3 / a2);
        a4 = math::square(a3);
        a4 = (a2 * a4);
        a1 = (a1 - a4);
        a4 = input_vec[3];
        a5 = input_vec[6];
        a6 = (a5 * a3);
        a4 = (a4 - a6);
        a4 = (a4 / a1);
        a6 = math::square(a4);
        a6 = (a1 * a6);
        a5 = (a5 / a2);
        a7 = math::square(a5);
        a7 = (a2 * a7);
        a6 = (a6 + a7);
        a0 = (a0 - a6);
        a6 = (1. / a0);
        output_vec[0] = a6;
        a6 = (a4 / a0);
        a7 = (-a6);
        output_vec[1] = a7;
        a8 = (a4 * a3);
        a8 = (a5 - a8);
        a8 = (a8 / a0);
        a0 = (-a8);
        output_vec[2] = a0;
        output_vec[3] = a7;
        a7 = (1. / a1);
        a6 = (a4 * a6);
        a7 = (a7 + a6);
        output_vec[4] = a7;
        a4 = (a4 * a8);
        a1 = (a3 / a1);
        a4 = (a4 - a1);
        output_vec[5] = a4;
        output_vec[6] = a0;
        output_vec[7] = a4;
        a2 = (1. / a2);
        a3 = (a3 * a4);
        a5 = (a5 * a8);
        a3 = (a3 - a5);
        a2 = (a2 - a3);
        output_vec[8] = a2;
      }
    };
  } // namespace internal
} // namespace pinocchio

#endif // ifndef __pinocchio_math_details_matrix_inversion_3x3_hpp__
