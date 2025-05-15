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

        assert(is_symmetric(matrix));

        const auto & input_vec = matrix.reshaped();
        auto & matrix_inverse = matrix_inverse_.const_cast_derived();
        auto output_vec = matrix_inverse.reshaped();

        Scalar a00, a01, a02, a03, a04, a05, a06, a07, a08, a09, a10;
        a00 = input_vec[0];
        a01 = input_vec[3];
        a02 = input_vec[6];
        a03 = input_vec[8];
        a03 = math::sqrt(a03);
        a02 = (a02 / a03);
        a04 = input_vec[7];
        a04 = (a04 / a03);
        a05 = (a02 * a04);
        a01 = (a01 - a05);
        a05 = input_vec[4];
        a06 = math::square(a04);
        a05 = (a05 - a06);
        a05 = math::sqrt(a05);
        a01 = (a01 / a05);
        a06 = math::square(a01);
        a07 = math::square(a02);
        a06 = (a06 + a07);
        a00 = (a00 - a06);
        a00 = math::sqrt(a00);
        a06 = (1. / a00);
        a06 = (a06 / a00);
        output_vec[0] = a06;
        a06 = (a01 / a05);
        a06 = (a06 / a00);
        a06 = (a06 / a00);
        a07 = (-a06);
        output_vec[1] = a07;
        a08 = (a02 / a03);
        a09 = (a04 / a03);
        a09 = (a09 / a05);
        a10 = (a01 * a09);
        a08 = (a08 - a10);
        a08 = (a08 / a00);
        a08 = (a08 / a00);
        a00 = (-a08);
        output_vec[2] = a00;
        output_vec[3] = a07;
        a07 = (1. / a05);
        a06 = (a01 * a06);
        a07 = (a07 + a06);
        a07 = (a07 / a05);
        output_vec[4] = a07;
        a01 = (a01 * a08);
        a01 = (a01 - a09);
        a01 = (a01 / a05);
        output_vec[5] = a01;
        output_vec[6] = a00;
        output_vec[7] = a01;
        a00 = (1. / a03);
        a04 = (a04 * a01);
        a02 = (a02 * a08);
        a04 = (a04 - a02);
        a00 = (a00 - a04);
        a00 = (a00 / a03);
        output_vec[8] = a00;
      }
    };
  } // namespace internal
} // namespace pinocchio

#endif // ifndef __pinocchio_math_details_matrix_inversion_3x3_hpp__
