//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_math_details_matrix_inversion_4x4_hpp__
#define __pinocchio_math_details_matrix_inversion_4x4_hpp__

#include "pinocchio/math/matrix.hpp"

namespace pinocchio
{
  namespace internal
  {
    template<>
    struct MatrixInversionCodeGeneratedImpl<4>
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

        Scalar a00, a01, a02, a03, a04, a05, a06, a07, a08, a09, a10, a11;
        Scalar a12, a13, a14, a15, a16;
        a00 = input_vec[0];
        a01 = input_vec[5];
        a02 = input_vec[10];
        a03 = input_vec[15];
        a04 = input_vec[14];
        a04 = (a04 / a03);
        a05 = math::square(a04);
        a05 = (a03 * a05);
        a02 = (a02 - a05);
        a05 = input_vec[9];
        a06 = input_vec[13];
        a07 = (a06 * a04);
        a05 = (a05 - a07);
        a05 = (a05 / a02);
        a07 = math::square(a05);
        a07 = (a02 * a07);
        a06 = (a06 / a03);
        a08 = math::square(a06);
        a08 = (a03 * a08);
        a07 = (a07 + a08);
        a01 = (a01 - a07);
        a07 = input_vec[4];
        a08 = input_vec[8];
        a09 = input_vec[12];
        a10 = (a09 * a04);
        a08 = (a08 - a10);
        a10 = (a08 * a05);
        a11 = (a09 * a06);
        a10 = (a10 + a11);
        a07 = (a07 - a10);
        a07 = (a07 / a01);
        a10 = math::square(a07);
        a10 = (a01 * a10);
        a08 = (a08 / a02);
        a11 = math::square(a08);
        a11 = (a02 * a11);
        a10 = (a10 + a11);
        a09 = (a09 / a03);
        a11 = math::square(a09);
        a11 = (a03 * a11);
        a10 = (a10 + a11);
        a00 = (a00 - a10);
        a10 = (1. / a00);
        output_vec[0] = a10;
        a10 = (a07 / a00);
        a11 = (-a10);
        output_vec[1] = a11;
        a12 = (a07 * a05);
        a12 = (a08 - a12);
        a12 = (a12 / a00);
        a13 = (-a12);
        output_vec[2] = a13;
        a14 = (a05 * a04);
        a14 = (a06 - a14);
        a15 = (a07 * a14);
        a16 = (a08 * a04);
        a15 = (a15 + a16);
        a15 = (a09 - a15);
        a15 = (a15 / a00);
        a00 = (-a15);
        output_vec[3] = a00;
        output_vec[4] = a11;
        a11 = (1. / a01);
        a10 = (a07 * a10);
        a11 = (a11 + a10);
        output_vec[5] = a11;
        a11 = (a07 * a12);
        a10 = (a05 / a01);
        a11 = (a11 - a10);
        output_vec[6] = a11;
        a07 = (a07 * a15);
        a14 = (a14 / a01);
        a07 = (a07 - a14);
        output_vec[7] = a07;
        output_vec[8] = a13;
        output_vec[9] = a11;
        a13 = (1. / a02);
        a11 = (a05 * a11);
        a12 = (a08 * a12);
        a11 = (a11 - a12);
        a13 = (a13 - a11);
        output_vec[10] = a13;
        a02 = (a04 / a02);
        a05 = (a05 * a07);
        a08 = (a08 * a15);
        a05 = (a05 - a08);
        a02 = (a02 + a05);
        a05 = (-a02);
        output_vec[11] = a05;
        output_vec[12] = a00;
        output_vec[13] = a07;
        output_vec[14] = a05;
        a03 = (1. / a03);
        a06 = (a06 * a07);
        a09 = (a09 * a15);
        a06 = (a06 - a09);
        a04 = (a04 * a02);
        a06 = (a06 - a04);
        a03 = (a03 - a06);
        output_vec[15] = a03;
      }
    };
  } // namespace internal
} // namespace pinocchio

#endif // ifndef __pinocchio_math_details_matrix_inversion_4x4_hpp__
