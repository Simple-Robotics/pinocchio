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
      static EIGEN_STRONG_INLINE void
      run(const Eigen::MatrixBase<M1> & matrix, const Eigen::MatrixBase<M2> & matrix_inverse_)
      {
        typedef typename M1::Scalar Scalar;

        assert(is_symmetric(matrix));

        const auto & input_vec = matrix.reshaped();
        auto & matrix_inverse = matrix_inverse_.const_cast_derived();
        auto output_vec = matrix_inverse.reshaped();

        Scalar a00, a01, a02, a03, a04, a05, a06, a07, a08, a09, a10, a11;
        Scalar a12, a13, a14, a15, a16, a17, a18, a19;

        a00 = input_vec[0];
        a01 = input_vec[4];
        a02 = input_vec[8];
        a03 = input_vec[12];
        a04 = input_vec[15];
        a04 = math::sqrt(a04);
        a04 = (Scalar(1) / a04);
        a03 = (a03 * a04);
        a05 = input_vec[14];
        a05 = (a05 * a04);
        a06 = (a03 * a05);
        a02 = (a02 - a06);
        a06 = input_vec[10];
        a07 = math::square(a05);
        a06 = (a06 - a07);
        a06 = math::sqrt(a06);
        a06 = (Scalar(1) / a06);
        a02 = (a02 * a06);
        a07 = input_vec[9];
        a08 = input_vec[13];
        a08 = (a08 * a04);
        a09 = (a08 * a05);
        a07 = (a07 - a09);
        a07 = (a07 * a06);
        a09 = (a02 * a07);
        a10 = (a03 * a08);
        a09 = (a09 + a10);
        a01 = (a01 - a09);
        a09 = input_vec[5];
        a10 = math::square(a07);
        a11 = math::square(a08);
        a10 = (a10 + a11);
        a09 = (a09 - a10);
        a09 = math::sqrt(a09);
        a09 = (Scalar(1) / a09);
        a01 = (a01 * a09);
        a10 = math::square(a01);
        a11 = math::square(a02);
        a10 = (a10 + a11);
        a11 = math::square(a03);
        a10 = (a10 + a11);
        a00 = (a00 - a10);
        a00 = math::sqrt(a00);
        a00 = (Scalar(1) / a00);
        a10 = (a00 * a00);
        output_vec[0] = a10;
        a10 = (a01 * a09);
        a10 = (a10 * a00);
        a10 = (a10 * a00);
        a11 = (-a10);
        output_vec[1] = a11;
        a12 = (a02 * a06);
        a13 = (a07 * a06);
        a13 = (a13 * a09);
        a14 = (a01 * a13);
        a12 = (a12 - a14);
        a12 = (a12 * a00);
        a12 = (a12 * a00);
        a14 = (-a12);
        output_vec[2] = a14;
        a15 = (a03 * a04);
        a16 = (a08 * a04);
        a17 = (a05 * a04);
        a17 = (a17 * a06);
        a18 = (a07 * a17);
        a16 = (a16 - a18);
        a16 = (a16 * a09);
        a18 = (a01 * a16);
        a19 = (a02 * a17);
        a18 = (a18 + a19);
        a15 = (a15 - a18);
        a15 = (a15 * a00);
        a15 = (a15 * a00);
        a00 = (-a15);
        output_vec[3] = a00;
        output_vec[4] = a11;
        a11 = a09;
        a10 = (a01 * a10);
        a11 = (a11 + a10);
        a11 = (a11 * a09);
        output_vec[5] = a11;
        a11 = (a01 * a12);
        a11 = (a11 - a13);
        a11 = (a11 * a09);
        output_vec[6] = a11;
        a01 = (a01 * a15);
        a01 = (a01 - a16);
        a01 = (a01 * a09);
        output_vec[7] = a01;
        output_vec[8] = a14;
        output_vec[9] = a11;
        a14 = a06;
        a11 = (a07 * a11);
        a12 = (a02 * a12);
        a11 = (a11 - a12);
        a14 = (a14 - a11);
        a14 = (a14 * a06);
        output_vec[10] = a14;
        a07 = (a07 * a01);
        a02 = (a02 * a15);
        a07 = (a07 - a02);
        a17 = (a17 + a07);
        a17 = (a17 * a06);
        a06 = (-a17);
        output_vec[11] = a06;
        output_vec[12] = a00;
        output_vec[13] = a01;
        output_vec[14] = a06;
        a06 = (Scalar(1) * a04);
        a08 = (a08 * a01);
        a03 = (a03 * a15);
        a08 = (a08 - a03);
        a05 = (a05 * a17);
        a08 = (a08 - a05);
        a06 = (a06 - a08);
        a06 = (a06 * a04);
        output_vec[15] = a06;
      }
    };
  } // namespace internal
} // namespace pinocchio

#endif // ifndef __pinocchio_math_details_matrix_inversion_4x4_hpp__
