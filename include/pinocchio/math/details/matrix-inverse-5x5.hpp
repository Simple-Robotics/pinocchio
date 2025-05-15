//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_math_details_matrix_inversion_5x5_hpp__
#define __pinocchio_math_details_matrix_inversion_5x5_hpp__

#include "pinocchio/math/matrix.hpp"

namespace pinocchio
{
  namespace internal
  {
    template<>
    struct MatrixInversionCodeGeneratedImpl<5>
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
        Scalar a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23;
        Scalar a24, a25, a26, a27, a28, a29;
        a00 = input_vec[0];
        a01 = input_vec[5];
        a02 = input_vec[10];
        a03 = input_vec[15];
        a04 = input_vec[20];
        a05 = input_vec[24];
        a05 = math::sqrt(a05);
        a04 = (a04 / a05);
        a06 = input_vec[23];
        a06 = (a06 / a05);
        a07 = (a04 * a06);
        a03 = (a03 - a07);
        a07 = input_vec[18];
        a08 = math::square(a06);
        a07 = (a07 - a08);
        a07 = math::sqrt(a07);
        a03 = (a03 / a07);
        a08 = input_vec[17];
        a09 = input_vec[22];
        a09 = (a09 / a05);
        a10 = (a09 * a06);
        a08 = (a08 - a10);
        a08 = (a08 / a07);
        a10 = (a03 * a08);
        a11 = (a04 * a09);
        a10 = (a10 + a11);
        a02 = (a02 - a10);
        a10 = input_vec[12];
        a11 = math::square(a08);
        a12 = math::square(a09);
        a11 = (a11 + a12);
        a10 = (a10 - a11);
        a10 = math::sqrt(a10);
        a02 = (a02 / a10);
        a11 = input_vec[11];
        a12 = input_vec[16];
        a13 = input_vec[21];
        a13 = (a13 / a05);
        a14 = (a13 * a06);
        a12 = (a12 - a14);
        a12 = (a12 / a07);
        a14 = (a12 * a08);
        a15 = (a13 * a09);
        a14 = (a14 + a15);
        a11 = (a11 - a14);
        a11 = (a11 / a10);
        a14 = (a02 * a11);
        a15 = (a03 * a12);
        a14 = (a14 + a15);
        a15 = (a04 * a13);
        a14 = (a14 + a15);
        a01 = (a01 - a14);
        a14 = input_vec[6];
        a15 = math::square(a11);
        a16 = math::square(a12);
        a15 = (a15 + a16);
        a16 = math::square(a13);
        a15 = (a15 + a16);
        a14 = (a14 - a15);
        a14 = math::sqrt(a14);
        a01 = (a01 / a14);
        a15 = math::square(a01);
        a16 = math::square(a02);
        a15 = (a15 + a16);
        a16 = math::square(a03);
        a15 = (a15 + a16);
        a16 = math::square(a04);
        a15 = (a15 + a16);
        a00 = (a00 - a15);
        a00 = math::sqrt(a00);
        a15 = (1. / a00);
        a15 = (a15 / a00);
        output_vec[0] = a15;
        a15 = (a01 / a14);
        a15 = (a15 / a00);
        a15 = (a15 / a00);
        a16 = (-a15);
        output_vec[1] = a16;
        a17 = (a02 / a10);
        a18 = (a11 / a10);
        a18 = (a18 / a14);
        a19 = (a01 * a18);
        a17 = (a17 - a19);
        a17 = (a17 / a00);
        a17 = (a17 / a00);
        a19 = (-a17);
        output_vec[2] = a19;
        a20 = (a03 / a07);
        a21 = (a12 / a07);
        a22 = (a08 / a07);
        a22 = (a22 / a10);
        a23 = (a11 * a22);
        a21 = (a21 - a23);
        a21 = (a21 / a14);
        a23 = (a01 * a21);
        a24 = (a02 * a22);
        a23 = (a23 + a24);
        a20 = (a20 - a23);
        a20 = (a20 / a00);
        a20 = (a20 / a00);
        a23 = (-a20);
        output_vec[3] = a23;
        a24 = (a04 / a05);
        a25 = (a13 / a05);
        a26 = (a09 / a05);
        a27 = (a06 / a05);
        a27 = (a27 / a07);
        a28 = (a08 * a27);
        a26 = (a26 - a28);
        a26 = (a26 / a10);
        a28 = (a11 * a26);
        a29 = (a12 * a27);
        a28 = (a28 + a29);
        a25 = (a25 - a28);
        a25 = (a25 / a14);
        a28 = (a01 * a25);
        a29 = (a02 * a26);
        a28 = (a28 + a29);
        a29 = (a03 * a27);
        a28 = (a28 + a29);
        a24 = (a24 - a28);
        a24 = (a24 / a00);
        a24 = (a24 / a00);
        a00 = (-a24);
        output_vec[4] = a00;
        output_vec[5] = a16;
        a16 = (1. / a14);
        a15 = (a01 * a15);
        a16 = (a16 + a15);
        a16 = (a16 / a14);
        output_vec[6] = a16;
        a16 = (a01 * a17);
        a16 = (a16 - a18);
        a16 = (a16 / a14);
        output_vec[7] = a16;
        a18 = (a01 * a20);
        a18 = (a18 - a21);
        a18 = (a18 / a14);
        output_vec[8] = a18;
        a01 = (a01 * a24);
        a01 = (a01 - a25);
        a01 = (a01 / a14);
        output_vec[9] = a01;
        output_vec[10] = a19;
        output_vec[11] = a16;
        a19 = (1. / a10);
        a16 = (a11 * a16);
        a17 = (a02 * a17);
        a16 = (a16 - a17);
        a19 = (a19 - a16);
        a19 = (a19 / a10);
        output_vec[12] = a19;
        a19 = (a11 * a18);
        a16 = (a02 * a20);
        a19 = (a19 - a16);
        a22 = (a22 + a19);
        a22 = (a22 / a10);
        a19 = (-a22);
        output_vec[13] = a19;
        a11 = (a11 * a01);
        a02 = (a02 * a24);
        a11 = (a11 - a02);
        a26 = (a26 + a11);
        a26 = (a26 / a10);
        a10 = (-a26);
        output_vec[14] = a10;
        output_vec[15] = a23;
        output_vec[16] = a18;
        output_vec[17] = a19;
        a19 = (1. / a07);
        a18 = (a12 * a18);
        a20 = (a03 * a20);
        a18 = (a18 - a20);
        a22 = (a08 * a22);
        a18 = (a18 - a22);
        a19 = (a19 - a18);
        a19 = (a19 / a07);
        output_vec[18] = a19;
        a12 = (a12 * a01);
        a03 = (a03 * a24);
        a12 = (a12 - a03);
        a08 = (a08 * a26);
        a12 = (a12 - a08);
        a27 = (a27 + a12);
        a27 = (a27 / a07);
        a07 = (-a27);
        output_vec[19] = a07;
        output_vec[20] = a00;
        output_vec[21] = a01;
        output_vec[22] = a10;
        output_vec[23] = a07;
        a07 = (1. / a05);
        a13 = (a13 * a01);
        a04 = (a04 * a24);
        a13 = (a13 - a04);
        a09 = (a09 * a26);
        a13 = (a13 - a09);
        a06 = (a06 * a27);
        a13 = (a13 - a06);
        a07 = (a07 - a13);
        a07 = (a07 / a05);
        output_vec[24] = a07;
      }
    };
  } // namespace internal
} // namespace pinocchio

#endif // ifndef __pinocchio_math_details_matrix_inversion_5x5_hpp__
