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
        Scalar a24, a25;
        a00 = input_vec[0];
        a01 = input_vec[6];
        a02 = input_vec[12];
        a03 = input_vec[18];
        a04 = input_vec[24];
        a05 = input_vec[23];
        a05 = (a05 / a04);
        a06 = math::square(a05);
        a06 = (a04 * a06);
        a03 = (a03 - a06);
        a06 = input_vec[17];
        a07 = input_vec[22];
        a08 = (a07 * a05);
        a06 = (a06 - a08);
        a06 = (a06 / a03);
        a08 = math::square(a06);
        a08 = (a03 * a08);
        a07 = (a07 / a04);
        a09 = math::square(a07);
        a09 = (a04 * a09);
        a08 = (a08 + a09);
        a02 = (a02 - a08);
        a08 = input_vec[11];
        a09 = input_vec[16];
        a10 = input_vec[21];
        a11 = (a10 * a05);
        a09 = (a09 - a11);
        a11 = (a09 * a06);
        a12 = (a10 * a07);
        a11 = (a11 + a12);
        a08 = (a08 - a11);
        a08 = (a08 / a02);
        a11 = math::square(a08);
        a11 = (a02 * a11);
        a09 = (a09 / a03);
        a12 = math::square(a09);
        a12 = (a03 * a12);
        a11 = (a11 + a12);
        a10 = (a10 / a04);
        a12 = math::square(a10);
        a12 = (a04 * a12);
        a11 = (a11 + a12);
        a01 = (a01 - a11);
        a11 = input_vec[5];
        a12 = input_vec[10];
        a13 = input_vec[15];
        a14 = input_vec[20];
        a15 = (a14 * a05);
        a13 = (a13 - a15);
        a15 = (a13 * a06);
        a16 = (a14 * a07);
        a15 = (a15 + a16);
        a12 = (a12 - a15);
        a15 = (a12 * a08);
        a16 = (a13 * a09);
        a15 = (a15 + a16);
        a16 = (a14 * a10);
        a15 = (a15 + a16);
        a11 = (a11 - a15);
        a11 = (a11 / a01);
        a15 = math::square(a11);
        a15 = (a01 * a15);
        a12 = (a12 / a02);
        a16 = math::square(a12);
        a16 = (a02 * a16);
        a15 = (a15 + a16);
        a13 = (a13 / a03);
        a16 = math::square(a13);
        a16 = (a03 * a16);
        a15 = (a15 + a16);
        a14 = (a14 / a04);
        a16 = math::square(a14);
        a16 = (a04 * a16);
        a15 = (a15 + a16);
        a00 = (a00 - a15);
        a15 = (1. / a00);
        output_vec[0] = a15;
        a15 = (a11 / a00);
        a16 = (-a15);
        output_vec[1] = a16;
        a17 = (a11 * a08);
        a17 = (a12 - a17);
        a17 = (a17 / a00);
        a18 = (-a17);
        output_vec[2] = a18;
        a19 = (a08 * a06);
        a19 = (a09 - a19);
        a20 = (a11 * a19);
        a21 = (a12 * a06);
        a20 = (a20 + a21);
        a20 = (a13 - a20);
        a20 = (a20 / a00);
        a21 = (-a20);
        output_vec[3] = a21;
        a22 = (a06 * a05);
        a22 = (a07 - a22);
        a23 = (a08 * a22);
        a24 = (a09 * a05);
        a23 = (a23 + a24);
        a23 = (a10 - a23);
        a24 = (a11 * a23);
        a25 = (a12 * a22);
        a24 = (a24 + a25);
        a25 = (a13 * a05);
        a24 = (a24 + a25);
        a24 = (a14 - a24);
        a24 = (a24 / a00);
        a00 = (-a24);
        output_vec[4] = a00;
        output_vec[5] = a16;
        a16 = (1. / a01);
        a15 = (a11 * a15);
        a16 = (a16 + a15);
        output_vec[6] = a16;
        a16 = (a11 * a17);
        a15 = (a08 / a01);
        a16 = (a16 - a15);
        output_vec[7] = a16;
        a15 = (a11 * a20);
        a19 = (a19 / a01);
        a15 = (a15 - a19);
        output_vec[8] = a15;
        a11 = (a11 * a24);
        a23 = (a23 / a01);
        a11 = (a11 - a23);
        output_vec[9] = a11;
        output_vec[10] = a18;
        output_vec[11] = a16;
        a18 = (1. / a02);
        a16 = (a08 * a16);
        a17 = (a12 * a17);
        a16 = (a16 - a17);
        a18 = (a18 - a16);
        output_vec[12] = a18;
        a18 = (a06 / a02);
        a16 = (a08 * a15);
        a17 = (a12 * a20);
        a16 = (a16 - a17);
        a18 = (a18 + a16);
        a16 = (-a18);
        output_vec[13] = a16;
        a22 = (a22 / a02);
        a08 = (a08 * a11);
        a12 = (a12 * a24);
        a08 = (a08 - a12);
        a22 = (a22 + a08);
        a08 = (-a22);
        output_vec[14] = a08;
        output_vec[15] = a21;
        output_vec[16] = a15;
        output_vec[17] = a16;
        a16 = (1. / a03);
        a15 = (a09 * a15);
        a20 = (a13 * a20);
        a15 = (a15 - a20);
        a18 = (a06 * a18);
        a15 = (a15 - a18);
        a16 = (a16 - a15);
        output_vec[18] = a16;
        a03 = (a05 / a03);
        a09 = (a09 * a11);
        a13 = (a13 * a24);
        a09 = (a09 - a13);
        a06 = (a06 * a22);
        a09 = (a09 - a06);
        a03 = (a03 + a09);
        a09 = (-a03);
        output_vec[19] = a09;
        output_vec[20] = a00;
        output_vec[21] = a11;
        output_vec[22] = a08;
        output_vec[23] = a09;
        a04 = (1. / a04);
        a10 = (a10 * a11);
        a14 = (a14 * a24);
        a10 = (a10 - a14);
        a07 = (a07 * a22);
        a10 = (a10 - a07);
        a05 = (a05 * a03);
        a10 = (a10 - a05);
        a04 = (a04 - a10);
        output_vec[24] = a04;
      }
    };
  } // namespace internal
} // namespace pinocchio

#endif // ifndef __pinocchio_math_details_matrix_inversion_5x5_hpp__
