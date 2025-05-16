//
// Copyright (c) 2019-2025 INRIA
//

#ifndef __pinocchio_math_matrix_inverse_hpp__
#define __pinocchio_math_matrix_inverse_hpp__

#include "pinocchio/math/matrix.hpp"

namespace pinocchio
{
  namespace internal
  {
    struct MatrixInversionDefaultImpl
    {
      template<typename M1, typename M2>
      static EIGEN_STRONG_INLINE void
      run(const Eigen::MatrixBase<M1> & matrix, const Eigen::MatrixBase<M2> & matrix_inverse)
      {
        matrix_inverse.const_cast_derived().noalias() = matrix.inverse();
      }
    };

    struct MatrixInversionDynamicMatrixImpl
    {
      template<typename M1, typename M2>
      static EIGEN_STRONG_INLINE void
      run(const Eigen::MatrixBase<M1> & matrix, const Eigen::MatrixBase<M2> & matrix_inverse)
      {
        assert(is_symmetric(matrix));
        auto & matrix_inverse_ = matrix_inverse.const_cast_derived();
        matrix_inverse_.setIdentity();
        matrix.llt().solveInPlace(matrix_inverse_);
      }
    };

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

    template<int RowsAtCompileTime, int ColsAtCompileTime = RowsAtCompileTime>
    struct MatrixInversionImpl : MatrixInversionDynamicMatrixImpl
    {
    };

    template<>
    struct MatrixInversionImpl<1> : MatrixInversionDefaultImpl
    {
    };
    template<>
    struct MatrixInversionImpl<2> : MatrixInversionDefaultImpl
    {
    };
    template<>
    struct MatrixInversionImpl<3> : MatrixInversionDefaultImpl
    {
    };
    template<>
    struct MatrixInversionImpl<4> : MatrixInversionDefaultImpl
    {
    };

    template<
      typename InputMatrix,
      bool is_floating_point = pinocchio::is_floating_point<typename InputMatrix::Scalar>::value>
    struct MatrixInversion
    : MatrixInversionImpl<InputMatrix::RowsAtCompileTime, InputMatrix::ColsAtCompileTime>
    {
    };

    template<typename InputMatrix>
    struct MatrixInversion<InputMatrix, false>
    {
      template<typename M1, typename M2>
      static EIGEN_STRONG_INLINE void
      run(const Eigen::MatrixBase<M1> & matrix, const Eigen::MatrixBase<M2> & matrix_inverse)
      {
        inverse(matrix, matrix_inverse.const_cast_derived());
      }
    };

  } // namespace internal

  template<typename M1, typename M2>
  EIGEN_STRONG_INLINE void matrix_inversion(
    const Eigen::MatrixBase<M1> & matrix, const Eigen::MatrixBase<M2> & matrix_inverse)
  {
    internal::MatrixInversion<M1>::run(matrix, matrix_inverse.const_cast_derived());
  }

  template<typename M1, typename M2>
  EIGEN_STRONG_INLINE void matrix_inversion_code_generated(
    const Eigen::MatrixBase<M1> & matrix, const Eigen::MatrixBase<M2> & matrix_inverse)
  {
    typedef internal::MatrixInversionCodeGeneratedImpl<M1::RowsAtCompileTime, M1::ColsAtCompileTime>
      Runner;
    Runner::run(matrix, matrix_inverse.const_cast_derived());
  }
} // namespace pinocchio

#include "pinocchio/math/details/matrix-inverse-1x1.hpp"
#include "pinocchio/math/details/matrix-inverse-2x2.hpp"
#include "pinocchio/math/details/matrix-inverse-3x3.hpp"
#include "pinocchio/math/details/matrix-inverse-4x4.hpp"
#include "pinocchio/math/details/matrix-inverse-5x5.hpp"
#include "pinocchio/math/details/matrix-inverse-6x6.hpp"
#include "pinocchio/math/details/matrix-inverse-7x7.hpp"

#endif // ifndef __pinocchio_math_matrix_inverse_hpp__
