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
    struct MatrixInversionEigenDefaultImpl
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
        typedef typename M1::RealScalar RealScalar;
        assert(is_symmetric(matrix, math::sqrt(dummy_precision<RealScalar>())));

        auto & matrix_inverse_ = matrix_inverse.const_cast_derived();

        matrix_inverse_.setIdentity();
#ifdef PINOCCHIO_MAC_ARM64
        matrix.ldlt().solveInPlace(matrix_inverse_);
#else
        matrix.llt().solveInPlace(matrix_inverse_);
#endif
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

namespace pinocchio
{
  namespace internal
  {

    template<int RowsAtCompileTime, int ColsAtCompileTime = RowsAtCompileTime>
    struct MatrixInversionImpl : MatrixInversionDynamicMatrixImpl
    {
    };

#define SET_MATRIX_INVERSION_FOR(size, Impl)                                                       \
  template<>                                                                                       \
  struct MatrixInversionImpl<size> : Impl                                                          \
  {                                                                                                \
  };

    // For size lower than 4, we can use the spezialized inverse of Eigen
    SET_MATRIX_INVERSION_FOR(1, MatrixInversionEigenDefaultImpl)
    SET_MATRIX_INVERSION_FOR(2, MatrixInversionEigenDefaultImpl)
    SET_MATRIX_INVERSION_FOR(3, MatrixInversionEigenDefaultImpl)
    SET_MATRIX_INVERSION_FOR(4, MatrixInversionEigenDefaultImpl)

    // For size in [5,12], we can use code generated impl
    SET_MATRIX_INVERSION_FOR(5, MatrixInversionCodeGeneratedImpl<5>)
    SET_MATRIX_INVERSION_FOR(6, MatrixInversionCodeGeneratedImpl<6>)
    SET_MATRIX_INVERSION_FOR(7, MatrixInversionCodeGeneratedImpl<7>)
    SET_MATRIX_INVERSION_FOR(8, MatrixInversionCodeGeneratedImpl<8>)
    SET_MATRIX_INVERSION_FOR(9, MatrixInversionCodeGeneratedImpl<9>)
    SET_MATRIX_INVERSION_FOR(10, MatrixInversionCodeGeneratedImpl<10>)
    SET_MATRIX_INVERSION_FOR(11, MatrixInversionCodeGeneratedImpl<11>)
    SET_MATRIX_INVERSION_FOR(12, MatrixInversionCodeGeneratedImpl<12>)

#undef SET_MATRIX_INVERSION_FOR

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

#endif // ifndef __pinocchio_math_matrix_inverse_hpp__
