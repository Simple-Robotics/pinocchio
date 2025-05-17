//
// Copyright (c) 2019-2025 INRIA
//

#ifndef __pinocchio_math_matrix_inverse_hpp__
#define __pinocchio_math_matrix_inverse_hpp__

#include "pinocchio/math/matrix.hpp"
#include "pinocchio/math/matrix-inverse-code-generated.hpp"

namespace pinocchio
{
  namespace internal
  {
    template<int Size = ::Eigen::Dynamic>
    struct MatrixInversionEigenDefaultImpl
    {
      template<typename M1, typename M2>
      static EIGEN_STRONG_INLINE void
      run(const Eigen::MatrixBase<M1> & matrix, const Eigen::MatrixBase<M2> & matrix_inverse)
      {
        matrix_inverse.const_cast_derived().noalias() =
          matrix.template block<Size, Size>(0, 0).inverse();
      }
    };

    template<>
    struct MatrixInversionEigenDefaultImpl<::Eigen::Dynamic>
    {
      template<typename M1, typename M2>
      static EIGEN_STRONG_INLINE void
      run(const Eigen::MatrixBase<M1> & matrix, const Eigen::MatrixBase<M2> & matrix_inverse)
      {
        matrix_inverse.const_cast_derived().noalias() = matrix.inverse();
      }
    };

    template<int RowsAtCompileTime, int ColsAtCompileTime = RowsAtCompileTime>
    struct MatrixInversionImpl;

#define SET_MATRIX_INVERSION_FOR(size, Impl)                                                       \
  template<>                                                                                       \
  struct MatrixInversionImpl<size> : Impl                                                          \
  {                                                                                                \
  };

    // For size lower than 4, we can use the spezialized inverse of Eigen
    SET_MATRIX_INVERSION_FOR(1, MatrixInversionEigenDefaultImpl<1>)
    SET_MATRIX_INVERSION_FOR(2, MatrixInversionEigenDefaultImpl<2>)
    SET_MATRIX_INVERSION_FOR(3, MatrixInversionEigenDefaultImpl<3>)
    SET_MATRIX_INVERSION_FOR(4, MatrixInversionEigenDefaultImpl<4>)

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

    struct MatrixInversionDynamicMatrixImpl
    {
      template<typename M1, typename M2>
      static EIGEN_STRONG_INLINE void
      run(const Eigen::MatrixBase<M1> & matrix, const Eigen::MatrixBase<M2> & matrix_inverse)
      {
        PINOCCHIO_MAYBE_UNUSED typedef typename M1::RealScalar RealScalar;
        assert(is_symmetric(matrix, math::sqrt(dummy_precision<RealScalar>())));

#define CASE_RUN_SPECIFIC_MATRIX_INVERSE_FOR_SIZE(size)                                            \
case size:                                                                                         \
  MatrixInversionImpl<size>::run(matrix.derived(), matrix_inverse.const_cast_derived());           \
  break;

        switch (matrix.rows())
        {
          CASE_RUN_SPECIFIC_MATRIX_INVERSE_FOR_SIZE(1)
          CASE_RUN_SPECIFIC_MATRIX_INVERSE_FOR_SIZE(2)
          CASE_RUN_SPECIFIC_MATRIX_INVERSE_FOR_SIZE(3)
          CASE_RUN_SPECIFIC_MATRIX_INVERSE_FOR_SIZE(4)
          CASE_RUN_SPECIFIC_MATRIX_INVERSE_FOR_SIZE(5)
          CASE_RUN_SPECIFIC_MATRIX_INVERSE_FOR_SIZE(6)
          CASE_RUN_SPECIFIC_MATRIX_INVERSE_FOR_SIZE(7)
          CASE_RUN_SPECIFIC_MATRIX_INVERSE_FOR_SIZE(8)
          CASE_RUN_SPECIFIC_MATRIX_INVERSE_FOR_SIZE(9)
          CASE_RUN_SPECIFIC_MATRIX_INVERSE_FOR_SIZE(10)
          CASE_RUN_SPECIFIC_MATRIX_INVERSE_FOR_SIZE(11)
          CASE_RUN_SPECIFIC_MATRIX_INVERSE_FOR_SIZE(12)
        default:
          generic(matrix.derived(), matrix_inverse.const_cast_derived());
          break;
        }

#undef CASE_RUN_SPECIFIC_MATRIX_INVERSE_FOR_SIZE
      }

      template<typename M1, typename M2>
      static EIGEN_STRONG_INLINE void
      generic(const Eigen::MatrixBase<M1> & matrix, const Eigen::MatrixBase<M2> & matrix_inverse)
      {
        auto & matrix_inverse_ = matrix_inverse.const_cast_derived();

        matrix_inverse_.setIdentity();
#ifdef PINOCCHIO_MAC_ARM64
        matrix.ldlt().solveInPlace(matrix_inverse_);
#else
        matrix.llt().solveInPlace(matrix_inverse_);
#endif
      }
    };

    template<int RowsAtCompileTime, int ColsAtCompileTime>
    struct MatrixInversionImpl : MatrixInversionDynamicMatrixImpl
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
} // namespace pinocchio

#endif // ifndef __pinocchio_math_matrix_inverse_hpp__
