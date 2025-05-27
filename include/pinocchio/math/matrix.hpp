//
// Copyright (c) 2018-2025 INRIA
// Copyright (c) 2016-2018 CNRS
//

#ifndef __pinocchio_math_matrix_hpp__
#define __pinocchio_math_matrix_hpp__

#include "pinocchio/macros.hpp"
#include "pinocchio/math/fwd.hpp"
#include "pinocchio/utils/static-if.hpp"

#include <boost/type_traits.hpp>

namespace pinocchio
{

  template<typename Derived>
  inline bool hasNaN(const Eigen::DenseBase<Derived> & m)
  {
    return !((m.derived().array() == m.derived().array()).all());
  }

  template<typename Matrix1, typename Matrix2>
  bool isApproxOrZero(
    const Eigen::MatrixBase<Matrix1> & mat1,
    const Eigen::MatrixBase<Matrix2> & mat2,
    const typename Matrix1::RealScalar & prec =
      Eigen::NumTraits<typename Matrix1::RealScalar>::dummy_precision())
  {
    const bool mat1_is_zero = mat1.isZero(prec);
    const bool mat2_is_zero = mat2.isZero(prec);

    const bool mat1_is_approx_mat2 = mat1.isApprox(mat2, prec);
    return mat1_is_approx_mat2 || (mat1_is_zero && mat2_is_zero);
  }

  namespace internal
  {
    template<
      typename MatrixLike,
      bool value = is_floating_point<typename MatrixLike::Scalar>::value>
    struct isZeroAlgo
    {
      typedef typename MatrixLike::Scalar Scalar;
      typedef typename MatrixLike::RealScalar RealScalar;

      static bool run(
        const Eigen::MatrixBase<MatrixLike> & mat,
        const RealScalar & prec = Eigen::NumTraits<Scalar>::dummy_precision())
      {
        return mat.isZero(prec);
      }
    };

    template<typename MatrixLike>
    struct isZeroAlgo<MatrixLike, false>
    {
      typedef typename MatrixLike::Scalar Scalar;
      typedef typename MatrixLike::RealScalar RealScalar;

      static bool run(
        const Eigen::MatrixBase<MatrixLike> & /*vec*/,
        const RealScalar & prec = Eigen::NumTraits<Scalar>::dummy_precision())
      {
        PINOCCHIO_UNUSED_VARIABLE(prec);
        return true;
      }
    };
  } // namespace internal

  template<typename MatrixLike>
  inline bool isZero(
    const Eigen::MatrixBase<MatrixLike> & m,
    const typename MatrixLike::RealScalar & prec =
      Eigen::NumTraits<typename MatrixLike::Scalar>::dummy_precision())
  {
    return internal::isZeroAlgo<MatrixLike>::run(m, prec);
  }

  template<typename M1, typename M2>
  struct MatrixMatrixProduct
  {
#if EIGEN_VERSION_AT_LEAST(3, 2, 90)
    typedef typename Eigen::Product<M1, M2> type;
#else
    typedef typename Eigen::ProductReturnType<M1, M2>::Type type;
#endif
  };

  template<typename Scalar, typename Matrix>
  struct ScalarMatrixProduct
  {
#if EIGEN_VERSION_AT_LEAST(3, 3, 0)
    typedef Eigen::CwiseBinaryOp<
      EIGEN_CAT(EIGEN_CAT(Eigen::internal::scalar_, product), _op) < Scalar,
      typename Eigen::internal::traits<Matrix>::Scalar>,
      const typename Eigen::internal::plain_constant_type<Matrix, Scalar>::type,
      const Matrix > type;
#elif EIGEN_VERSION_AT_LEAST(3, 2, 90)
    typedef Eigen::CwiseUnaryOp<Eigen::internal::scalar_multiple_op<Scalar>, const Matrix> type;
#else
    typedef const Eigen::CwiseUnaryOp<Eigen::internal::scalar_multiple_op<Scalar>, const Matrix>
      type;
#endif
  };

  template<typename Matrix, typename Scalar>
  struct MatrixScalarProduct
  {
#if EIGEN_VERSION_AT_LEAST(3, 3, 0)
    typedef Eigen::CwiseBinaryOp<
      EIGEN_CAT(EIGEN_CAT(Eigen::internal::scalar_, product), _op) <
        typename Eigen::internal::traits<Matrix>::Scalar,
      Scalar>,
      const Matrix,
      const typename Eigen::internal::plain_constant_type<Matrix, Scalar>::type > type;
#elif EIGEN_VERSION_AT_LEAST(3, 2, 90)
    typedef Eigen::CwiseUnaryOp<Eigen::internal::scalar_multiple_op<Scalar>, const Matrix> type;
#else
    typedef const Eigen::CwiseUnaryOp<Eigen::internal::scalar_multiple_op<Scalar>, const Matrix>
      type;
#endif
  };

  namespace internal
  {
    template<
      typename MatrixLike,
      bool value = is_floating_point<typename MatrixLike::Scalar>::value>
    struct isUnitaryAlgo
    {
      typedef typename MatrixLike::Scalar Scalar;
      typedef typename MatrixLike::RealScalar RealScalar;

      static bool run(
        const Eigen::MatrixBase<MatrixLike> & mat,
        const RealScalar & prec = Eigen::NumTraits<Scalar>::dummy_precision())
      {
        return mat.isUnitary(prec);
      }
    };

    template<typename MatrixLike>
    struct isUnitaryAlgo<MatrixLike, false>
    {
      typedef typename MatrixLike::Scalar Scalar;
      typedef typename MatrixLike::RealScalar RealScalar;

      static bool run(
        const Eigen::MatrixBase<MatrixLike> & /*vec*/,
        const RealScalar & prec = Eigen::NumTraits<Scalar>::dummy_precision())
      {
        PINOCCHIO_UNUSED_VARIABLE(prec);
        return true;
      }
    };
  } // namespace internal

  ///
  /// \brief Check whether the input matrix is Unitary within the given precision.
  ///
  /// \param[in] mat Input matrix
  /// \param[in] prec Required precision
  ///
  /// \returns true if mat is unitary within the precision prec
  ///
  template<typename MatrixLike>
  inline bool isUnitary(
    const Eigen::MatrixBase<MatrixLike> & mat,
    const typename MatrixLike::RealScalar & prec =
      Eigen::NumTraits<typename MatrixLike::Scalar>::dummy_precision())
  {
    return internal::isUnitaryAlgo<MatrixLike>::run(mat, prec);
  }

  namespace internal
  {
    template<
      typename VectorLike,
      bool value = is_floating_point<typename VectorLike::Scalar>::value>
    struct isNormalizedAlgo
    {
      typedef typename VectorLike::Scalar Scalar;
      typedef typename VectorLike::RealScalar RealScalar;

      static bool run(
        const Eigen::MatrixBase<VectorLike> & vec,
        const RealScalar & prec = Eigen::NumTraits<RealScalar>::dummy_precision())
      {
        return math::fabs(static_cast<RealScalar>(vec.norm() - RealScalar(1))) <= prec;
      }
    };

    template<typename VectorLike>
    struct isNormalizedAlgo<VectorLike, false>
    {
      typedef typename VectorLike::Scalar Scalar;
      typedef typename VectorLike::RealScalar RealScalar;

      static bool run(
        const Eigen::MatrixBase<VectorLike> & /*vec*/,
        const RealScalar & prec = Eigen::NumTraits<RealScalar>::dummy_precision())
      {
        PINOCCHIO_UNUSED_VARIABLE(prec);
        return true;
      }
    };
  } // namespace internal

  ///
  /// \brief Check whether the input vector is Normalized within the given precision.
  ///
  /// \param[in] vec Input vector
  /// \param[in] prec Required precision
  ///
  /// \returns true if vec is normalized within the precision prec.
  ///
  template<typename VectorLike>
  inline bool isNormalized(
    const Eigen::MatrixBase<VectorLike> & vec,
    const typename VectorLike::RealScalar & prec =
      Eigen::NumTraits<typename VectorLike::Scalar>::dummy_precision())
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike);
    return internal::isNormalizedAlgo<VectorLike>::run(vec, prec);
  }

  namespace internal
  {
    template<
      typename VectorLike,
      bool value = is_floating_point<typename VectorLike::Scalar>::value>
    struct normalizeAlgo
    {
      static void run(const Eigen::MatrixBase<VectorLike> & vec)
      {
        return vec.const_cast_derived().normalize();
      }
    };

    template<typename VectorLike>
    struct normalizeAlgo<VectorLike, false>
    {
      static void run(const Eigen::MatrixBase<VectorLike> & vec)
      {
        using namespace internal;
        typedef typename VectorLike::RealScalar RealScalar;
        typedef typename VectorLike::Scalar Scalar;
        const RealScalar z = vec.squaredNorm();
        const Scalar sqrt_z = if_then_else(GT, z, Scalar(0), math::sqrt(z), Scalar(1));
        vec.const_cast_derived() /= sqrt_z;
      }
    };
  } // namespace internal

  ///
  /// \brief Normalize the input vector.
  ///
  /// \param[in] vec Input vector
  ///
  template<typename VectorLike>
  inline void normalize(const Eigen::MatrixBase<VectorLike> & vec)
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike);
    internal::normalizeAlgo<VectorLike>::run(vec.const_cast_derived());
  }

  namespace internal
  {
    template<typename Scalar>
    struct CallCorrectMatrixInverseAccordingToScalar
    {
      template<typename MatrixIn, typename MatrixOut>
      static void
      run(const Eigen::MatrixBase<MatrixIn> & m_in, const Eigen::MatrixBase<MatrixOut> & dest)
      {
        MatrixOut & dest_ = PINOCCHIO_EIGEN_CONST_CAST(MatrixOut, dest);
        dest_.noalias() = m_in.inverse();
      }
    };

  } // namespace internal

  template<typename MatrixIn, typename MatrixOut>
  inline void
  inverse(const Eigen::MatrixBase<MatrixIn> & m_in, const Eigen::MatrixBase<MatrixOut> & dest)
  {
    MatrixOut & dest_ = PINOCCHIO_EIGEN_CONST_CAST(MatrixOut, dest);
    internal::CallCorrectMatrixInverseAccordingToScalar<typename MatrixIn::Scalar>::run(
      m_in, dest_);
  }

  namespace internal
  {
    template<
      typename MatrixLike,
      bool value = is_floating_point<typename MatrixLike::Scalar>::value>
    struct isSymmetricAlgo
    {
      typedef typename MatrixLike::Scalar Scalar;
      typedef typename MatrixLike::RealScalar RealScalar;

      static bool run(
        const Eigen::MatrixBase<MatrixLike> & mat,
        const RealScalar & prec = Eigen::NumTraits<RealScalar>::dummy_precision())
      {
        if (mat.rows() != mat.cols())
          return false;
        return (mat - mat.transpose()).isZero(prec);
      }
    };

    template<typename MatrixLike>
    struct isSymmetricAlgo<MatrixLike, false>
    {
      typedef typename MatrixLike::Scalar Scalar;
      typedef typename MatrixLike::RealScalar RealScalar;

      static bool run(
        const Eigen::MatrixBase<MatrixLike> & /*mat*/,
        const RealScalar & prec = Eigen::NumTraits<RealScalar>::dummy_precision())
      {
        PINOCCHIO_UNUSED_VARIABLE(prec);
        return true;
      }
    };
  } // namespace internal

  ///
  /// \brief Check whether the input matrix is symmetric within the given precision.
  ///
  /// \param[in] mat Input matrix
  /// \param[in] prec Required precision
  ///
  /// \returns true if mat is symmetric within the precision prec.
  ///
  template<typename MatrixLike>
  inline bool isSymmetric(
    const Eigen::MatrixBase<MatrixLike> & mat,
    const typename MatrixLike::RealScalar & prec =
      Eigen::NumTraits<typename MatrixLike::Scalar>::dummy_precision())
  {
    return internal::isSymmetricAlgo<MatrixLike>::run(mat, prec);
  }

  namespace internal
  {
    template<typename XprType, typename DestType, typename Weak = void>
    struct evalToImpl
    {
      static void run(const XprType & xpr, DestType & dest)
      {
        xpr.evalTo(dest);
      }
    };

    template<typename X1, typename X2, typename DenseDerived>
    struct evalToImpl<Eigen::Product<X1, X2>, DenseDerived, void>
    {
      typedef Eigen::MatrixBase<DenseDerived> DestType;
      typedef Eigen::Product<X1, X2> XprType;
      static void run(const XprType & xpr, DestType & dest)
      {
        dest.noalias() = xpr;
      }
    };

  } // namespace internal

  namespace internal
  {
    template<
      typename M1,
      typename M2,
      bool value = is_floating_point<typename M1::Scalar>::value
                   && is_floating_point<typename M2::Scalar>::value>
    struct arrayCompareAllAlgo
    {
      static bool run(
        const Eigen::MatrixBase<M1> & m1,
        const Eigen::MatrixBase<M2> & m2,
        internal::ComparisonOperators op)
      {
        switch (op)
        {
        case LT:
          return (m1.array() < m2.array()).all();
        case LE:
          return (m1.array() <= m2.array()).all();
        case EQ:
          return (m1.array() == m2.array()).all();
        case GE:
          return (m1.array() >= m2.array()).all();
        case GT:
          return (m1.array() > m2.array()).all();
        default:
          return false;
        }
      }
    };

    template<typename M1, typename M2>
    struct arrayCompareAllAlgo<M1, M2, false>
    {

      static bool run(
        const Eigen::MatrixBase<M1> & m1,
        const Eigen::MatrixBase<M2> & m2,
        internal::ComparisonOperators op)
      {
        PINOCCHIO_UNUSED_VARIABLE(m1);
        PINOCCHIO_UNUSED_VARIABLE(m2);
        PINOCCHIO_UNUSED_VARIABLE(op);
        return true;
      }
    };
  } // namespace internal

  ///
  /// \brief Compare all the elements of two matrices
  ///
  /// \param[in] m1 Input matrix
  /// \param[in] m2 Input matrix
  /// \param[in] op Comparison operation
  ///
  /// \returns true if op is true on all m1 and m2 elements
  ///
  template<typename M1, typename M2>
  inline bool arrayCompareAll(
    const Eigen::MatrixBase<M1> & m1,
    const Eigen::MatrixBase<M2> & m2,
    internal::ComparisonOperators op)
  {
    return internal::arrayCompareAllAlgo<M1, M2>::run(m1, m2, op);
  }

  namespace internal
  {
    template<
      typename M,
      typename Lower,
      typename Upper,
      typename Res,
      bool value = is_floating_point<typename M::Scalar>::value
                   && is_floating_point<typename Lower::Scalar>::value
                   && is_floating_point<typename Upper::Scalar>::value
                   && is_floating_point<typename Res::Scalar>::value>
    struct arrayBoundAlgo
    {
      static void run(
        const Eigen::MatrixBase<M> & x,
        const Eigen::MatrixBase<Lower> & lower,
        const Eigen::MatrixBase<Upper> & upper,
        const Eigen::MatrixBase<Res> & res)
      {
        res.const_cast_derived() = x.array().max(lower.array()).min(upper.array());
      }
    };

    template<typename M, typename Lower, typename Upper, typename Res>
    struct arrayBoundAlgo<M, Lower, Upper, Res, false>
    {

      static void run(
        const Eigen::MatrixBase<M> & x,
        const Eigen::MatrixBase<Lower> & lower,
        const Eigen::MatrixBase<Upper> & upper,
        const Eigen::MatrixBase<Res> & res)
      {
        PINOCCHIO_UNUSED_VARIABLE(lower);
        PINOCCHIO_UNUSED_VARIABLE(upper);
        res.const_cast_derived() = x;
      }
    };
  } // namespace internal

  ///
  /// \brief Bound a matrix between upper and lower bound
  ///
  /// \param[in] x Input matrix
  /// \param[in] lower lower bound matrix
  /// \param[in] upper upper bound matrix
  /// \param[out] res result
  ///
  template<typename M, typename Lower, typename Upper, typename Res>
  inline void arrayBound(
    const Eigen::MatrixBase<M> & x,
    const Eigen::MatrixBase<Lower> & lower,
    const Eigen::MatrixBase<Upper> & upper,
    const Eigen::MatrixBase<Res> & res)
  {
    internal::arrayBoundAlgo<M, Lower, Upper, Res>::run(x, lower, upper, res);
  }

  namespace internal
  {
    template<
      typename M,
      typename Res,
      bool value = is_floating_point<typename M::Scalar>::value
                   && is_floating_point<typename Res::Scalar>::value>
    struct arrayMinAlgo
    {
      static void run(
        const Eigen::MatrixBase<M> & x,
        typename Eigen::MatrixBase<M>::Scalar min,
        const Eigen::MatrixBase<Res> & res)
      {
        res.const_cast_derived() = x.array().min(min);
      }
    };

    template<typename M, typename Res>
    struct arrayMinAlgo<M, Res, false>
    {

      static void run(
        const Eigen::MatrixBase<M> & x,
        typename Eigen::MatrixBase<M>::Scalar min,
        const Eigen::MatrixBase<Res> & res)
      {
        PINOCCHIO_UNUSED_VARIABLE(min);
        res.const_cast_derived() = x;
      }
    };
  } // namespace internal

  ///
  /// \brief Apply an upper bound to a matrix
  ///
  /// \param[in] x Input matrix
  /// \param[in] min Matrix upper bound
  /// \param[out] res result
  ///
  template<typename M, typename Res>
  inline void arrayMin(
    const Eigen::MatrixBase<M> & x,
    typename Eigen::MatrixBase<M>::Scalar min,
    const Eigen::MatrixBase<Res> & res)
  {
    internal::arrayMinAlgo<M, Res>::run(x, min, res);
  }

  template<typename XprType, typename DestType>
  inline void evalTo(const XprType & xpr, DestType & dest)
  {
    internal::evalToImpl<XprType, DestType>::run(xpr, dest);
  }

  template<typename Matrix>
  Eigen::Ref<Matrix> make_ref(const Eigen::PlainObjectBase<Matrix> & mat)
  {
    typedef Eigen::Ref<Matrix> ReturnType;
    return ReturnType(mat.const_cast_derived());
  }

  /// \brief Helper to make the matrix symmetric
  ///
  /// \param[in,out] mat Input matrix to symmetrize.
  /// \param[in] mode Part of the matrix to symmetrize : Eigen::Upper or Eigen::Lower
  template<typename Matrix>
  void make_symmetric(const Eigen::MatrixBase<Matrix> & mat, const int mode = Eigen::Upper)
  {
    PINOCCHIO_CHECK_INPUT_ARGUMENT(mode == Eigen::Upper || mode == Eigen::Lower);

    if (mode == Eigen::Upper)
    {
      mat.const_cast_derived().template triangularView<Eigen::StrictlyLower>() =
        mat.transpose().template triangularView<Eigen::StrictlyLower>();
    }
    else if (mode == Eigen::Lower)
    {
      mat.const_cast_derived().template triangularView<Eigen::StrictlyUpper>() =
        mat.transpose().template triangularView<Eigen::StrictlyUpper>();
    }
  }

  ///
  /// \brief Helper to check whether the input matrix is square.
  ///
  /// \param[in] mat Input matrix to check whether it is square.
  ///
  template<typename Matrix>
  EIGEN_STRONG_INLINE bool is_square(const Eigen::MatrixBase<Matrix> & mat)
  {
    return mat.rows() == mat.cols();
  }

  ///
  /// \brief Helper to check whether the input matrix is symmetric.
  ///
  /// \param[in] mat Input matrix to check symmetry.
  /// \param[in] prec Numerical precision of the check (optional).
  ///
  template<typename Matrix>
  bool is_symmetric(
    const Eigen::MatrixBase<Matrix> & mat,
    const typename Matrix::Scalar & prec =
      Eigen::NumTraits<typename Matrix::Scalar>::dummy_precision())
  {
    if (!is_square(mat))
      return false;

    return mat.reshaped().isApprox(mat.transpose().reshaped(), prec);
  }

  /// \brief Enforce the symmetry of the input matrix
  template<typename Matrix>
  void enforceSymmetry(const Eigen::MatrixBase<Matrix> & mat_)
  {
    PINOCCHIO_CHECK_SQUARE_MATRIX(mat_);

    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix) PlainMatrix;
    typedef typename Matrix::Scalar Scalar;
    typedef Eigen::Map<PlainMatrix> MapMatrix;

    auto & mat = mat_.const_cast_derived();
    MapMatrix tmp = MapMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, mat.rows(), mat.rows()));

    tmp = 0.5 * (mat + mat.transpose());
    mat = tmp;

    assert(mat == mat.transpose());
  }

  template<typename Matrix>
  typename PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix) make_copy(const Eigen::MatrixBase<Matrix> & mat)
  {
    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix) ReturnType;
    return ReturnType(mat);
  }

  template<typename T>
  struct is_eigen_ref : std::false_type
  {
  };

  template<typename PlainObjectType, int Options, typename StrideType>
  struct is_eigen_ref<Eigen::Ref<PlainObjectType, Options, StrideType>> : std::true_type
  {
  };
} // namespace pinocchio

#endif // #ifndef __pinocchio_math_matrix_hpp__
