//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_math_gram_schmidt_orthonormalisation_hpp__
#define __pinocchio_math_gram_schmidt_orthonormalisation_hpp__

#include "pinocchio/math/fwd.hpp"

namespace pinocchio
{
  ///  \brief Perform the Gram-Schmidt orthogonalization on the input/output vector for a given
  /// input basis
  ///
  ///  \param[in] basis Orthonormal basis.
  ///  \param[in,out] vec Vector to orthonomarlize wrt the input basis.
  ///  \param[in] threshold Only perform the orthonormalization if the scalar product between the
  /// current column and the input vector is above the given threshold.
  ///
  template<typename MatrixType, typename VectorType>
  void orthogonalization(
    const Eigen::MatrixBase<MatrixType> & basis,
    const Eigen::MatrixBase<VectorType> & vec_,
    const typename MatrixType::Scalar & threshold = 0)
  {
    typedef typename VectorType::Scalar Scalar;
    VectorType & vec = vec_.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(basis.rows(), vec.size());
    //    assert((basis.transpose() * basis).isIdentity() && "The input basis is not orthonormal.");

    for (Eigen::DenseIndex col_id = 0; col_id < basis.cols(); ++col_id)
    {
      const auto col = basis.col(col_id);
      const Scalar alpha = col.dot(vec);
      if (math::fabs(alpha) >= threshold) // only perform the orthonormalization if the scalar
                                          // product is above a certain threshold
        vec -= alpha * col;
    }

    //    std::cout << "basis.transpose() * vec: " << (basis.transpose() *
    //    vec).eval().array().abs().maxCoeff() << std::endl; if(!(basis.transpose() *
    //    vec).isZero(1e-10))
    //    {
    //      const auto res = (basis.transpose() * vec).eval();
    //      for(Eigen::Index i = 0; i < res.size(); ++i)
    //      {
    //        std::cout << "i - " << i << " value: " << res[i] << std::endl;
    //      }
    //      std::cout << "vec: " << vec.norm() << std::endl;
    //    }

    //    if (threshold == 0)
    //      assert((basis.transpose() * vec).isZero(1e-10));
  }

  ///  \brief Perform the orthonormalization of the input matrix via the Gram-Schmidt procedure.
  ///
  ///  \param[in,out] matrix_ Matrix to orthonormalize.
  ///  \param[in] threshold Only perform the orthonormalization if the scalar product between the
  /// current column and the input vector is above the given threshold.
  ///
  template<typename MatrixType>
  void orthonormalization(
    const Eigen::MatrixBase<MatrixType> & matrix_,
    const typename MatrixType::Scalar & threshold = 0)
  {
    auto & matrix = matrix_.const_cast_derived();

    for (Eigen::Index i = 0; i < matrix.cols(); ++i)
    {
      auto vec = matrix.col(i);
      if (i > 0)
        orthogonalization(matrix.leftCols(i), vec);

      const auto vec_norm = vec.norm();
      vec /= vec_norm;
    }
  }

  ///  \brief Check whether the input matrix is orthonormal up to a given input precision.
  ///
  ///  \param[in,out] matrix_ Matrix to orthonormalize.
  ///  \param[in] prec Numerical precision of the check (optional).
  ///
  template<typename MatrixType>
  bool isOrthonormal(
    const Eigen::MatrixBase<MatrixType> & matrix,
    const typename MatrixType::Scalar & prec =
      Eigen::NumTraits<typename MatrixType::Scalar>::dummy_precision())
  {
    return (matrix.transpose() * matrix).isIdentity(prec);
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_math_gram_schmidt_orthonormalisation_hpp__
