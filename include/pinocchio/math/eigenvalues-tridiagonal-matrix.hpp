//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_math_eigenvalues_tridiagonal_matrix_hpp__
#define __pinocchio_math_eigenvalues_tridiagonal_matrix_hpp__

#include "pinocchio/math/fwd.hpp"

namespace pinocchio
{
  template<typename Scalar, int Options> struct TridiagonalSymmetricMatrixTpl;
  
  ///
  /// \brief Computes the spectrum[m1:m2] of the input tridiagonal matrix up to precision eps
  ///
  /// \param[in] tridiagonal_mat a Tridiagonal Symmetric matrix
  /// \param[in] m1 the index of the first eigenvalue to compute (lowest)
  /// \param[in] m2 the index of the last eigenvalue to compute (largest)
  /// \param[in] eps tolerance in the estimate of the eigenvalues
  ///
  /// \details This functions implements the seminal work of W. BARTH, R. S. MARTIN and J. H. WILKINSON which can be downloaded at https://link.springer.com/content/pdf/10.1007/BF02162154.pdf
  /// \remarks This function proceeds to some minimal memory allocation for efficiency
  ///
  template<typename Scalar, int Options>
  Eigen::Matrix<Scalar,Eigen::Dynamic,1,Options>
  computeSpectrum(const TridiagonalSymmetricMatrixTpl<Scalar,Options> & tridiagonal_mat,
                  const Eigen::DenseIndex m1,
                  const Eigen::DenseIndex m2,
                  Scalar eps = 1e-4)
  {
    typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1,Options> ReturnType;
    typedef TridiagonalSymmetricMatrixTpl<Scalar,Options> TridiagonalSymmetricMatrix;
    typedef typename TridiagonalSymmetricMatrix::CoeffVectorType CoeffVectorType;
    
    PINOCCHIO_CHECK_INPUT_ARGUMENT(m1 <= m2, "m1 should be lower than m2.");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(m1 >= 0, "m1 should be greater than 0.");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(m2 >= 0, "m2 should be greater than 0.");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(m2 <= tridiagonal_mat.rows(), "m2 should be lower than the size of the tridiagonal matrix.");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(eps >= Scalar(0), "eps should be greater than 0.")
    
    const Eigen::DenseIndex n = tridiagonal_mat.rows();
    const Eigen::DenseIndex dm = m2 - m1 + 1;
    const Scalar relfeh = 2 * Eigen::NumTraits<Scalar>::epsilon();
    
    assert((Scalar(1) + relfeh) > Scalar(1));
    
    const auto & alphas = tridiagonal_mat.diagonal();
    const auto & betas_ = tridiagonal_mat.subDiagonal();
    CoeffVectorType betas_abs = CoeffVectorType::Zero(n);
    betas_abs.array().tail(n-1) = betas_.array().abs();
    CoeffVectorType betas_square = betas_abs.array().square();
    
    Scalar
    xmin = alphas[n-1] - betas_abs[n-1],
    xmax = alphas[n-1] + betas_abs[n-1];
    
    for(Eigen::DenseIndex i = n-2; i >= 0; --i)
    {
      const Scalar h = betas_abs[i] + betas_abs[i+1];
      xmax = math::max(alphas[i] + h,xmax);
      xmin = math::min(alphas[i] - h,xmin);
    }
    
    
    Scalar eps2 = relfeh * ((xmin + xmax > 0) ? xmax : xmin);
    eps2 = 0.5 * eps + 7 * eps2;
    
    // Inner block
    Scalar x0 = xmax;
    ReturnType spectrum = ReturnType::Zero(n);
    auto & x = spectrum;
    CoeffVectorType wu = CoeffVectorType::Zero(n);
    
    x.segment(m1,dm).fill(xmax);
    wu.segment(m1,dm).fill(xmin);
    
//    Eigen::DenseIndex z = 0;
    // Loop for the kth eigenvalue
    
    for(Eigen::DenseIndex k = m2; k >= m1; --k)
    {
      Scalar xu = xmin;
      for(Eigen::DenseIndex i = k; i >= m1; --i)
      {
        if (xu <= wu[i])
        {
          xu = wu[i];
          x0 = math::min(x0,x[k]);
          while ((x0 - xu) > (2 * relfeh * (math::fabs(xu) + math::fabs(x0)) + eps))
          {
//            z++;
            Scalar x1 = 0.5 * (xu + x0);
            Eigen::DenseIndex a = -1;
            Scalar q = 1.;
            for (Eigen::DenseIndex j = 0; j < n; ++j)
            {
              const Scalar dq = q != Scalar(0) ? betas_square[j] / q : betas_abs[j] / relfeh;
              q = alphas[j] - x1 - dq;
              if (q < Scalar(0)) a++;
            } // for
            if (a < k)
            {
              if (a < m1) 
              {
                xu = wu[m1] = x1;
              }
              else
              {
                xu = wu[a+1] = x1;
                x[a] = math::min(x[a],x1);
              }
            }
            else
            {
              x0 = x1;
            }
          } // end while
          x[k] = 0.5 * (xu + x0);
        }
      }
    }

    return spectrum;
  }
  
  ///
  /// \brief Computes the full spectrum of the input tridiagonal matrix up to precision eps
  ///
  /// \param[in] tridiagonal_mat a Tridiagonal Symmetric matrix
  /// \param[in] eps tolerance in the estimate of the eigenvalues
  ///
  /// \details This functions implements the seminal work of W. BARTH, R. S. MARTIN and J. H. WILKINSON which can be downloaded at https://link.springer.com/content/pdf/10.1007/BF02162154.pdf
  /// \remarks This function proceeds to some minimal memory allocation for efficiency
  ///
  template<typename Scalar, int Options>
  Eigen::Matrix<Scalar,Eigen::Dynamic,1,Options>
  computeSpectrum(const TridiagonalSymmetricMatrixTpl<Scalar,Options> & tridiagonal_mat,
                  Scalar eps = 1e-4)
  {
    return computeSpectrum(tridiagonal_mat,0,tridiagonal_mat.cols()-1,eps);
  }
  
  ///
  /// \brief Computes the kth eigenvalue associated with the input tridiagonal matrix up to precision eps
  ///
  /// \param[in] tridiagonal_mat a Tridiagonal Symmetric matrix
  /// \param[in] eigenvalue_index index of the eigenvalue to compute
  /// \param[in] eps tolerance in the estimate of the eigenvalues
  ///
  /// \returns The kth eigenvalue
  /// \see computeSpectrum
  template<typename Scalar, int Options>
  Scalar
  computeEigenvalue(const TridiagonalSymmetricMatrixTpl<Scalar,Options> & tridiagonal_mat,
                    const Eigen::DenseIndex eigenvalue_index,
                    Scalar eps = 1e-4)
  {
    return computeSpectrum(tridiagonal_mat,eigenvalue_index,eigenvalue_index,eps)[eigenvalue_index];
  }
}

#endif //#ifndef __pinocchio_math_eigenvalues_tridiagonal_matrix_hpp__
