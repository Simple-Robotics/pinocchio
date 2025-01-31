//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_math_lanczos_decomposition_hpp__
#define __pinocchio_math_lanczos_decomposition_hpp__

#include "pinocchio/math/fwd.hpp"
#include "pinocchio/math/tridiagonal-matrix.hpp"
#include "pinocchio/math/gram-schmidt-orthonormalisation.hpp"

namespace pinocchio
{

  /// \brief Compute the largest eigenvalues and the associated principle eigenvector via power
  /// iteration
  template<typename _Matrix>
  struct LanczosDecompositionTpl
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(_Matrix) PlainMatrix;
    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(typename PlainMatrix::ColXpr) Vector;
    typedef typename _Matrix::Scalar Scalar;
    enum
    {
      Options = _Matrix::Options
    };
    typedef TridiagonalSymmetricMatrixTpl<Scalar, Options> TridiagonalMatrix;

    /// \brief Default constructor for the Lanczos decomposition from an input matrix
    template<typename MatrixLikeType>
    LanczosDecompositionTpl(const MatrixLikeType & mat, const Eigen::DenseIndex decomposition_size)
    : m_Qs(mat.rows(), decomposition_size)
    , m_Ts(decomposition_size)
    , m_A_times_q(mat.rows())
    , m_rank(-1)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(mat.rows() == mat.cols(), "The input matrix is not square.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        decomposition_size >= 1, "The size of the decomposition should be greater than one.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        decomposition_size <= mat.rows(),
        "The size of the decomposition should not be larger than the number of rows.");

      compute(mat);
    }

    /// \brief Constructor for the Lanczos decomposition from given sizes.
    /// Compute must be called afterwards.
    LanczosDecompositionTpl(
      const Eigen::DenseIndex size, const Eigen::DenseIndex decomposition_size)
    : m_Qs(size, decomposition_size)
    , m_Ts(decomposition_size)
    , m_A_times_q(size)
    , m_rank(-1)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        decomposition_size >= 1, "The size of the decomposition should be greater than one.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        decomposition_size <= size,
        "The size of the decomposition should not be larger than the number of rows.");
    }

    bool operator==(const LanczosDecompositionTpl & other) const
    {
      if (this == &other)
        return true;
      return m_Qs == other.m_Qs && m_Ts == other.m_Ts && m_rank == other.m_rank;
    }

    bool operator!=(const LanczosDecompositionTpl & other) const
    {
      return !(*this == other);
    }

    ///
    /// \brief Computes the Lanczos decomposition of the input matrix A
    ///
    /// \param[in] A The matrix to decompose
    ///
    template<typename MatrixLikeType>
    void compute(const MatrixLikeType & A)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(A.rows() == A.cols(), "The input matrix is not square.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        A.rows() == m_Qs.rows(), "The input matrix is not of correct size.");

      const Eigen::DenseIndex num_cols = A.cols();
      const Eigen::DenseIndex decomposition_size = m_Ts.cols();
      auto & alphas = m_Ts.diagonal();
      auto & betas = m_Ts.subDiagonal();

      const Scalar prec = 2 * Eigen::NumTraits<Scalar>::epsilon();

      m_Qs.setIdentity();

      m_Ts.setZero();
      m_rank = 1;
      for (Eigen::DenseIndex k = 0; k < decomposition_size; ++k)
      {
        const auto q = m_Qs.col(k);
        m_A_times_q.noalias() = A * q;
        alphas[k] = q.dot(m_A_times_q);

        if (k < decomposition_size - 1)
        {
          auto q_next = m_Qs.col(k + 1);
          m_A_times_q -= alphas[k] * q;
          if (k > 0)
          {
            const auto q_previous = m_Qs.col(k - 1);
            m_A_times_q -= betas[k - 1] * q_previous;
          }

          // Perform Gram-Schmidt orthogonalization procedure.
          if (k > 0)
            orthonormalisation(m_Qs.leftCols(k + 1), m_A_times_q);
          assert(m_Qs.leftCols(k + 1).cols() == k + 1);

          // Compute beta
          betas[k] = m_A_times_q.norm();
          if (betas[k] <= prec)
          {
            // Pick a new arbitrary vector
            bool found_new_base_vector = false;

            Scalar q_next_norm = -1; //= q_next.norm();

            for (Eigen::DenseIndex j = 0; j < num_cols; ++j)
            {
              const Eigen::DenseIndex base_col_id = (k + 1 + j) % num_cols;
              if (base_col_id == 0)
                continue; // The first column of Qs is the first unit vector.

              typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(decltype(q_next)) VectorPlain;
              q_next = VectorPlain::Unit(num_cols, base_col_id);
              orthonormalisation(m_Qs.leftCols(k + 1), q_next);
              q_next_norm = q_next.norm();
              if (q_next_norm > prec)
              {
                found_new_base_vector = true;
                break;
              }
            }

            assert(found_new_base_vector && "Issue with picking a new arbitrary vector.");

            q_next /= q_next_norm;
            betas[k] = 0.;
            m_rank++;
          }
          else
          {
            q_next.noalias() = m_A_times_q / betas[k];
            m_rank++;
          }
        }
      }

      m_Qs.rightCols(decomposition_size - m_rank).setZero();
    }

    ///
    /// \brief Computes the residual associated with the decomposition, namely, the quantity \f$ A
    /// Q_s - Q_s T_s \f$
    ///
    /// \param[in] A the matrix that have been decomposed.
    ///
    /// \returns The residual of the decomposition
    ///
    template<typename MatrixLikeType>
    PlainMatrix computeDecompositionResidual(const MatrixLikeType & A) const
    {
      const Eigen::DenseIndex last_col_id = m_Ts.cols() - 1;
      const auto & alphas = m_Ts.diagonal();
      const auto & betas = m_Ts.subDiagonal();

      PlainMatrix residual = A * m_Qs;
      residual -= (m_Qs * m_Ts).eval();

      const auto q = m_Qs.col(last_col_id);

      auto & tmp_vec = m_A_times_q; // use m_A_times_q as a temporary vector
      tmp_vec.noalias() = A * q;
      tmp_vec -= alphas[last_col_id] * q;
      if (last_col_id > 0)
        tmp_vec -= betas[last_col_id - 1] * m_Qs.col(last_col_id - 1);

      residual.col(last_col_id) -= tmp_vec;

      return residual;
    }

    /// \brief Returns the tridiagonal matrix associated with the Lanczos decomposition
    const TridiagonalMatrix & Ts() const
    {
      return m_Ts;
    }
    /// \brief Returns the tridiagonal matrix associated with the Lanczos decomposition
    TridiagonalMatrix & Ts()
    {
      return m_Ts;
    }

    /// \brief Returns the orthogonal basis associated with the Lanczos decomposition
    const PlainMatrix & Qs() const
    {
      return m_Qs;
    }
    /// \brief Returns the orthogonal basis associated with the Lanczos decomposition
    PlainMatrix & Qs()
    {
      return m_Qs;
    }

    /// \brief Returns the size
    Eigen::DenseIndex size() const
    {
      return m_Qs.rows();
    }

    /// \brief Returns the size of the decomposition
    Eigen::DenseIndex decompositionSize() const
    {
      return m_Ts.rows();
    }

  protected:
    PlainMatrix m_Qs;
    TridiagonalMatrix m_Ts;
    mutable Vector m_A_times_q;
    Eigen::DenseIndex m_rank;
  };
} // namespace pinocchio

#endif // ifndef __pinocchio_math_lanczos_decomposition_hpp__
