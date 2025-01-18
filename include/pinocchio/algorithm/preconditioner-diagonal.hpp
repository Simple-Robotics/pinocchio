//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_preconditioner_diagonal_hpp__
#define __pinocchio_algorithm_preconditioner_diagonal_hpp__

#include "pinocchio/algorithm/preconditioner-base.hpp"
namespace pinocchio
{

  ///
  /// \brief Pre-Conditioning a constraint problem is done purely for numerical reasons.
  /// We have a problem that looks like $$ min_{x \in K} x^T G x + g^T x $$.
  /// The pre-conditioner is a diagonal matrix P.
  /// We write x = P * x_bar (hence x_bar = P^{-1} * x), g_bar = Pg and G_bar = P*G*P, such that
  /// the problem now becomes:
  /// $$ min_{x_bar \in K} x_bar^T G_bar x + g_bar^T x $$.
  //
  /// \note We call the original problem working on (x, g, G) the **unscaled** problem.
  /// We call the new problem working on (x_bar, g_bar, G_bar) the **scaled** problem.
  ///
  template<typename PreconditionerVectorLike>
  struct PreconditionerDiagonal
  : PreconditionerBase<PreconditionerDiagonal<PreconditionerVectorLike>>
  {

    /// \brief Default constructor takes a vector.
    explicit PreconditionerDiagonal(const PreconditionerVectorLike & diagonal)
    : m_preconditioner_diagonal(diagonal)
    , m_preconditioner_square(diagonal)
    {
      m_preconditioner_square.array() *= diagonal.array();
    }

    /// \brief Performs the scale operation to go from x to x_bar: x_bar = P^{-1} * x.
    template<typename MatrixIn, typename MatrixOut>
    void
    scale(const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & x_bar) const
    {
      auto & x_bar_ = x_bar.const_cast_derived();
      x_bar_.array() = x.array() / m_preconditioner_diagonal.array();
    }

    /// \brief see \ref scale
    template<typename MatrixIn>
    void scaleInPlace(const Eigen::MatrixBase<MatrixIn> & x) const
    {
      auto & x_ = x.const_cast_derived();
      x_.array() = x.array() / m_preconditioner_diagonal.array();
    }

    /// \brief Performs the unscale operation to go from x_bar to x: x = P * x_bar.
    template<typename MatrixIn, typename MatrixOut>
    void
    unscale(const Eigen::MatrixBase<MatrixIn> & x_bar, const Eigen::MatrixBase<MatrixOut> & x) const
    {
      auto & x_ = x.const_cast_derived();
      x_.array() = x_bar.array() * m_preconditioner_diagonal.array();
    }

    /// \brief see \ref \unscale
    template<typename MatrixIn>
    void unscaleInPlace(const Eigen::MatrixBase<MatrixIn> & x) const
    {
      auto & x_ = x.const_cast_derived();
      x_.array() *= m_preconditioner_diagonal.array();
    }

    /// \brief Performs the scale operation twice: x_bar = P^{-2} * x.
    template<typename MatrixIn, typename MatrixOut>
    void scaleSquare(
      const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & x_bar) const
    {
      auto & x_bar_ = x_bar.const_cast_derived();
      x_bar_.array() = x.array() / m_preconditioner_square.array();
    }

    /// \brief Performs the unscale operation twice: x = P * x_bar.
    template<typename MatrixIn, typename MatrixOut>
    void unscaleSquare(
      const Eigen::MatrixBase<MatrixIn> & x_bar, const Eigen::MatrixBase<MatrixOut> & x) const
    {
      auto & x_ = x.const_cast_derived();
      x_.array() = x_bar.array() * m_preconditioner_square.array();
    }

    Eigen::DenseIndex rows() const
    {
      return m_preconditioner_diagonal.size();
    }
    Eigen::DenseIndex cols() const
    {
      return m_preconditioner_diagonal.size();
    }

    void setDiagonal(const PreconditionerVectorLike & x)
    {
      m_preconditioner_diagonal = x;
      m_preconditioner_square.array() = x.array() * x.array();
    }

    const PreconditionerVectorLike & getDiagonal() const
    {
      return m_preconditioner_diagonal;
    }

  protected:
    PreconditionerVectorLike m_preconditioner_diagonal;
    PreconditionerVectorLike m_preconditioner_square;

  }; // struct PreconditionerDiagonal

} // namespace pinocchio

#endif // #ifndef __pinocchio_algorithm_preconditioner_diagonal_hpp__
