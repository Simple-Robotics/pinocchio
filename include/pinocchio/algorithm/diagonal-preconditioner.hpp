//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_diagonal_preconditioner_hpp__
#define __pinocchio_algorithm_diagonal_preconditioner_hpp__

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
  template<typename VectorLike>
  struct DiagonalPreconditioner : PreconditionerBase<DiagonalPreconditioner<VectorLike>>
  {

    /// \brief Default constructor takes a vector.
    /// @param diagonal Vector composing the diagonal of the preconditioner
    template<typename InputVector>
    explicit DiagonalPreconditioner(const Eigen::MatrixBase<InputVector> & diagonal)
    : m_diagonal(diagonal)
    , m_squared_diagonal(diagonal)
    {
      typedef typename VectorLike::Scalar Scalar;
      PINOCCHIO_CHECK_INPUT_ARGUMENT((diagonal.array() >= Scalar(0)).all());
      m_squared_diagonal.array() *= diagonal.array();
    }

    /// @brief Default constructor from a given size.
    /// @param size Size of the preconditioner
    explicit DiagonalPreconditioner(const Eigen::Index size)
    : m_diagonal(VectorLike::Ones(size))
    , m_squared_diagonal(VectorLike::Ones(size))
    {
    }

    /// \brief Construct an identity preconditioner
    /// @param size Size of the preconditioner
    static DiagonalPreconditioner
    Identity(const Eigen::Index size){return DiagonalPreconditioner(size)}

    /// \brief Move constructor
    DiagonalPreconditioner(DiagonalPreconditioner && other)
    : m_diagonal(std::move(other.m_diagonal))
    , m_squared_diagonal(std::move(other.m_squared_diagonal))
    {
    }

    /// \brief Performs the scale operation to go from x to x_bar: x_bar = P^{-1} * x.
    template<typename MatrixIn, typename MatrixOut>
    void
    scale(const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & x_bar) const
    {
      auto & x_bar_ = x_bar.const_cast_derived();
      x_bar_ = x;
      scaleInPlace(x_bar_);
    }

    /// \brief see \ref scale
    template<typename MatrixIn>
    void scaleInPlace(const Eigen::MatrixBase<MatrixIn> & x) const
    {
      auto & x_ = x.const_cast_derived();
      x_.array() = x.array() / m_diagonal.array();
    }

    /// \brief Performs the unscale operation to go from x_bar to x: x = P * x_bar.
    template<typename MatrixIn, typename MatrixOut>
    void
    unscale(const Eigen::MatrixBase<MatrixIn> & x_bar, const Eigen::MatrixBase<MatrixOut> & x) const
    {
      auto & x_ = x.const_cast_derived();
      x_ = x_bar;
      unscaleInPlace(x_);
    }

    /// \brief see \ref \unscale
    template<typename MatrixIn>
    void unscaleInPlace(const Eigen::MatrixBase<MatrixIn> & x) const
    {
      auto & x_ = x.const_cast_derived();
      x_.array() *= m_diagonal.array();
    }

    /// \brief Performs the scale operation twice: x_bar = P^{-2} * x.
    template<typename MatrixIn, typename MatrixOut>
    void scaleSquare(
      const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & x_bar) const
    {
      auto & x_bar_ = x_bar.const_cast_derived();
      x_bar_.array() = x.array() / m_squared_diagonal.array();
    }

    /// \brief Performs the unscale operation twice: x = P * x_bar.
    template<typename MatrixIn, typename MatrixOut>
    void unscaleSquare(
      const Eigen::MatrixBase<MatrixIn> & x_bar, const Eigen::MatrixBase<MatrixOut> & x) const
    {
      auto & x_ = x.const_cast_derived();
      x_.array() = x_bar.array() * m_squared_diagonal.array();
    }

    Eigen::DenseIndex rows() const
    {
      return m_diagonal.size();
    }
    Eigen::DenseIndex cols() const
    {
      return m_diagonal.size();
    }

    void setDiagonal(const VectorLike & x)
    {
      m_diagonal = x;
      m_squared_diagonal.array() = x.array() * x.array();
    }

    const VectorLike & getDiagonal() const
    {
      return m_diagonal;
    }

  protected:
    VectorLike m_diagonal;
    VectorLike m_squared_diagonal;

  }; // struct DiagonalPreconditioner

} // namespace pinocchio

#endif // #ifndef __pinocchio_algorithm_diagonal_preconditioner_hpp__
