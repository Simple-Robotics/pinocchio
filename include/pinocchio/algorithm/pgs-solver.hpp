//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_algorithm_pgs_solver_hpp__
#define __pinocchio_algorithm_pgs_solver_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/contact-solver-base.hpp"

namespace pinocchio
{

  template<typename Scalar>
  struct PGSContactSolverTpl;
  typedef PGSContactSolverTpl<context::Scalar> PGSContactSolver;

  /// \brief Projected Gauss Siedel solver
  template<typename _Scalar>
  struct PGSContactSolverTpl : ContactSolverBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ContactSolverBaseTpl<Scalar> Base;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;

    explicit PGSContactSolverTpl(const int problem_size)
    : Base(problem_size)
    , x(problem_size)
    , x_previous(problem_size)
    , y(problem_size)
    {
    }

    ///
    /// \brief Solve the constrained problem composed of problem data (G,g,constraint_sets) and
    /// starting from the initial guess.
    ///
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the contact problem.
    /// \param[in] g Free contact acceleration or velicity associted with the contact problem.
    /// \param[in] constraint_models Vector of constraint models.
    /// \param[in] x Initial guess solution of the problem.
    /// \param[in] over_relax Optional over relaxation value, default to 1.
    ///
    /// \returns True if the problem has converged.
    template<
      typename MatrixLike,
      typename VectorLike,
      template<typename T>
      class Holder,
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename VectorLikeOut>
    bool solve(
      const MatrixLike & G,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<
      Holder<const ConstraintModel>,
      ConstraintModelAllocator> & constraint_models,
      // const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const Eigen::DenseBase<VectorLikeOut> & x,
      const Scalar over_relax = Scalar(1));

    ///
    /// \brief Solve the constrained problem composed of problem data (G,g,constraint_sets) and
    /// starting from the initial guess.
    ///
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the contact problem.
    /// \param[in] g Free contact acceleration or velicity associted with the contact problem.
    /// \param[in] constraint_models Vector of constraint models.
    /// \param[in,out] x Initial guess and output solution of the problem
    /// \param[in] over_relax Over relaxation value
    ///
    /// \returns True if the problem has converged.
    template<
      typename MatrixLike,
      typename VectorLike,
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename VectorLikeInitialGuess>
    bool solve(
      const MatrixLike & G,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const Eigen::DenseBase<VectorLikeInitialGuess> & x,
      const Scalar over_relax = Scalar(1));

    /// \returns the primal solution of the problem
    const VectorXs & getPrimalSolution() const
    {
      return x;
    }
    /// \returns the dual solution of the problem
    const VectorXs & getDualSolution() const
    {
      return y;
    }

  protected:
    /// \brief Previous temporary value of the optimum.
    VectorXs x, x_previous;
    /// \brief Value of the dual variable
    VectorXs y;
#ifdef PINOCCHIO_WITH_HPP_FCL
    using Base::timer;
#endif // PINOCCHIO_WITH_HPP_FCL

  }; // struct PGSContactSolverTpl
} // namespace pinocchio

#include "pinocchio/algorithm/pgs-solver.hxx"

#endif // ifndef __pinocchio_algorithm_pgs_solver_hpp__
