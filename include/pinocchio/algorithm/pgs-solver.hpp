//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_algorithm_pgs_solver_hpp__
#define __pinocchio_algorithm_pgs_solver_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/contact-solver-base.hpp"
#include "pinocchio/algorithm/delassus-operator-dense.hpp"
#include <boost/optional.hpp>

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
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Ref<const VectorXs> RefConstVectorXs;

    typedef typename Base::SolverStats SolverStats;

    explicit PGSContactSolverTpl(const int problem_size)
    : Base(problem_size)
    , x(problem_size)
    , x_previous(problem_size)
    , y(problem_size)
    , y_to_pos(problem_size)
    , time_scaling_acc_to_constraints(VectorXs::Zero(problem_size))
    , time_scaling_constraints_to_pos(VectorXs::Zero(problem_size))
    , gs(VectorXs::Zero(problem_size))
    , stats()
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
      typename VectorLike,
      template<typename T> class Holder,
      typename ConstraintModel,
      typename ConstraintModelAllocator>
    bool solve(
      const DelassusOperatorDense & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<Holder<const ConstraintModel>, ConstraintModelAllocator> &
        constraint_models,
      const Scalar dt,
      const boost::optional<RefConstVectorXs> x_guess = boost::none,
      const Scalar over_relax = Scalar(1),
      const bool solve_ncp = true,
      const bool stat_record = false);

    ///
    /// \brief Solve the constrained problem composed of problem data (G,g,constraint_sets) and
    /// starting from the initial guess.
    ///
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the contact problem.
    /// \param[in] g Free contact acceleration or velicity associted with the contact problem.
    /// \param[in] constraint_models Vector of constraint models.
    /// \param[in] x_guess Initial guess and output solution of the problem
    /// \param[in] over_relax Over relaxation value
    ///
    /// \returns True if the problem has converged.
    template<typename VectorLike, typename ConstraintModel, typename ConstraintModelAllocator>
    bool solve(
      const DelassusOperatorDense & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const Scalar dt,
      const boost::optional<RefConstVectorXs> x_guess = boost::none,
      const Scalar over_relax = Scalar(1),
      const bool solve_ncp = true,
      const bool stat_record = false);

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

    SolverStats & getStats()
    {
      return stats;
    }

  protected:
    /// \brief Previous temporary value of the optimum.
    VectorXs x, x_previous;
    /// \brief Value of the dual variable
    VectorXs y, y_to_pos;

    /// \brief Time scaling vector for constraints
    VectorXs time_scaling_acc_to_constraints, time_scaling_constraints_to_pos;
    /// \brief Vector g divided by time scaling (g / time_scaling_acc_to_constraints)
    VectorXs gs;

#ifdef PINOCCHIO_WITH_HPP_FCL
    using Base::timer;
#endif // PINOCCHIO_WITH_HPP_FCL

    /// \brief Stats of the solver
    SolverStats stats;

  }; // struct PGSContactSolverTpl
} // namespace pinocchio

#include "pinocchio/algorithm/pgs-solver.hxx"

#endif // ifndef __pinocchio_algorithm_pgs_solver_hpp__
