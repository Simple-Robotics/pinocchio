//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_algorithm_admm_solver_hpp__
#define __pinocchio_algorithm_admm_solver_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/math/fwd.hpp"
#include "pinocchio/math/comparison-operators.hpp"
#include "pinocchio/math/eigenvalues.hpp"

#include "pinocchio/algorithm/contact-solver-base.hpp"
#include "pinocchio/algorithm/delassus-operator-base.hpp"

#include "pinocchio/math/lanczos-decomposition.hpp"

#include <boost/optional.hpp>

namespace pinocchio
{
  template<typename Scalar>
  struct ADMMContactSolverTpl;
  typedef ADMMContactSolverTpl<context::Scalar> ADMMContactSolver;

  template<typename _Scalar>
  struct ADMMSpectralUpdateRuleTpl
  {
    typedef _Scalar Scalar;

    ADMMSpectralUpdateRuleTpl(
      const Scalar ratio_primal_dual, const Scalar L, const Scalar m, const Scalar rho_power_factor)
    : ratio_primal_dual(ratio_primal_dual)
    , rho_increment(std::pow(L / m, rho_power_factor))
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(m > Scalar(0), "m should be positive.");
    }

    Scalar getRatioPrimalDual() const
    {
      return ratio_primal_dual;
    }
    void setRatioPrimalDual(const Scalar ratio_primal_dual)
    {
      this->ratio_primal_dual = ratio_primal_dual;
    }

    Scalar getRhoIncrement() const
    {
      return rho_increment;
    }
    void setRhoIncrement(const Scalar cond, const Scalar rho_power_factor)
    {
      rho_increment = std::pow(cond, rho_power_factor);
    }

    bool eval(const Scalar primal_feasibility, const Scalar dual_feasibility, Scalar & rho) const
    {
      bool rho_has_changed = false;
      if (primal_feasibility > ratio_primal_dual * dual_feasibility)
      {
        rho *= rho_increment;
        //        rho *= math::pow(cond,rho_power_factor);
        //        rho_power += rho_power_factor;
        rho_has_changed = true;
      }
      else if (dual_feasibility > ratio_primal_dual * primal_feasibility)
      {
        rho /= rho_increment;
        //        rho *= math::pow(cond,-rho_power_factor);
        //        rho_power -= rho_power_factor;
        rho_has_changed = true;
      }

      return rho_has_changed;
    }

    /// \brief Compute the penalty ADMM value from the current largest and lowest eigenvalues and
    /// the scaling spectral factor.
    static inline Scalar computeRho(const Scalar L, const Scalar m, const Scalar rho_power)
    {
      const Scalar cond = L / m;
      const Scalar rho = math::sqrt(L * m) * math::pow(cond, rho_power);
      return rho;
    }

    /// \brief Compute the  scaling spectral factor of the ADMM penalty term from the current
    /// largest and lowest eigenvalues and the ADMM penalty term.
    static inline Scalar computeRhoPower(const Scalar L, const Scalar m, const Scalar rho)
    {
      const Scalar cond = L / m;
      const Scalar sqtr_L_m = math::sqrt(L * m);
      const Scalar rho_power = math::log(rho / sqtr_L_m) / math::log(cond);
      return rho_power;
    }

  protected:
    Scalar ratio_primal_dual;
    Scalar rho_increment;
  };

  template<typename _Scalar>
  struct ADMMLinearUpdateRuleTpl
  {
    typedef _Scalar Scalar;

    ADMMLinearUpdateRuleTpl(
      const Scalar ratio_primal_dual, const Scalar increase_factor, const Scalar decrease_factor)
    : ratio_primal_dual(ratio_primal_dual)
    , increase_factor(increase_factor)
    , decrease_factor(decrease_factor)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        increase_factor > Scalar(1), "increase_factor should be greater than one.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        decrease_factor > Scalar(1), "decrease_factor should be greater than one.");
    }

    ADMMLinearUpdateRuleTpl(const Scalar ratio_primal_dual, const Scalar factor)
    : ratio_primal_dual(ratio_primal_dual)
    , increase_factor(factor)
    , decrease_factor(factor)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(factor > Scalar(1), "factor should be greater than one.");
    }

    bool eval(const Scalar primal_feasibility, const Scalar dual_feasibility, Scalar & rho) const
    {
      bool rho_has_changed = false;
      if (primal_feasibility > ratio_primal_dual * dual_feasibility)
      {
        rho *= increase_factor;
        rho_has_changed = true;
      }
      else if (dual_feasibility > ratio_primal_dual * primal_feasibility)
      {
        rho /= decrease_factor;
        rho_has_changed = true;
      }

      return rho_has_changed;
    }

  protected:
    Scalar ratio_primal_dual;
    Scalar increase_factor, decrease_factor;
  };

  enum class ADMMUpdateRule : char
  {
    SPECTRAL = 'S',
    LINEAR = 'L',
  };

  template<typename Scalar>
  union ADMMUpdateRuleContainerTpl {
    ADMMUpdateRuleContainerTpl()
    : dummy(){};
    ADMMSpectralUpdateRuleTpl<Scalar> spectral_rule;
    ADMMLinearUpdateRuleTpl<Scalar> linear_rule;

  protected:
    struct Dummy
    {
      Dummy(){};
    };

    Dummy dummy{};
  };

  template<typename _Scalar>
  struct PINOCCHIO_UNSUPPORTED_MESSAGE("The API will change towards more flexibility")
    ADMMContactSolverTpl : ContactSolverBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ContactSolverBaseTpl<_Scalar> Base;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef const Eigen::Ref<const VectorXs> ConstRefVectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
    typedef LanczosDecompositionTpl<VectorXs> LanczosAlgo;

    using Base::problem_size;

    //    struct SolverParameters
    //    {
    //      explicit SolverParameters(const int problem_dim)
    //      : rho_power(Scalar(0.2))
    //      , ratio_primal_dual(Scalar(10))
    //      , mu_prox
    //      {
    //
    //      }
    //
    //      /// \brief Rho solver ADMM
    //      boost::optional<Scalar> rho;
    //      /// \brief Power value associated to rho. This quantity will be automatically updated.
    //      Scalar rho_power;
    //      /// \brief Ratio primal/dual
    //      Scalar ratio_primal_dual;
    //      /// \brief Proximal value
    //      Scalar mu_prox;
    //
    //      /// \brief Largest eigenvalue
    //      boost::optional<Scalar> L_value;
    //      /// \brief Largest eigenvector
    //      boost::optional<VectorXs> L_vector;
    //    };
    //
    struct SolverStats
    {
      explicit SolverStats(const int max_it)
      : it(0)
      , cholesky_update_count(0)
      {
        primal_feasibility.reserve(size_t(max_it));
        dual_feasibility.reserve(size_t(max_it));
        dual_feasibility_ncp.reserve(size_t(max_it));
        complementarity.reserve(size_t(max_it));
        rho.reserve(size_t(max_it));
      }

      void reset()
      {
        primal_feasibility.clear();
        dual_feasibility.clear();
        complementarity.clear();
        dual_feasibility_ncp.clear();
        rho.clear();
        it = 0;
        cholesky_update_count = 0;
      }

      size_t size() const
      {
        return primal_feasibility.size();
      }

      ///  \brief Number of total iterations.
      int it;

      ///  \brief Number of Cholesky updates.
      int cholesky_update_count;

      /// \brief History of primal feasibility values.
      std::vector<Scalar> primal_feasibility;

      /// \brief History of dual feasibility values.
      std::vector<Scalar> dual_feasibility;
      std::vector<Scalar> dual_feasibility_ncp;

      /// \brief History of complementarity values.
      std::vector<Scalar> complementarity;

      /// \brief History of rho values.
      std::vector<Scalar> rho;
    };
    //
    //    struct SolverResults
    //    {
    //      explicit SolverResults(const int problem_dim, const int max_it)
    //      : L_vector(problem_dim)
    //
    //      /// \brief Largest eigenvalue
    //      Scalar L_value;
    //      /// \brief Largest eigenvector
    //      VectorXs L_vector;
    //
    //      SolverStats stats;
    //    };

    explicit ADMMContactSolverTpl(
      int problem_dim,
      Scalar mu_prox = Scalar(1e-6),
      Scalar tau = Scalar(0.5),
      Scalar rho_power = Scalar(0.2),
      Scalar rho_power_factor = Scalar(0.05),
      Scalar linear_update_rule_factor = Scalar(2),
      Scalar ratio_primal_dual = Scalar(10),
      int max_it_largest_eigen_value_solver = 20)
    : Base(problem_dim)
    , is_initialized(false)
    , mu_prox(mu_prox)
    , tau(tau)
    , rho(10.)
    , rho_power(rho_power)
    , rho_power_factor(rho_power_factor)
    , linear_update_rule_factor(linear_update_rule_factor)
    , ratio_primal_dual(ratio_primal_dual)
    , max_it_largest_eigen_value_solver(max_it_largest_eigen_value_solver)
    , x_(VectorXs::Zero(problem_dim))
    , y_(VectorXs::Zero(problem_dim))
    , x_previous(VectorXs::Zero(problem_dim))
    , y_previous(VectorXs::Zero(problem_dim))
    , z_previous(VectorXs::Zero(problem_dim))
    , z_(VectorXs::Zero(problem_dim))
    , s_(VectorXs::Zero(problem_dim))
    , rhs(problem_dim)
    , primal_feasibility_vector(VectorXs::Zero(problem_dim))
    , dual_feasibility_vector(VectorXs::Zero(problem_dim))
    , stats(Base::max_it)
    {
    }

    /// \brief Set the ADMM penalty value.
    void setRho(const Scalar rho)
    {
      this->rho = rho;
    }
    /// \brief Get the ADMM penalty value.
    Scalar getRho() const
    {
      return rho;
    }

    /// \brief Set the power associated to the problem conditionning.
    void setRhoPower(const Scalar rho_power)
    {
      this->rho_power = rho_power;
    }
    /// \brief Get the power associated to the problem conditionning.
    Scalar getRhoPower() const
    {
      return rho_power;
    }

    /// \brief Set the power factor associated to the problem conditionning.
    void setRhoPowerFactor(const Scalar rho_power_factor)
    {
      this->rho_power_factor = rho_power_factor;
    }
    /// \brief Get the value of the increase/decrease factor associated to the problem
    /// conditionning.
    Scalar getRhoPowerFactor() const
    {
      return rho_power_factor;
    }

    /// \brief Set the update factor of the Linear update rule
    void setLinearUpdateRuleFactor(const Scalar linear_update_rule_factor)
    {
      this->linear_update_rule_factor = linear_update_rule_factor;
    }
    /// \brief Get the value of the increase/decrease factor of the Linear update rule
    Scalar getLinearUpdateRuleFactor() const
    {
      return linear_update_rule_factor;
    }

    /// \brief Set the tau linear scaling factor.
    void setTau(const Scalar tau)
    {
      this->tau = tau;
    }
    /// \brief Get the tau linear scaling factor.
    Scalar getTau() const
    {
      return tau;
    }

    /// \brief Set the proximal value.
    void setProximalValue(const Scalar mu)
    {
      this->mu_prox = mu;
    }
    /// \brief Get the proximal value.
    Scalar getProximalValue() const
    {
      return mu_prox;
    }

    /// \brief Set the primal/dual ratio.
    void setRatioPrimalDual(const Scalar ratio_primal_dual)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        ratio_primal_dual > 0., "The ratio primal/dual should be positive strictly");
      this->ratio_primal_dual = ratio_primal_dual;
    }
    /// \brief Get the primal/dual ratio.
    Scalar getRatioPrimalDual() const
    {
      return ratio_primal_dual;
    }

    ///  \returns the number of updates of the Cholesky factorization due to rho updates.
    int getCholeskyUpdateCount() const
    {
      return cholesky_update_count;
    }

    ///
    /// \brief Solve the constraint problem composed of problem data (G,g,constraint_models) and
    /// starting from the initial guess.
    ///
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the constraint problem.
    /// \param[in] g Free constraint acceleration or velicity associted with the constraint problem.
    /// \param[in] constraint_models Vector of constraints.
    /// \param[in] R Proximal regularization value associated to the compliant constraints
    /// (corresponds to the lowest non-zero). \param[in] primal_guess Optional initial guess of the
    /// primal solution (constrained forces). \param[in] dual_guess Optinal Initial guess of the
    /// dual solution (constrained velocities). \param[in] solve_ncp whether to solve the NCP (true)
    /// or CCP (false) \param[in] compute_largest_eigen_values run power iteration algorithm
    /// \param[in] admm_update_rule update rule for ADMM (linear or spectral)
    /// \param[in] stat_record record solver metrics
    ///
    /// \returns True if the problem has converged.
    template<
      typename DelassusDerived,
      typename VectorLike,
      template<typename T>
      class Holder,
      typename ConstraintModel,
      typename ConstraintAllocator,
      typename VectorLikeR>
    bool solve(
      DelassusOperatorBase<DelassusDerived> & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<Holder<const ConstraintModel>, ConstraintAllocator> & constraint_models,
      const Eigen::MatrixBase<VectorLikeR> & R,
      const boost::optional<ConstRefVectorXs> primal_guess = boost::none,
      const boost::optional<ConstRefVectorXs> dual_guess = boost::none,
      bool solve_ncp = true,
      bool compute_largest_eigen_values = true,
      ADMMUpdateRule admm_update_rule = ADMMUpdateRule::SPECTRAL,
      bool stat_record = false);

    ///
    /// \brief Solve the constraint problem composed of problem data (G,g,constraint_models) and
    /// starting from the initial guess.
    ///
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the constraint problem.
    /// \param[in] g Free constraint acceleration or velicity associted with the constraint problem.
    /// \param[in] constraint_models Vector of constraints.
    /// \param[in] R Proximal regularization value associated to the compliant constraints
    /// (corresponds to the lowest non-zero). \param[in] primal_guess Optional initial guess of the
    /// primal solution (constrained forces). \param[in] dual_guess Optinal Initial guess of the
    /// dual solution (constrained velocities). \param[in] solve_ncp whether to solve the NCP (true)
    /// or CCP (false) \param[in] compute_largest_eigen_values run power iteration algorithm
    /// \param[in] admm_update_rule update rule for ADMM (linear or spectral)
    /// \param[in] stat_record record solver metrics
    ///
    /// \returns True if the problem has converged.
    template<
      typename DelassusDerived,
      typename VectorLike,
      typename ConstraintModel,
      typename ConstraintAllocator,
      typename VectorLikeR>
    bool solve(
      DelassusOperatorBase<DelassusDerived> & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<ConstraintModel, ConstraintAllocator> & constraint_models,
      const Eigen::MatrixBase<VectorLikeR> & R,
      const boost::optional<ConstRefVectorXs> primal_guess = boost::none,
      const boost::optional<ConstRefVectorXs> dual_guess = boost::none,
      bool solve_ncp = true,
      bool compute_largest_eigen_values = true,
      ADMMUpdateRule admm_update_rule = ADMMUpdateRule::SPECTRAL,
      bool stat_record = false)
    {
      typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
      typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

      WrappedConstraintModelVector wrapped_constraint_models(
        constraint_models.cbegin(), constraint_models.cend());

      return solve(
        delassus, g, wrapped_constraint_models, R, primal_guess, dual_guess, solve_ncp,
        compute_largest_eigen_values, admm_update_rule, stat_record);
    }

    ///
    /// \brief Solve the constraint problem composed of problem data (G,g,constraint_models) and
    /// starting from the initial guess.
    ///
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the constraint problem.
    /// \param[in] g Free constraint acceleration or velicity associted with the constraint problem.
    /// \param[in] constraint_models Vector of constraints.
    /// \param[in] primal_guess Optional initial guess of the primal solution (constrained forces).
    /// \param[in] solve_ncp whether to solve the NCP (true) or CCP (false)
    ///
    /// \returns True if the problem has converged.
    template<
      typename DelassusDerived,
      typename VectorLike,
      template<typename T>
      class Holder,
      typename ConstraintModel,
      typename ConstraintAllocator,
      typename VectorLikeOut>
    bool solve(
      DelassusOperatorBase<DelassusDerived> & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<Holder<const ConstraintModel>, ConstraintAllocator> & constraint_models,
      const Eigen::DenseBase<VectorLikeOut> & primal_guess,
      bool solve_ncp = true)
    {
      return solve(
        delassus.derived(), g.derived(), constraint_models, VectorXs::Zero(problem_size),
        primal_guess.const_cast_derived(), boost::none, solve_ncp);
    }

    ///
    /// \brief Solve the constraint problem composed of problem data (G,g,constraint_models) and
    /// starting from the initial guess.
    ///
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the constraint problem.
    /// \param[in] g Free constraint acceleration or velicity associted with the constraint problem.
    /// \param[in] constraint_models Vector of constraints.
    /// \param[in] primal_guess Optional initial guess of the primal solution (constrained forces).
    /// \param[in] solve_ncp whether to solve the NCP (true) or CCP (false)
    ///
    /// \returns True if the problem has converged.
    template<
      typename DelassusDerived,
      typename VectorLike,
      typename ConstraintModel,
      typename ConstraintAllocator,
      typename VectorLikeOut>
    bool solve(
      DelassusOperatorBase<DelassusDerived> & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<ConstraintModel, ConstraintAllocator> & constraint_models,
      const Eigen::DenseBase<VectorLikeOut> & primal_guess,
      bool solve_ncp = true)
    {
      typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
      typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

      WrappedConstraintModelVector wrapped_constraint_models(
        constraint_models.cbegin(), constraint_models.cend());

      return solve(delassus, g, wrapped_constraint_models, primal_guess, solve_ncp);
    }

    /// \returns the primal solution of the problem
    const VectorXs & getPrimalSolution() const
    {
      return y_;
    }
    /// \returns the dual solution of the problem
    const VectorXs & getDualSolution() const
    {
      return z_;
    }
    /// \returns the complementarity shift
    const VectorXs & getComplementarityShift() const
    {
      return s_;
    }

    SolverStats & getStats()
    {
      return stats;
    }

  protected:
    bool is_initialized;

    /// \brief proximal value
    Scalar mu_prox;

    /// \brief Linear scaling of the ADMM penalty term
    Scalar tau;
    /// \brief Penalty term associated to the ADMM.
    Scalar rho;

    // Set of parameters associated with the Spectral update rule
    /// \brief Power value associated to rho. This quantity will be automatically updated.
    Scalar rho_power;
    /// \brief Update factor for the primal/dual update of rho.
    Scalar rho_power_factor;

    // Set of parameters associated with the Linear update rule
    /// \brief value of the increase/decrease factor
    Scalar linear_update_rule_factor;

    ///  \brief Ratio primal/dual
    Scalar ratio_primal_dual;

    /// \brief Maximum number of iterarions called for the power iteration algorithm
    int max_it_largest_eigen_value_solver;

    /// \brief Primal variables (corresponds to the contact forces)
    VectorXs x_, y_;
    /// \brief Previous value of y.
    VectorXs x_previous, y_previous, z_previous;
    /// \brief Dual varible of the ADMM (corresponds to the contact velocity or acceleration).
    VectorXs z_;
    /// \brief De Saxé shift
    VectorXs s_;

    VectorXs rhs, primal_feasibility_vector, dual_feasibility_vector;

    int cholesky_update_count;

    /// \brief Stats of the solver
    SolverStats stats;

#ifdef PINOCCHIO_WITH_HPP_FCL
    using Base::timer;
#endif // PINOCCHIO_WITH_HPP_FCL
  }; // struct ADMMContactSolverTpl
} // namespace pinocchio

#include "pinocchio/algorithm/admm-solver.hxx"

#endif // ifndef __pinocchio_algorithm_admm_solver_hpp__
