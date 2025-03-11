//
// Copyright (c) 2022-2025 INRIA
//

#ifndef __pinocchio_algorithm_admm_solver_hxx__
#define __pinocchio_algorithm_admm_solver_hxx__

#include <limits>

#include "pinocchio/algorithm/contact-solver-utils.hpp"
#include "pinocchio/algorithm/constraints/coulomb-friction-cone.hpp"
#include "pinocchio/algorithm/constraints/visitors/constraint-model-visitor.hpp"
#include "pinocchio/algorithm/delassus-operator-preconditioned.hpp"

#include "pinocchio/tracy.hpp"

namespace pinocchio
{

  template<typename DriftVectorLike, typename Scalar>
  struct ZeroInitialGuessMaxConstraintViolationVisitor
  : visitors::ConstraintUnaryVisitorBase<
      ZeroInitialGuessMaxConstraintViolationVisitor<DriftVectorLike, Scalar>>
  {
    using ArgsType = boost::fusion::vector<const DriftVectorLike &, Scalar &>;
    using Base = visitors::ConstraintUnaryVisitorBase<
      ZeroInitialGuessMaxConstraintViolationVisitor<DriftVectorLike, Scalar>>;

    template<typename ConstraintModel>
    static void algo(
      const ConstraintModelBase<ConstraintModel> & cmodel,
      const DriftVectorLike & drift,
      Scalar & max_violation)
    {
      return algo_impl(cmodel.set(), drift, max_violation);
    }

    template<typename VectorLike>
    static void algo_impl(
      const CoulombFrictionConeTpl<Scalar> & set,
      const Eigen::MatrixBase<VectorLike> & drift,
      Scalar & max_violation)
    {
      PINOCCHIO_UNUSED_VARIABLE(set);
      const Scalar violation = -drift.coeff(2);
      if (violation > max_violation)
      {
        max_violation = violation;
      }
    }

    template<typename VectorLike>
    static void algo_impl(
      const UnboundedSetTpl<Scalar> & set,
      const Eigen::MatrixBase<VectorLike> & drift,
      Scalar & max_violation)
    {
      PINOCCHIO_UNUSED_VARIABLE(set);
      const Scalar violation = drift.template lpNorm<Eigen::Infinity>();
      if (violation > max_violation)
      {
        max_violation = violation;
      }
    }

    template<typename VectorLike>
    static void algo_impl(
      const BoxSetTpl<Scalar> & set,
      const Eigen::MatrixBase<VectorLike> & drift,
      Scalar & max_violation)
    {
      PINOCCHIO_UNUSED_VARIABLE(set);
      const Scalar violation = drift.template lpNorm<Eigen::Infinity>();
      if (violation > max_violation)
      {
        max_violation = violation;
      }
    }

    template<typename VectorLike>
    static void algo_impl(
      const JointLimitConstraintConeTpl<Scalar> & set,
      const Eigen::MatrixBase<VectorLike> & drift,
      Scalar & max_violation)
    {
      Scalar negative_orthant_violation = 0;
      if (set.getNegativeOrthant().size() > 0)
      {
        negative_orthant_violation =
          math::max(Scalar(0), drift.head(set.getNegativeOrthant().size()).maxCoeff());
      }
      Scalar positive_orthant_violation = 0;
      if (set.getPositiveOrthant().size() > 0)
      {
        positive_orthant_violation =
          -math::min(Scalar(0), drift.tail(set.getPositiveOrthant().size()).minCoeff());
      }
      const Scalar violation = math::max(negative_orthant_violation, positive_orthant_violation);
      if (violation > max_violation)
      {
        max_violation = violation;
      }
    }

    template<typename ConstraintSet, typename VectorLike>
    static void algo_impl(
      const ConstraintSet & set,
      const Eigen::MatrixBase<VectorLike> & drift,
      Scalar & max_violation)
    {
      // do nothing
      PINOCCHIO_UNUSED_VARIABLE(set);
      PINOCCHIO_UNUSED_VARIABLE(drift);
      PINOCCHIO_UNUSED_VARIABLE(max_violation);
    }

    /// ::run for individual constraints
    template<typename ConstraintModel>
    static void run(
      const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
      const DriftVectorLike & drift,
      Scalar & max_violation)
    {
      algo(cmodel.derived(), drift, max_violation);
    }

    /// ::run for constraints variant
    template<int Options, template<typename S, int O> class ConstraintCollectionTpl>
    static void run(
      const pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const DriftVectorLike & drift,
      Scalar & max_violation)
    {
      ArgsType args(drift, max_violation);
      // Note: Base::run will call `algo` of this visitor
      Base::run(cmodel.derived(), args);
    }
  }; // struct ZeroInitialGuessMaxConstraintViolationVisitor

  template<
    template<typename T> class Holder,
    typename ConstraintModel,
    typename ConstraintModelAllocator,
    typename VectorLikeIn>
  typename ConstraintModel::Scalar computeZeroInitialGuessMaxConstraintViolation(
    const std::vector<Holder<const ConstraintModel>, ConstraintModelAllocator> & constraint_models,
    const Eigen::DenseBase<VectorLikeIn> & drift)
  {
    PINOCCHIO_TRACY_ZONE_SCOPED_N("computeZeroInitialGuessMaxConstraintViolation");
    Eigen::DenseIndex cindex = 0;

    using SegmentType = typename VectorLikeIn::ConstSegmentReturnType;
    using Scalar = typename ConstraintModel::Scalar;

    Scalar max_violation = Scalar(0);
    for (const ConstraintModel & cmodel : constraint_models)
    {
      const auto csize = cmodel.size();

      SegmentType drift_segment = drift.segment(cindex, csize);
      typedef ZeroInitialGuessMaxConstraintViolationVisitor<SegmentType, Scalar> Algo;

      Algo::run(cmodel, drift_segment, max_violation);

      cindex += csize;
    }
    return max_violation;
  }

  template<typename _Scalar>
  template<
    typename DelassusDerived,
    typename VectorLike,
    template<typename T> class Holder,
    typename ConstraintModel,
    typename ConstraintModelAllocator>
  bool ADMMContactSolverTpl<_Scalar>::solve(
    DelassusOperatorBase<DelassusDerived> & _delassus,
    const Eigen::MatrixBase<VectorLike> & g,
    const std::vector<Holder<const ConstraintModel>, ConstraintModelAllocator> & constraint_models,
    const Scalar dt,
    const boost::optional<RefConstVectorXs> preconditioner,
    const boost::optional<RefConstVectorXs> primal_guess,
    const boost::optional<RefConstVectorXs> dual_guess,
    const bool solve_ncp,
    const ADMMUpdateRule admm_update_rule,
    const bool stat_record)

  {
    // Unused for now
    PINOCCHIO_UNUSED_VARIABLE(dual_guess);

    using namespace internal;
    typedef ADMMSpectralUpdateRuleTpl<Scalar> ADMMSpectralUpdateRule;
    typedef ADMMLinearUpdateRuleTpl<Scalar> ADMMLinearUpdateRule;

    typedef ADMMUpdateRuleContainerTpl<Scalar> ADMMUpdateRuleContainer;

    typedef DelassusOperatorPreconditionedTpl<DelassusDerived, DiagonalPreconditioner>
      DelassusOperatorPreconditioned;
    DelassusDerived & delassus = _delassus.derived();

    const Scalar mu_R = delassus.getCompliance().minCoeff();
    PINOCCHIO_CHECK_INPUT_ARGUMENT(dt >= Scalar(0), "dt should be positive.");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(tau <= Scalar(1) && tau > Scalar(0), "tau should lie in ]0,1].");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(mu_prox >= 0, "mu_prox should be positive.");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(mu_R >= Scalar(0), "R should be a positive vector.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(delassus.getCompliance().size(), problem_size);

    // First, we initialize the primal and dual variables
    int it = 0;
    cholesky_update_count = 0;

    Scalar complementarity, dx_bar_norm, dy_bar_norm, dz_bar_norm, //
      primal_feasibility, dual_feasibility;

    if (stat_record)
    {
      stats.reserve(this->max_it);
      stats.reset();
    }

    // Then, we get the time_scaling_acc_to_constraints T from the constraints to construct gs,
    // which is g time-scaled depending on the formulation of each constraint: gs = T^{-1} * g. The
    // idea is that if we formulate a given constraint at the position/velocity/acceleration level,
    // we want to measure constraint satisfaction for this constraint at the same
    // position/velocity/acceleration level.
    // However, to take admm steps, we work at the (force, acceleration) level for all constraints.
    // In short:
    // -> gs is used to perform optimization steps (we typically work on min_x x^TGx + gs^Tx).
    // -> time_scaling_acc_to_constraints is used to check for constraint satisfaction (we typically
    // want TGx + g = 0). It allows to go from accelerations to the units of each constraints. This
    // way, x and y are always forces expressed in N.
    // -> time_scaling_constraints_to_pos similarly allows to go from the units of the constraints
    // to positions in m. Warning: this constraints time-scaling has (a priori) nothing to do with
    // the pre-conditioner.
    getTimeScalingFromAccelerationToConstraints(
      constraint_models, dt, time_scaling_acc_to_constraints);
    getTimeScalingFromConstraintsToPosition(
      time_scaling_acc_to_constraints, dt, time_scaling_constraints_to_pos);
    gs = g.array() / time_scaling_acc_to_constraints.array();
    const Scalar g_pos_norm_inf =
      (g.cwiseProduct(time_scaling_constraints_to_pos)).template lpNorm<Eigen::Infinity>();

    // Initialize De Saxé shift to 0
    // For the CCP, there is no shift
    // For the NCP, the shift will be initialized using z
    s_.setZero();

    // Initial update of the variables
    // Init x
    if (primal_guess)
    {
      x_ = primal_guess.get();
      PINOCCHIO_CHECK_ARGUMENT_SIZE(x_.size(), problem_size);
    }
    else
    {
      x_.setZero();
    }

    // Retrieve the pre-conditioner
    if (preconditioner)
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(preconditioner_.rows(), problem_size);
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        preconditioner_.getDiagonal().minCoeff() > Scalar(0),
        "Preconditioner should be a strictly positive vector.");
      preconditioner_.setDiagonal(preconditioner.get());
    }
    else
    {
      preconditioner_.setDiagonal(VectorXs::Ones(problem_size));
    }

    // Init y
    computeConeProjection(constraint_models, x_, y_);

    // Init z -> z_ = (G + R) * y_ + g
    delassus.applyOnTheRight(y_, z_);
    z_ += gs;
    z_ -= y_.cwiseProduct(delassus.getDamping());
    if (solve_ncp)
    {
      computeDeSaxeCorrection(constraint_models, z_, s_);
      z_ += s_; // Add De Saxé shift
    }

    // Computing the convergence criterion of the initial guess
    primal_feasibility = 0; // always feasible because y is projected

    // complementarity of the initial guess
    // NB: complementarity is computed between a force y_ (in N) and z_ which unit is that of the
    // constraint formulation level.
    rhs =
      z_.array() * time_scaling_acc_to_constraints.array(); // back to constraint formulation level
    complementarity = computeConicComplementarity(constraint_models, rhs, y_);

    // dual feasibility is computed in "position" on the z_ variable (and not on z_bar_).
    dual_feasibility_vector = rhs;
    computeDualConeProjection(constraint_models, rhs, rhs);
    dual_feasibility_vector -= rhs;
    dual_feasibility = dual_feasibility_vector.template lpNorm<Eigen::Infinity>();
    const Scalar absolute_residual_warm_start = math::max(complementarity, dual_feasibility);

    dual_feasibility_vector.array() *= time_scaling_constraints_to_pos.array();
    dual_feasibility = dual_feasibility_vector.template lpNorm<Eigen::Infinity>();
    this->absolute_residual = math::max(complementarity, dual_feasibility);

    // Checking if the initial guess is better than 0.
    // if instead of the x_ initial guess, x_ is set to 0, then z_ = g.
    // -> we check how much constraints violation is induced by using g as the dual variable.
    // note: here we work with g and not gs, because we check for constraints violation at the
    // formulation level of each constraints.
    const Scalar absolute_residual_zero_guess =
      computeZeroInitialGuessMaxConstraintViolation(constraint_models, g);

    if (absolute_residual_zero_guess < absolute_residual_warm_start)
    { // If true, this means that the zero value initial guess leads a better feasibility in the
      // sense of the constraints satisfaction.
      // So we set the primal variables to the 0 initial guess and the dual variable to g.
      x_.setZero();
      y_.setZero();
      z_ = gs;
      if (solve_ncp)
      {
        {
          PINOCCHIO_TRACY_ZONE_SCOPED_N(
            "ADMMContactSolverTpl::solve - second computeDeSaxeCorrection");
          computeDeSaxeCorrection(constraint_models, z_, s_);
        }
        z_ += s_; // Add De Saxé shift
      }
      rhs = z_.array() * time_scaling_acc_to_constraints.array();
      dual_feasibility_vector = rhs;
      computeDualConeProjection(constraint_models, rhs, rhs);
      dual_feasibility_vector -= rhs; // Dual feasibility vector for the new null guess
      dual_feasibility_vector.array() *= time_scaling_constraints_to_pos.array();
      // We set the new convergence criterion
      this->absolute_residual = absolute_residual_zero_guess;
    }
    // We test convergence
    bool abs_prec_reached = this->absolute_residual < this->absolute_precision;

    if (!abs_prec_reached)
    { // the initial guess is not solution of the problem so we run the ADMM algorithm
      // Applying the preconditioner to work on a problem with a better scaling
      DelassusOperatorPreconditioned G_bar(_delassus, preconditioner_);
      rhs = VectorXs::Constant(this->problem_size, mu_prox);
      G_bar.updateDamping(rhs);      // G_bar =  P*(G+R)*P + mu_prox*Id
      scaleDualSolution(gs, g_bar_); // g_bar = P * gs
      scalePrimalSolution(x_, x_bar_);
      scalePrimalSolution(y_, y_bar_);
      scaleDualSolution(z_, z_bar_);

      // Setup ADMM update rules:
      // Before running ADMM, we compute the largest and smallest eigenvalues of delassus in order
      // to be able to use a spectral update rule for the proximal parameter (rho)
      // TODO should we evaluate the eigenvalues of G or Gbar ?
      Scalar L, m, rho;
      ADMMUpdateRuleContainer admm_update_rule_container;
      switch (admm_update_rule)
      {
      case (ADMMUpdateRule::SPECTRAL): {
        if (this->problem_size > 1)
        {
          PINOCCHIO_TRACY_ZONE_SCOPED_N("ADMMContactSolverTpl::solve - lanczos");
          m = rhs.minCoeff();
          this->lanczos_decomposition.compute(G_bar);
          L = ::pinocchio::computeLargestEigenvalue(this->lanczos_decomposition.Ts(), 1e-8);
#ifndef NDEBUG
          const bool enforce_symmetry = true;
          Eigen::MatrixXd delassus = G_bar.matrix(enforce_symmetry);
          Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(delassus);
          Eigen::VectorXd eigvals = solver.eigenvalues();
          // Scalar true_m = eigvals.minCoeff();
          Scalar true_L = eigvals.maxCoeff();
          //          if (true_m > 0)
          //          {
          //            assert(
          //              math::fabs((true_m - m) / math::max(true_m, m)) < 0.01
          //              && "true_m and m are too far apart.");
          //          }
          assert(
            math::fabs((true_L - L) / math::max(true_L, L)) < 0.01
            && "true_L and L are too far apart.");
#endif // NDEBUG
        }
        else
        {
          // TODO adapt this with Gbar
          typedef Eigen::Matrix<Scalar, 1, 1> Vector1;
          const Vector1 G = delassus * Vector1::Constant(1);
          m = L = G.coeff(0);
        }
        admm_update_rule_container.spectral_rule =
          ADMMSpectralUpdateRule(ratio_primal_dual, L, m, rho_power_factor);
        rho = ADMMSpectralUpdateRule::computeRho(L, m, rho_power);
        break;
      }
      case (ADMMUpdateRule::LINEAR):
        admm_update_rule_container.linear_rule =
          ADMMLinearUpdateRule(ratio_primal_dual, linear_update_rule_factor);
        rho = this->rho; // use the rho value stored in the solver.
        break;
      case (ADMMUpdateRule::CONSTANT):
        rho = this->rho; // use the rho value stored in the solver.
        break;
      }

      // clamp the rho
      rho = math::max(1e-8, rho);

      PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();

      // Update the cholesky decomposition
      Scalar prox_value = mu_prox + tau * rho;
      rhs = VectorXs::Constant(this->problem_size, prox_value);
      G_bar.updateDamping(rhs);
      Scalar old_prox_value = prox_value;
      cholesky_update_count = 1;

      is_initialized = true;

      // End of Initialization phase
      abs_prec_reached = false;
      bool rel_prec_reached = false;

      Scalar x_bar_norm_inf = x_bar_.template lpNorm<Eigen::Infinity>();
      Scalar y_bar_norm_inf = y_bar_.template lpNorm<Eigen::Infinity>();
      Scalar z_bar_norm_inf = z_bar_.template lpNorm<Eigen::Infinity>();
      Scalar x_bar_previous_norm_inf = x_bar_norm_inf;
      Scalar y_bar_previous_norm_inf = y_bar_norm_inf;
      Scalar z_bar_previous_norm_inf = z_bar_norm_inf;
      it = 1;
#ifdef PINOCCHIO_WITH_HPP_FCL
      timer.start();
#endif // PINOCCHIO_WITH_HPP_FCL
      for (; it <= Base::max_it; ++it)
      {

        x_bar_previous = x_bar_;
        y_bar_previous = y_bar_;
        z_bar_previous = z_bar_;
        complementarity = Scalar(0);

        if (solve_ncp)
        {
          // s-update
          computeDeSaxeCorrection(constraint_models, z_bar_, s_bar_);
        }

        // x-update
        rhs = -(g_bar_ + s_bar_ - (rho * tau) * y_bar_ - mu_prox * x_bar_ - z_bar_);
        G_bar.solveInPlace(rhs);
        x_bar_ = rhs;

        // y-update
        rhs -= z_bar_ / (tau * rho);
        computeScaledConeProjection(constraint_models, rhs, preconditioner_.getDiagonal(), y_bar_);

        // z-update
        z_bar_ -= (tau * rho) * (x_bar_ - y_bar_);

        // check termination criteria
        primal_feasibility_vector_bar = x_bar_ - y_bar_;

        {
          VectorXs & dx_bar = rhs;
          dx_bar = x_bar_ - x_bar_previous;
          dx_bar_norm =
            dx_bar.template lpNorm<Eigen::Infinity>(); // check relative progress on x_bar
          dual_feasibility_vector_bar = mu_prox * dx_bar;
        }

        {
          VectorXs & dy_bar = rhs;
          dy_bar = y_bar_ - y_bar_previous;
          dy_bar_norm =
            dy_bar.template lpNorm<Eigen::Infinity>(); // check relative progress on y_bar
          dual_feasibility_vector_bar += (tau * rho) * dy_bar;
        }

        {
          VectorXs & dz_bar = rhs;
          dz_bar = z_bar_ - z_bar_previous;
          dz_bar_norm =
            dz_bar.template lpNorm<Eigen::Infinity>(); // check relative progress on z_bar
        }

        // We unscale the quantities to work with stopping criterion from the original (unscaled)
        // problem
        unscalePrimalSolution(primal_feasibility_vector_bar, primal_feasibility_vector);
        primal_feasibility = primal_feasibility_vector.template lpNorm<Eigen::Infinity>();
        unscaleDualSolution(dual_feasibility_vector_bar, dual_feasibility_vector);
        const Scalar dual_feasibility_admm =
          dual_feasibility_vector.template lpNorm<Eigen::Infinity>();
        dual_feasibility_vector.array() *= time_scaling_acc_to_constraints.array();
        const Scalar dual_feasibility_constraint =
          dual_feasibility_vector.template lpNorm<Eigen::Infinity>();
        dual_feasibility_vector.array() *= time_scaling_constraints_to_pos.array();
        dual_feasibility = dual_feasibility_vector.template lpNorm<Eigen::Infinity>();
        unscalePrimalSolution(y_bar_, y_);
        unscaleDualSolution(z_bar_, z_);
        rhs = z_.array() * time_scaling_acc_to_constraints.array();
        complementarity = computeConicComplementarity(constraint_models, rhs, y_);

        if (stat_record)
        {
          VectorXs tmp(rhs);
          G_bar.applyOnTheRight(y_bar_, rhs);
          rhs += g_bar_ - prox_value * y_bar_;
          unscaleDualSolution(rhs, tmp);
          if (solve_ncp)
          {
            computeDeSaxeCorrection(constraint_models, tmp, rhs);
            tmp += rhs;
          }

          tmp.array() *=
            time_scaling_acc_to_constraints.array(); // back to constraint formulation level
          rhs = tmp;
          internal::computeDualConeProjection(constraint_models, rhs, rhs);
          tmp -= rhs;

          Scalar dual_feasibility_ncp = tmp.template lpNorm<Eigen::Infinity>();

          stats.primal_feasibility.push_back(primal_feasibility);
          stats.dual_feasibility.push_back(dual_feasibility);
          stats.dual_feasibility_admm.push_back(dual_feasibility_admm);
          stats.dual_feasibility_ncp.push_back(dual_feasibility_ncp);
          stats.dual_feasibility_constraint.push_back(dual_feasibility_constraint);
          stats.complementarity.push_back(complementarity);
          stats.rho.push_back(rho);
        }

        // Checking stopping residual
        const Scalar x_norm_inf = x_.template lpNorm<Eigen::Infinity>();
        const Scalar y_norm_inf = y_.template lpNorm<Eigen::Infinity>();
        const Scalar z_norm_inf = z_.template lpNorm<Eigen::Infinity>();
        if (
          check_expression_if_real<Scalar, false>(complementarity <= this->absolute_precision)
          && check_expression_if_real<Scalar, false>(
            dual_feasibility
            <= this->absolute_precision
                 + this->relative_precision * math::max(g_pos_norm_inf, z_norm_inf))
          && check_expression_if_real<Scalar, false>(
            primal_feasibility <= this->absolute_precision
                                    + this->relative_precision * math::max(x_norm_inf, y_norm_inf)))
          abs_prec_reached = true;
        else
          abs_prec_reached = false;

        x_bar_norm_inf = x_bar_.template lpNorm<Eigen::Infinity>();
        y_bar_norm_inf = y_bar_.template lpNorm<Eigen::Infinity>();
        z_bar_norm_inf = z_bar_.template lpNorm<Eigen::Infinity>();
        if (
          check_expression_if_real<Scalar, false>(
            dx_bar_norm
            <= this->relative_precision * math::max(x_bar_norm_inf, x_bar_previous_norm_inf))
          && check_expression_if_real<Scalar, false>(
            dy_bar_norm
            <= this->relative_precision * math::max(y_bar_norm_inf, y_bar_previous_norm_inf))
          && check_expression_if_real<Scalar, false>(
            dz_bar_norm
            <= this->relative_precision * math::max(z_bar_norm_inf, z_bar_previous_norm_inf)))
          rel_prec_reached = true;
        else
          rel_prec_reached = false;

        if (abs_prec_reached || rel_prec_reached)
          break;

        // Apply rho according to the primal_dual_ratio
        bool update_delassus_factorization = false;
        switch (admm_update_rule)
        {
        case (ADMMUpdateRule::SPECTRAL):
          update_delassus_factorization = admm_update_rule_container.spectral_rule.eval(
            primal_feasibility, dual_feasibility_constraint, rho);
          break;
        case (ADMMUpdateRule::LINEAR):
          update_delassus_factorization = admm_update_rule_container.linear_rule.eval(
            primal_feasibility, dual_feasibility_constraint, rho);
          break;
        case (ADMMUpdateRule::CONSTANT):
          break;
        }

        // clamp rho
        rho = math::max(1e-8, rho);

        // Account for potential update of rho
        if (update_delassus_factorization)
        {
          prox_value = mu_prox + tau * rho;
          if (old_prox_value != prox_value)
          {
            rhs = VectorXs::Constant(this->problem_size, prox_value);
            G_bar.updateDamping(rhs);
            cholesky_update_count++;
            old_prox_value = prox_value;
          }
        }

        x_bar_previous_norm_inf = x_bar_norm_inf;
        y_bar_previous_norm_inf = y_bar_norm_inf;
        z_bar_previous_norm_inf = z_bar_norm_inf;
      } // end ADMM main for loop

      unscalePrimalSolution(x_bar_, x_);
      unscalePrimalSolution(y_bar_, y_);
      unscaleDualSolution(z_bar_, z_);
      unscaleDualSolution(s_bar_, s_);

      this->relative_residual = math::max(
        dx_bar_norm / math::max(x_bar_norm_inf, x_bar_previous_norm_inf),
        dy_bar_norm / math::max(y_bar_norm_inf, y_bar_previous_norm_inf));
      this->relative_residual = math::max(
        this->relative_residual, dz_bar_norm / math::max(z_bar_norm_inf, z_bar_previous_norm_inf));
      this->absolute_residual =
        math::max(primal_feasibility, math::max(complementarity, dual_feasibility));

      // Save values of spectral update rule
      if (admm_update_rule == ADMMUpdateRule::SPECTRAL)
      {
        this->rho_power = ADMMSpectralUpdateRule::computeRhoPower(L, m, rho);
        this->rho = rho;
      }
    }
    PINOCCHIO_EIGEN_MALLOC_ALLOWED();

#ifdef PINOCCHIO_WITH_HPP_FCL
    timer.stop();
#endif // PINOCCHIO_WITH_HPP_FCL
    //

    this->it = it;
    // we time-rescale dual solution and desaxe correction
    // so that z_ and s_ are back at the constraints formulations levels
    z_constraint_ = z_.array() * time_scaling_acc_to_constraints.array();
    s_constraint_ = s_.array() * time_scaling_acc_to_constraints.array();

    if (stat_record)
    {
      stats.it = it;
      stats.cholesky_update_count = cholesky_update_count;
    }

    if (abs_prec_reached)
      return true;

    return false;
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_admm_solver_hxx__
