//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_algorithm_admm_solver_hxx__
#define __pinocchio_algorithm_admm_solver_hxx__

#include <limits>
#include <iomanip>
#include <iostream>

#include "pinocchio/algorithm/contact-solver-utils.hpp"
#include "pinocchio/algorithm/constraints/coulomb-friction-cone.hpp"
#include "pinocchio/algorithm/constraints/visitors/constraint-model-visitor.hpp"

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

    template<typename ConstraintSet, typename VectorLike>
    static void algo_impl(
      const ConstraintSet & set,
      const Eigen::MatrixBase<VectorLike> & drift,
      Scalar & max_violation)
    {
      // TODO: for other sets, how should we compute?
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
    typename ConstraintModelAllocator,
    typename VectorLikeR>
  bool ADMMContactSolverTpl<_Scalar>::solve(
    DelassusOperatorBase<DelassusDerived> & _delassus,
    const Eigen::MatrixBase<VectorLike> & g,
    const std::vector<Holder<const ConstraintModel>, ConstraintModelAllocator> & constraint_models,
    const Eigen::MatrixBase<VectorLikeR> & R,
    const boost::optional<ConstRefVectorXs> primal_guess,
    const boost::optional<ConstRefVectorXs> dual_guess,
    bool solve_ncp,
    bool compute_largest_eigen_values,
    ADMMUpdateRule admm_update_rule,
    bool stat_record)

  {
    using namespace internal;
    typedef ADMMSpectralUpdateRuleTpl<Scalar> ADMMSpectralUpdateRule;
    typedef ADMMLinearUpdateRuleTpl<Scalar> ADMMLinearUpdateRule;

    typedef ADMMUpdateRuleContainerTpl<Scalar> ADMMUpdateRuleContainer;
    DelassusDerived & delassus = _delassus.derived();

    const Scalar mu_R = R.minCoeff();
    PINOCCHIO_CHECK_INPUT_ARGUMENT(tau <= Scalar(1) && tau > Scalar(0), "tau should lie in ]0,1].");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(mu_prox >= 0, "mu_prox should be positive.");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(mu_R >= Scalar(0), "R should be a positive vector.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(R.size(), problem_size);

    // First, we initialize the primal and dual variables
    int it = 0;
    Scalar complementarity, dx_norm, dy_norm, dz_norm, //
      primal_feasibility, dual_feasibility_ncp, dual_feasibility;

    // we add the compliance to the delassus
    rhs = R + VectorXs::Constant(this->problem_size, mu_prox);
    delassus.updateDamping(rhs);
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
    else if (!is_initialized)
    {
      x_.setZero();
    }
    else
    {
      x_ = y_; // takes the current value stored in the solver
    }

    // Init y
    computeConeProjection(constraint_models, x_, y_);

    // Init z
    if (dual_guess)
    {
      z_ = dual_guess.get();
      PINOCCHIO_CHECK_ARGUMENT_SIZE(z_.size(), problem_size);
    }
    else if (!is_initialized)
    {
      delassus.applyOnTheRight(y_, z_); // z = (G + R + mu_prox*Id)* y
      z_.noalias() += -mu_prox * y_ + g;
      if (solve_ncp)
      {
        computeComplementarityShift(constraint_models, z_, s_);
        z_ += s_; // Add De Saxé shift
      }
    }

    dual_feasibility_vector = z_;
    computeDualConeProjection(constraint_models, z_, z_);
    dual_feasibility_vector -= z_;

    // Checking if the initial guess is better than 0
    complementarity = computeConicComplementarity(
      constraint_models, z_, y_); // Complementarity of the initial guess
    // we Search for the max violation of the constraints
    // Note: max_violation is always <= 0
    const Scalar max_violation =
      computeZeroInitialGuessMaxConstraintViolation(constraint_models, g);

    if (max_violation < complementarity)
    { // If true, this means that the zero value initial guess leads a better feasibility in the
      // sense of the contact complementarity
      x_.setZero();
      y_.setZero();
      z_ = g;
      if (solve_ncp)
      {
        computeComplementarityShift(constraint_models, z_, s_);
        z_ += s_; // Add De Saxé shift
      }
      dual_feasibility_vector = z_;
      computeDualConeProjection(constraint_models, z_, z_);
      dual_feasibility_vector -= z_; // Dual feasibility vector for the new null guess
      complementarity = computeConicComplementarity(
        constraint_models, z_, y_); // Complementarity of the new null guess
    }
    // We compute the convergence criterion
    dual_feasibility = dual_feasibility_vector.template lpNorm<Eigen::Infinity>();
    this->absolute_residual = math::max(complementarity, dual_feasibility);
    bool abs_prec_reached = this->absolute_residual < this->absolute_precision;

    if (!abs_prec_reached)
    { // the initial guess is not solution of the problem so we run the ADMM algorithm
      // Before running ADMM, we compute the largest and smallest eigenvalues of delassus in order
      // to be able to use a spectral update rule for the proximal parameter (rho) Setup ADMM update
      // rules
      Scalar L, m, rho;
      ADMMUpdateRuleContainer admm_update_rule_container;
      switch (admm_update_rule)
      {
      case (ADMMUpdateRule::SPECTRAL):
        if (compute_largest_eigen_values)
        {
          power_iteration_algo.run(delassus);
        }
        m = mu_prox + mu_R;
        L = power_iteration_algo.largest_eigen_value;
        admm_update_rule_container.spectral_rule =
          ADMMSpectralUpdateRule(ratio_primal_dual, L, m, rho_power_factor);
        rho = ADMMSpectralUpdateRule::computeRho(L, m, rho_power);
        break;
      case (ADMMUpdateRule::LINEAR):
        admm_update_rule_container.linear_rule =
          ADMMLinearUpdateRule(ratio_primal_dual, linear_update_rule_factor);
        rho = this->rho; // use the rho value stored in the solver.
        break;
      }

      // clamp the rho
      rho = math::max(1e-8, rho);

      PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();

      // Update the cholesky decomposition
      Scalar prox_value = mu_prox + tau * rho;
      rhs = R + VectorXs::Constant(this->problem_size, prox_value);
      delassus.updateDamping(rhs);
      cholesky_update_count = 1;

      if (stat_record)
      {
        stats.reset();

        // Compute initial problem primal and dual feasibility
        primal_feasibility_vector = x_ - y_;
        primal_feasibility = primal_feasibility_vector.template lpNorm<Eigen::Infinity>();
      }

      is_initialized = true;

      // End of Initialization phase
      abs_prec_reached = false;
      bool rel_prec_reached = false;

      Scalar x_norm_inf = x_.template lpNorm<Eigen::Infinity>();
      Scalar y_norm_inf = y_.template lpNorm<Eigen::Infinity>();
      Scalar z_norm_inf = z_.template lpNorm<Eigen::Infinity>();
      Scalar x_previous_norm_inf = x_norm_inf;
      Scalar y_previous_norm_inf = y_norm_inf;
      Scalar z_previous_norm_inf = z_norm_inf;
      it = 1;
#ifdef PINOCCHIO_WITH_HPP_FCL
      timer.start();
#endif // PINOCCHIO_WITH_HPP_FCL
      for (; it <= Base::max_it; ++it)
      {

        x_previous = x_;
        y_previous = y_;
        z_previous = z_;
        complementarity = Scalar(0);

        if (solve_ncp)
        {
          // s-update
          computeComplementarityShift(constraint_models, z_, s_);
        }

        // x-update
        rhs = -(g + s_ - (rho * tau) * y_ - mu_prox * x_ - z_);
        delassus.solveInPlace(rhs);
        x_ = rhs;

        // y-update
        rhs -= z_ / (tau * rho);
        computeConeProjection(constraint_models, rhs, y_);

        // z-update
        z_ -= (tau * rho) * (x_ - y_);

        // check termination criteria
        primal_feasibility_vector = x_ - y_;

        {
          VectorXs & dx = rhs;
          dx = x_ - x_previous;
          dx_norm = dx.template lpNorm<Eigen::Infinity>(); // check relative progress on x
          dual_feasibility_vector.noalias() += mu_prox * dx;
        }

        {
          VectorXs & dy = rhs;
          dy = y_ - y_previous;
          dy_norm = dy.template lpNorm<Eigen::Infinity>(); // check relative progress on y
          dual_feasibility_vector.noalias() = (tau * rho) * dy;
        }

        {
          VectorXs & dz = rhs;
          dz = z_ - z_previous;
          dz_norm = dz.template lpNorm<Eigen::Infinity>(); // check relative progress on z
        }

        primal_feasibility = primal_feasibility_vector.template lpNorm<Eigen::Infinity>();
        dual_feasibility = dual_feasibility_vector.template lpNorm<Eigen::Infinity>();
        complementarity = computeConicComplementarity(constraint_models, z_, y_);
        //      complementarity = z_.dot(y_)/constraint_models.size();

        if (stat_record)
        {
          VectorXs tmp(rhs);
          delassus.applyOnTheRight(y_, rhs);
          rhs.noalias() += g - prox_value * y_;
          if (solve_ncp)
          {
            computeComplementarityShift(constraint_models, rhs, tmp);
            rhs.noalias() += tmp;
          }

          internal::computeDualConeProjection(constraint_models, rhs, tmp);
          tmp -= rhs;

          dual_feasibility_ncp = tmp.template lpNorm<Eigen::Infinity>();

          stats.primal_feasibility.push_back(primal_feasibility);
          stats.dual_feasibility.push_back(dual_feasibility);
          stats.dual_feasibility_ncp.push_back(dual_feasibility_ncp);
          stats.complementarity.push_back(complementarity);
          stats.rho.push_back(rho);
        }

        // Checking stopping residual
        if (
          check_expression_if_real<Scalar, false>(complementarity <= this->absolute_precision)
          && check_expression_if_real<Scalar, false>(dual_feasibility <= this->absolute_precision)
          && check_expression_if_real<Scalar, false>(
            primal_feasibility <= this->absolute_precision))
          abs_prec_reached = true;
        else
          abs_prec_reached = false;

        x_norm_inf = x_.template lpNorm<Eigen::Infinity>();
        y_norm_inf = y_.template lpNorm<Eigen::Infinity>();
        z_norm_inf = z_.template lpNorm<Eigen::Infinity>();
        if (
          check_expression_if_real<Scalar, false>(
            dx_norm <= this->relative_precision * math::max(x_norm_inf, x_previous_norm_inf))
          && check_expression_if_real<Scalar, false>(
            dy_norm <= this->relative_precision * math::max(y_norm_inf, y_previous_norm_inf))
          && check_expression_if_real<Scalar, false>(
            dz_norm <= this->relative_precision * math::max(z_norm_inf, z_previous_norm_inf)))
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
            primal_feasibility, dual_feasibility, rho);
          break;
        case (ADMMUpdateRule::LINEAR):
          update_delassus_factorization =
            admm_update_rule_container.linear_rule.eval(primal_feasibility, dual_feasibility, rho);
          ;
          break;
        }

        // clamp rho
        rho = math::max(1e-8, rho);

        // Account for potential update of rho
        if (update_delassus_factorization)
        {
          prox_value = mu_prox + tau * rho;
          rhs = R + VectorXs::Constant(this->problem_size, prox_value);
          delassus.updateDamping(rhs);
          cholesky_update_count++;
        }

        x_previous_norm_inf = x_norm_inf;
        y_previous_norm_inf = y_norm_inf;
        z_previous_norm_inf = z_norm_inf;
      } // end ADMM main for loop
      this->relative_residual = math::max(
        dx_norm / math::max(x_norm_inf, x_previous_norm_inf),
        dy_norm / math::max(y_norm_inf, y_previous_norm_inf));
      this->relative_residual =
        math::max(this->relative_residual, dz_norm / math::max(z_norm_inf, z_previous_norm_inf));
      this->absolute_residual =
        math::max(primal_feasibility, math::max(complementarity, dual_feasibility));

      // Save values
      this->rho_power = ADMMSpectralUpdateRule::computeRhoPower(L, m, rho);
      this->rho = rho;

      if (stat_record)
      {
        stats.it = it;
        stats.cholesky_update_count = cholesky_update_count;
      }
    }
    PINOCCHIO_EIGEN_MALLOC_ALLOWED();

#ifdef PINOCCHIO_WITH_HPP_FCL
    timer.stop();
#endif // PINOCCHIO_WITH_HPP_FCL
    //

    this->it = it;

    //    if(abs_prec_reached || rel_prec_reached)
    if (abs_prec_reached)
      return true;

    return false;
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_admm_solver_hxx__
