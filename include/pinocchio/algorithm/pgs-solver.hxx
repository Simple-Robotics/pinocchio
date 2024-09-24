//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_algorithm_pgs_solver_hxx__
#define __pinocchio_algorithm_pgs_solver_hxx__

#include <limits>
#include "pinocchio/algorithm/constraints/coulomb-friction-cone.hpp"

namespace pinocchio
{

  template<typename _ConstraintSet>
  struct PGSConstraintProjectionStep
  {
    typedef _ConstraintSet ConstraintSet;
    typedef typename ConstraintSet::Scalar Scalar;
    typedef Eigen::Matrix<Scalar, 3, 1> Vector3;

    PGSConstraintProjectionStep(const Scalar over_relax_value, const ConstraintSet & set)
    : over_relax_value(over_relax_value)
    , set(set)
    , complementarity(Scalar(-1))
    , dual_feasibility(Scalar(-1))
    , primal_feasibility(Scalar(0))
    {
    }

    ///
    /// \brief Perform a projection step associated with the PGS algorithm
    ///
    /// \param[in] G_block block asscociated with the current
    /// \param[in,out] primal_vector_ primal vector which will be update with the new estimate
    /// \param[in,out] dual_vector_ primal vector which will be update with the new estimate
    ///
    template<typename BlockType, typename PrimalVectorType, typename DualVectorType>
    void project(
      const Eigen::MatrixBase<BlockType> & G_block,
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector_,
      const Eigen::MatrixBase<DualVectorType> & dual_vector_) const
    {
      typedef Eigen::Matrix<Scalar, 2, 1> Vector2;

      auto & primal_vector = primal_vector_.const_cast_derived();
      auto & dual_vector = dual_vector_.const_cast_derived();

      // Normal update
      Scalar & fz = dual_vector.coeffRef(2);
      const Scalar fz_previous = fz;
      fz -= Scalar(this->over_relax_value / G_block.coeff(2, 2)) * primal_vector[2];
      fz = math::max(Scalar(0), fz);

      // Account for the fz updated value
      primal_vector += G_block.col(2) * (fz - fz_previous);

      // Tangential update
      const Scalar min_D_tangent = math::min(G_block.coeff(0, 0), G_block.coeff(1, 1));
      auto f_tangent = dual_vector.template head<2>();
      const Vector2 f_tangent_previous = f_tangent;

      assert(min_D_tangent > 0 && "min_D_tangent is zero");
      f_tangent -= this->over_relax_value / min_D_tangent * primal_vector.template head<2>();
      const Scalar f_tangent_norm = f_tangent.norm();

      const Scalar mu_fz = this->set.mu * fz;
      if (f_tangent_norm > mu_fz) // Project in the circle of radius mu_fz
        f_tangent *= mu_fz / f_tangent_norm;

      // Account for the f_tangent updated value
      primal_vector.noalias() += G_block.template leftCols<2>() * (f_tangent - f_tangent_previous);
    }

    /// \brief Compute the feasibility conditions associated with the optimization problem
    template<typename PrimalVectorType, typename DualVectorType>
    void computeFeasibility(
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector,
      const Eigen::MatrixBase<DualVectorType> & dual_vector)
    {
      // TODO(jcarpent): change primal_feasibility and dual_feasibility
      this->primal_feasibility =
        Scalar(0); // always zero as the dual variable belongs to the friction cone.
      this->complementarity = this->set.computeContactComplementarity(primal_vector, dual_vector);
      assert(this->complementarity >= Scalar(0) && "The complementarity should be positive");

      typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
      const Vector3 primal_vector_corrected =
        primal_vector + this->set.computeNormalCorrection(primal_vector);
      const Vector3 reprojection_residual =
        this->set.dual().project(primal_vector_corrected) - primal_vector_corrected;
      this->dual_feasibility = reprojection_residual.norm();
    }

    const Scalar over_relax_value;
    const ConstraintSet & set;
    Scalar complementarity;
    Scalar dual_feasibility;
    Scalar primal_feasibility;

  }; // PGSConstraintProjectionStep

  template<typename _Scalar>
  template<
    typename MatrixLike,
    typename VectorLike,
    typename ConstraintAllocator,
    typename VectorLikeOut>
  bool PGSContactSolverTpl<_Scalar>::solve(
    const MatrixLike & G,
    const Eigen::MatrixBase<VectorLike> & g,
    const std::vector<CoulombFrictionConeTpl<Scalar>, ConstraintAllocator> & cones,
    const Eigen::DenseBase<VectorLikeOut> & x_sol,
    const Scalar over_relax)

  {
    typedef CoulombFrictionConeTpl<Scalar> CoulombFrictionCone;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorX;

    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      over_relax < Scalar(2) && over_relax > Scalar(0), "over_relax should lie in ]0,2[.")
    PINOCCHIO_CHECK_ARGUMENT_SIZE(g.size(), this->getProblemSize());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(G.rows(), this->getProblemSize());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(G.cols(), this->getProblemSize());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(x_sol.size(), this->getProblemSize());

    const size_t nc = cones.size(); // num constraints

    int it = 0;
    PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();

#ifdef PINOCCHIO_WITH_HPP_FCL
    timer.start();
#endif // PINOCCHIO_WITH_HPP_FCL

    Scalar complementarity, proximal_metric, primal_feasibility, dual_feasibility;
    bool abs_prec_reached = false, rel_prec_reached = false;
    x = x_sol;
    Scalar x_previous_norm_inf = x.template lpNorm<Eigen::Infinity>();

    const Eigen::DenseIndex constraint_set_size_max = 3;
    VectorX velocity_storage(constraint_set_size_max); // tmp variable
    for (; it < this->max_it; ++it)
    {
      x_previous = x;
      complementarity = Scalar(0);
      dual_feasibility = Scalar(0);
      primal_feasibility = Scalar(0);
      for (size_t cone_id = 0; cone_id < nc; ++cone_id)
      {
        const CoulombFrictionCone & cone = cones[cone_id];
        const int constraint_set_size = cone.size();
        const Eigen::DenseIndex row_id = constraint_set_size * Eigen::DenseIndex(cone_id);

        const auto G_block = G.block(row_id, row_id, constraint_set_size, constraint_set_size);
        auto force = x.segment(row_id, constraint_set_size);

        auto velocity = velocity_storage.head(constraint_set_size);

        // Update primal variable
        velocity.noalias() =
          G.middleRows(row_id, constraint_set_size) * x + g.segment(row_id, constraint_set_size);

        PGSConstraintProjectionStep<CoulombFrictionCone> step(over_relax, cone);
        step.project(G_block, velocity, force);
        step.computeFeasibility(velocity, force);

        // Update problem feasibility
        complementarity = math::max(complementarity, step.complementarity);
        dual_feasibility = math::max(dual_feasibility, step.dual_feasibility);
        primal_feasibility = math::max(primal_feasibility, step.primal_feasibility);
      }

      // Checking stopping residual
      if (
        check_expression_if_real<Scalar, false>(complementarity <= this->absolute_precision)
        && check_expression_if_real<Scalar, false>(dual_feasibility <= this->absolute_precision))
        abs_prec_reached = true;
      else
        abs_prec_reached = false;

      proximal_metric = (x - x_previous).template lpNorm<Eigen::Infinity>();
      const Scalar x_norm_inf = x.template lpNorm<Eigen::Infinity>();
      if (check_expression_if_real<Scalar, false>(
            proximal_metric
            <= this->relative_precision * math::max(x_norm_inf, x_previous_norm_inf)))
        rel_prec_reached = true;
      else
        rel_prec_reached = false;

      if (abs_prec_reached || rel_prec_reached)
        break;

      x_previous_norm_inf = x_norm_inf;
    }

#ifdef PINOCCHIO_WITH_HPP_FCL
    timer.stop();
#endif // PINOCCHIO_WITH_HPP_FCL

    PINOCCHIO_EIGEN_MALLOC_ALLOWED();

    this->absolute_residual = math::max(complementarity, dual_feasibility);
    this->relative_residual = proximal_metric;
    this->it = it;
    x_sol.const_cast_derived() = x;

    if (abs_prec_reached || rel_prec_reached)
      return true;

    return false;
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_pgs_solver_hxx__
