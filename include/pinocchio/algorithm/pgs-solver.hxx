//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_algorithm_pgs_solver_hxx__
#define __pinocchio_algorithm_pgs_solver_hxx__

#include "pinocchio/algorithm/constraints/sets.hpp"
#include "pinocchio/algorithm/constraints/visitors/constraint-model-visitor.hpp"
#include "pinocchio/algorithm/contact-solver-utils.hpp"

namespace pinocchio
{

  template<typename _Scalar>
  struct PGSConstraintProjectionStepBase
  {
    typedef _Scalar Scalar;

    explicit PGSConstraintProjectionStepBase(const Scalar over_relax_value)
    : over_relax_value(over_relax_value)
    {
    }

    const Scalar over_relax_value;
    Scalar complementarity;
    Scalar dual_feasibility;
    Scalar primal_feasibility;
  }; // PGSConstraintProjectionBase

  template<typename ConstraintSet>
  struct PGSConstraintProjectionStep
  {
  };

  template<
    typename Scalar,
    typename BlockType,
    typename ForceType,
    typename VelocityType,
    typename DualToPosType,
    typename TimeScalingType1,
    typename TimeScalingType2>
  struct PGSConstraintProjectionStepVisitor
  : visitors::ConstraintUnaryVisitorBase<PGSConstraintProjectionStepVisitor<
      Scalar,
      BlockType,
      ForceType,
      VelocityType,
      DualToPosType,
      TimeScalingType1,
      TimeScalingType2>>
  , PGSConstraintProjectionStepBase<Scalar>
  {
    typedef boost::fusion::vector<
      const Scalar,
      const BlockType &,
      ForceType &,
      VelocityType &,
      DualToPosType &,
      const TimeScalingType1 &,
      const TimeScalingType2 &,
      Scalar &,
      Scalar &,
      Scalar &>
      ArgsType;
    typedef PGSConstraintProjectionStepBase<Scalar> Base;
    typedef visitors::ConstraintUnaryVisitorBase<PGSConstraintProjectionStepVisitor<
      Scalar,
      BlockType,
      ForceType,
      VelocityType,
      DualToPosType,
      TimeScalingType1,
      TimeScalingType2>>
      VisitorBase;

    explicit PGSConstraintProjectionStepVisitor(const Scalar over_relax_value)
    : Base(over_relax_value)
    {
    }

    template<typename ConstraintModel>
    static void algo(
      const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
      const Scalar over_relax_value,
      const Eigen::EigenBase<BlockType> & G_block,
      ForceType & force,
      VelocityType & velocity,
      DualToPosType & dual_to_pos,
      const TimeScalingType1 & time_scaling_acc_to_constraints,
      const TimeScalingType2 & time_scaling_constraints_to_pos,
      Scalar & complementarity,
      Scalar & primal_feasibility,
      Scalar & dual_feasibility)
    {
      typedef typename ConstraintModel::ConstraintSet ConstraintSet;

      PGSConstraintProjectionStep<ConstraintSet> step(
        over_relax_value,
        cmodel.derived().set()); // TODO(jcarpent): change cmodel.derived().set() -> cmodel.set()
      step.project(G_block.derived(), force.const_cast_derived(), velocity.const_cast_derived());
      dual_to_pos = velocity.array() * time_scaling_acc_to_constraints.array();
      dual_to_pos.array() *= time_scaling_constraints_to_pos.array();
      step.computeFeasibility(force, dual_to_pos);

      complementarity = step.complementarity;
      dual_feasibility = step.dual_feasibility;
      primal_feasibility = step.primal_feasibility;
    }

    using VisitorBase::run;
    template<typename ConstraintModel>
    void run(
      const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
      const Eigen::EigenBase<BlockType> & G_block,
      ForceType & force,
      VelocityType & velocity,
      DualToPosType & dual_to_pos,
      const TimeScalingType1 & time_scaling_acc_to_constraints,
      const TimeScalingType2 & time_scaling_constraints_to_pos)
    {
      algo(
        cmodel.derived(), this->over_relax_value, G_block.derived(), force, velocity, dual_to_pos,
        time_scaling_acc_to_constraints, time_scaling_constraints_to_pos, this->complementarity,
        this->primal_feasibility, this->dual_feasibility);
    }

    template<int Options, template<typename S, int O> class ConstraintCollectionTpl>
    void run(
      const pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const Eigen::EigenBase<BlockType> & G_block,
      ForceType & force,
      VelocityType & velocity,
      DualToPosType & dual_to_pos,
      const TimeScalingType1 & time_scaling_acc_to_constraints,
      const TimeScalingType2 & time_scaling_constraints_to_pos)
    {
      ArgsType args(
        this->over_relax_value, G_block.derived(), force, velocity, dual_to_pos,
        time_scaling_acc_to_constraints, time_scaling_constraints_to_pos, this->complementarity,
        this->primal_feasibility, this->dual_feasibility);
      this->run(cmodel.derived(), args);
    }
  }; // struct PGSConstraintProjectionStepVisitor

  template<typename _Scalar>
  struct PGSConstraintProjectionStep<CoulombFrictionConeTpl<_Scalar>>
  : PGSConstraintProjectionStepBase<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef CoulombFrictionConeTpl<Scalar> ConstraintSet;
    typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
    typedef PGSConstraintProjectionStepBase<Scalar> Base;

    PGSConstraintProjectionStep(const Scalar over_relax_value, const ConstraintSet & set)
    : Base(over_relax_value)
    , set(set)
    {
    }

    ///
    /// \brief Perform a projection step associated with the PGS algorithm
    ///
    /// \param[in] G_block block asscociated with the current
    /// \param[in,out] primal_vector_ primal vector which will be update with the new estimate
    /// \param[in,out] dual_vector_ dual vector which will be update with the new estimate
    ///
    template<typename BlockType, typename PrimalVectorType, typename DualVectorType>
    void project(
      const Eigen::EigenBase<BlockType> & G_block_,
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector_,
      const Eigen::MatrixBase<DualVectorType> & dual_vector_) const
    {
      typedef Eigen::Matrix<Scalar, 2, 1> Vector2;

      auto & G_block = G_block_.derived();
      auto & primal_vector = primal_vector_.const_cast_derived();
      auto & dual_vector = dual_vector_.const_cast_derived();

      // Normal update
      Scalar & fz = primal_vector.coeffRef(2);
      const Scalar fz_previous = fz;
      fz -= Scalar(
              this->over_relax_value
              / math::max(Eigen::NumTraits<Scalar>::dummy_precision(), G_block.coeff(2, 2)))
            * dual_vector[2];
      fz = math::max(Scalar(0), fz);

      // Account for the fz updated value
      dual_vector += G_block.col(2) * (fz - fz_previous);

      // Tangential update
      const Scalar min_D_tangent = math::max(
        Eigen::NumTraits<Scalar>::dummy_precision(),
        math::min(G_block.coeff(0, 0), G_block.coeff(1, 1)));
      auto f_tangent = primal_vector.template head<2>();
      const Vector2 f_tangent_previous = f_tangent;

      f_tangent -= this->over_relax_value / min_D_tangent * dual_vector.template head<2>();
      const Scalar f_tangent_norm = f_tangent.norm();

      const Scalar mu_fz = this->set.mu * fz;
      if (f_tangent_norm > mu_fz) // Project in the circle of radius mu_fz
      {
        assert(f_tangent_norm > 0 && "f_tangent_norm is zero");
        f_tangent *= mu_fz / f_tangent_norm;
      }

      // Account for the f_tangent updated value
      dual_vector.noalias() += G_block.template leftCols<2>() * (f_tangent - f_tangent_previous);
    }

    /// \brief Compute the feasibility conditions associated with the optimization problem
    template<typename PrimalVectorType, typename DualVectorType>
    void computeFeasibility(
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector,
      const Eigen::MatrixBase<DualVectorType> & dual_vector)
    {
      // The name should be inverted.
      this->primal_feasibility =
        Scalar(0); // always zero as the primal variable belongs to the friction cone.

      typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
      const Vector3 dual_vector_corrected =
        dual_vector + this->set.computeNormalCorrection(dual_vector);
      this->complementarity =
        this->set.computeConicComplementarity(dual_vector_corrected, primal_vector);
      assert(this->complementarity >= Scalar(0) && "The complementarity should be positive");
      const Vector3 reprojection_residual =
        this->set.dual().project(dual_vector_corrected) - dual_vector_corrected;
      this->dual_feasibility = reprojection_residual.norm();
    }

    const ConstraintSet & set;

  }; // PGSConstraintProjectionStep<CoulombFrictionConeTpl<_Scalar>>

  template<typename _Scalar, int _Options>
  struct PGSConstraintProjectionStep<UnboundedSetTpl<_Scalar, _Options>>
  : PGSConstraintProjectionStepBase<_Scalar>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef UnboundedSetTpl<Scalar, Options> ConstraintSet;
    typedef PGSConstraintProjectionStepBase<Scalar> Base;

    PGSConstraintProjectionStep(const Scalar over_relax_value, const ConstraintSet & set)
    : Base(over_relax_value)
    , set(set)
    {
    }

    ///
    /// \brief Perform a projection step associated with the PGS algorithm
    ///
    /// \param[in] G_block block asscociated with the current
    /// \param[in,out] primal_vector_ primal vector which will be update with the new estimate
    /// \param[in,out] dual_vector_ dual vector which will be update with the new estimate
    ///
    template<typename BlockType, typename PrimalVectorType, typename DualVectorType>
    void project(
      const Eigen::MatrixBase<BlockType> & G_block_,
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector_,
      const Eigen::MatrixBase<DualVectorType> & dual_vector_) const
    {

      const Eigen::DenseIndex size = set.size();
      auto & G_block = G_block_.derived();
      auto & primal_vector = primal_vector_.const_cast_derived();
      auto & dual_vector = dual_vector_.const_cast_derived();

      for (Eigen::DenseIndex i = 0; i < size; ++i)
      {
        Scalar d_primal_value =
          -this->over_relax_value * dual_vector[i]
          / math::max(Eigen::NumTraits<Scalar>::dummy_precision(), G_block.coeff(i, i));
        primal_vector[i] += d_primal_value;
        dual_vector.noalias() += G_block.col(i) * d_primal_value; // TODO: this could be optimized
      }
    }

    ///
    /// \brief Perform a projection step associated with the PGS algorithm
    ///
    /// \param[in] G_block block asscociated with the current
    /// \param[in,out] primal_vector_ primal vector which will be update with the new estimate
    /// \param[in,out] dual_vector_ dual vector which will be update with the new estimate
    ///
    template<typename BlockType, typename PrimalVectorType, typename DualVectorType>
    void project(
      const Eigen::EigenBase<BlockType> & G_block_,
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector_,
      const Eigen::MatrixBase<DualVectorType> & dual_vector_) const
    {

      const Eigen::DenseIndex size = set.size();
      auto & G_block = G_block_.derived();
      auto & primal_vector = primal_vector_.const_cast_derived();
      auto & dual_vector = dual_vector_.const_cast_derived();

      for (Eigen::DenseIndex i = 0; i < size; ++i)
      {
        Scalar d_primal_value =
          -this->over_relax_value * dual_vector[i]
          / math::max(Eigen::NumTraits<Scalar>::dummy_precision(), G_block.coeff(i, i));
        primal_vector[i] += d_primal_value;
        dual_vector += G_block.col(i) * d_primal_value; // TODO: this could be optimized using aloca
      }
    }

    /// \brief Compute the feasibility conditions associated with the optimization problem
    template<typename PrimalVectorType, typename DualVectorType>
    void computeFeasibility(
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector,
      const Eigen::MatrixBase<DualVectorType> & dual_vector)
    {
      this->primal_feasibility = Scalar(0);
      this->complementarity = primal_vector.dot(dual_vector);
      this->dual_feasibility = dual_vector.template lpNorm<Eigen::Infinity>();
    }

    const ConstraintSet & set;

  }; // PGSConstraintProjectionStep<UnboundedSetTpl<_Scalar,_Options>>

  template<typename _Scalar>
  struct PGSConstraintProjectionStep<BoxSetTpl<_Scalar>> : PGSConstraintProjectionStepBase<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef BoxSetTpl<Scalar> ConstraintSet;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef PGSConstraintProjectionStepBase<Scalar> Base;

    PGSConstraintProjectionStep(const Scalar over_relax_value, const ConstraintSet & set)
    : Base(over_relax_value)
    , set(set)
    {
    }

    ///
    /// \brief Perform a projection step associated with the PGS algorithm
    ///
    /// \param[in] G_block block asscociated with the current
    /// \param[in,out] primal_vector_ primal vector which will be update with the new estimate
    /// \param[in,out] dual_vector_ dual vector which will be update with the new estimate
    ///
    template<typename BlockType, typename PrimalVectorType, typename DualVectorType>
    void project(
      const Eigen::MatrixBase<BlockType> & G_block,
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector,
      const Eigen::MatrixBase<DualVectorType> & dual_vector) const
    {
      project_impl(
        this->set, this->over_relax_value, G_block.derived(), primal_vector.const_cast_derived(),
        dual_vector.const_cast_derived());
    }

    template<
      typename ConstraintSetType,
      typename BlockType,
      typename PrimalVectorType,
      typename DualVectorType>
    static void project_impl(
      const ConstraintSetType & set,
      const Scalar over_relax_value,
      const Eigen::MatrixBase<BlockType> & G_block_,
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector_,
      const Eigen::MatrixBase<DualVectorType> & dual_vector_)
    {
      const auto & G_block = G_block_.derived();
      auto & primal_vector = primal_vector_.const_cast_derived();
      auto & dual_vector = dual_vector_.const_cast_derived();

      assert(
        primal_vector.size() == dual_vector.size()
        && "The two primal/dual vectors should be of the same size.");
      assert(
        primal_vector.size() == set.size()
        && "The the primal vector should be of the same size than the box set.");
      assert(
        dual_vector.size() == set.size()
        && "The the dual vector should be of the same size than the box set.");

      const Eigen::DenseIndex size = set.size();

      for (Eigen::DenseIndex row_id = 0; row_id < size; ++row_id)
      {
        Scalar & value = primal_vector.coeffRef(row_id);
        const Scalar value_previous = value;
        value -=
          Scalar(
            over_relax_value
            / math::max(Eigen::NumTraits<Scalar>::dummy_precision(), G_block.coeff(row_id, row_id)))
          * dual_vector[row_id];
        value = set.rowiseProject(row_id, value);
        dual_vector.noalias() +=
          G_block.col(row_id)
          * Scalar(value - value_previous); // TODO optimize: we only need dual_vector[row_id] for
                                            // the update and not the full dual vector
      }
    }

    ///
    /// \brief Perform a projection step associated with the PGS algorithm
    ///
    /// \param[in] G_block block asscociated with the current
    /// \param[in,out] primal_vector_ primal vector which will be update with the new estimate
    /// \param[in,out] dual_vector_ dual vector which will be update with the new estimate
    ///
    template<typename BlockType, typename PrimalVectorType, typename DualVectorType>
    void project(
      const Eigen::EigenBase<BlockType> & G_block,
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector,
      const Eigen::MatrixBase<DualVectorType> & dual_vector) const
    {
      project_impl(
        this->set, this->over_relax_value, G_block.derived(), primal_vector.const_cast_derived(),
        dual_vector.const_cast_derived());
    }

    template<
      typename ConstraintSetType,
      typename BlockType,
      typename PrimalVectorType,
      typename DualVectorType>
    static void project_impl(
      const ConstraintSetType & set,
      const Scalar over_relax_value,
      const Eigen::EigenBase<BlockType> & G_block_, // for Sparse matrices
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector_,
      const Eigen::MatrixBase<DualVectorType> & dual_vector_)
    {
      const auto & G_block = G_block_.derived();
      auto & primal_vector = primal_vector_.const_cast_derived();
      auto & dual_vector = dual_vector_.const_cast_derived();

      assert(
        primal_vector.size() == dual_vector.size()
        && "The two primal/dual vectors should be of the same size.");
      assert(
        primal_vector.size() == set.size()
        && "The the primal vector should be of the same size than the box set.");
      assert(
        dual_vector.size() == set.size()
        && "The the dual vector should be of the same size than the box set.");

      const Eigen::DenseIndex size = set.size();

      for (Eigen::DenseIndex row_id = 0; row_id < size; ++row_id)
      {
        Scalar & value = primal_vector.coeffRef(row_id);
        const Scalar value_previous = value;
        value -=
          Scalar(
            over_relax_value
            / math::max(Eigen::NumTraits<Scalar>::dummy_precision(), G_block.coeff(row_id, row_id)))
          * dual_vector[row_id];
        value = set.rowiseProject(row_id, value);
        dual_vector += G_block.col(row_id) * Scalar(value - value_previous);
      }
    }

    /// \brief Compute the feasibility conditions associated with the optimization problem
    template<typename PrimalVectorType, typename DualVectorType>
    void computeFeasibility(
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector,
      const Eigen::MatrixBase<DualVectorType> & dual_vector)
    {
      this->primal_feasibility =
        Scalar(0); // always zero as the primal variable belongs to the constraint set.
      this->dual_feasibility =
        Scalar(0); // always zero as the dual variable belongs to the constraint set.

      const Eigen::DenseIndex size = set.size();
      Scalar complementarity = Scalar(0);

      const auto & lb = set.lb();
      const auto & ub = set.ub();
      for (Eigen::DenseIndex row_id = 0; row_id < size; ++row_id)
      {
        const Scalar dual_positive_part = math::max(Scalar(0), dual_vector[row_id]);
        const Scalar dual_negative_part = dual_positive_part - dual_vector[row_id];

        Scalar row_complementarity = dual_positive_part * (primal_vector[row_id] - lb[row_id]);
        row_complementarity =
          math::max(row_complementarity, dual_negative_part * (ub[row_id] - primal_vector[row_id]));
        complementarity = math::max(complementarity, row_complementarity);
      }
      this->complementarity = complementarity;
    }

    const ConstraintSet & set;

  }; // PGSConstraintProjectionStep<BoxSetTpl<_Scalar>>

  template<typename _Scalar>
  struct PGSConstraintProjectionStep<JointLimitConstraintConeTpl<_Scalar>>
  : PGSConstraintProjectionStepBase<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef JointLimitConstraintConeTpl<Scalar> ConstraintSet;
    typedef BoxSetTpl<Scalar> BoxSet;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef PGSConstraintProjectionStepBase<Scalar> Base;

    PGSConstraintProjectionStep(const Scalar over_relax_value, const ConstraintSet & set)
    : Base(over_relax_value)
    , set(set)
    {
    }

    ///
    /// \brief Perform a projection step associated with the PGS algorithm
    ///
    /// \param[in] G_block block asscociated with the current
    /// \param[in,out] primal_vector_ primal vector which will be update with the new estimate
    /// \param[in,out] dual_vector_ dual vector which will be update with the new estimate
    ///
    template<typename BlockType, typename PrimalVectorType, typename DualVectorType>
    void project(
      const Eigen::MatrixBase<BlockType> & G_block_,
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector_,
      const Eigen::MatrixBase<DualVectorType> & dual_vector_) const
    {
      PGSConstraintProjectionStep<BoxSet>::project_impl(
        this->set, this->over_relax_value, G_block_.derived(), primal_vector_.const_cast_derived(),
        dual_vector_.const_cast_derived());
    }

    ///
    /// \brief Perform a projection step associated with the PGS algorithm
    ///
    /// \param[in] G_block block asscociated with the current
    /// \param[in,out] primal_vector_ primal vector which will be update with the new estimate
    /// \param[in,out] dual_vector_ dual vector which will be update with the new estimate
    ///
    template<typename BlockType, typename PrimalVectorType, typename DualVectorType>
    void project(
      const Eigen::EigenBase<BlockType> & G_block_, // for Sparse matrices
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector_,
      const Eigen::MatrixBase<DualVectorType> & dual_vector_) const
    {
      PGSConstraintProjectionStep<BoxSet>::project_impl(
        this->set, this->over_relax_value, G_block_.derived(), primal_vector_.const_cast_derived(),
        dual_vector_.const_cast_derived());
    }

    /// \brief Compute the feasibility conditions associated with the optimization problem
    template<typename PrimalVectorType, typename DualVectorType>
    void computeFeasibility(
      const Eigen::MatrixBase<PrimalVectorType> & primal_vector,
      const Eigen::MatrixBase<DualVectorType> & dual_vector)
    {
      this->primal_feasibility =
        Scalar(0); // always zero as the primal variable belongs to the constraint set.

      const Eigen::DenseIndex size = set.size();
      Scalar complementarity = Scalar(0);
      Scalar dual_feasibility = Scalar(0);

      for (Eigen::DenseIndex row_id = 0; row_id < size; ++row_id)
      {
        const Scalar row_complementarity =
          math::fabs(Scalar(primal_vector[row_id] * dual_vector[row_id]));
        complementarity = math::max(complementarity, row_complementarity);

        const Scalar row_dual_feasibility =
          math::fabs(dual_vector[row_id] - set.dual().rowiseProject(row_id, dual_vector[row_id]));
        dual_feasibility = math::max(dual_feasibility, row_dual_feasibility);
      }
      this->complementarity = complementarity;
      this->dual_feasibility = dual_feasibility;
    }

    const ConstraintSet & set;

  }; // PGSConstraintProjectionStep<JointLimitConstraintConeTpl<_Scalar>>

  template<typename _Scalar>
  template<
    typename VectorLike,
    template<typename T> class Holder,
    typename ConstraintModel,
    typename ConstraintModelAllocator>
  bool PGSContactSolverTpl<_Scalar>::solve(
    const DelassusOperatorDense & delassus,
    const Eigen::MatrixBase<VectorLike> & g,
    const std::vector<Holder<const ConstraintModel>, ConstraintModelAllocator> & constraint_models,
    const Scalar dt,
    const boost::optional<RefConstVectorXs> x_guess,
    const Scalar over_relax,
    const bool solve_ncp,
    const bool stat_record)

  {
    const MatrixXs & G = delassus.undampedMatrix();
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      over_relax < Scalar(2) && over_relax > Scalar(0), "over_relax should lie in ]0,2[.")
    PINOCCHIO_CHECK_ARGUMENT_SIZE(g.size(), this->getProblemSize());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(G.rows(), this->getProblemSize());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(G.cols(), this->getProblemSize());
    if (x_guess)
    {
      x = x_guess.get();
      PINOCCHIO_CHECK_ARGUMENT_SIZE(x.size(), this->getProblemSize());
    }
    else
    {
      x.setZero();
    }

    const size_t nc = constraint_models.size(); // num constraints

    internal::getTimeScalingFromAccelerationToConstraints(
      constraint_models, dt, time_scaling_acc_to_constraints);
    internal::getTimeScalingFromConstraintsToPosition(
      time_scaling_acc_to_constraints, dt, time_scaling_constraints_to_pos);
    gs = g.array() / time_scaling_acc_to_constraints.array();

    int it = 1;
    PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();

#ifdef PINOCCHIO_WITH_HPP_FCL
    timer.start();
#endif // PINOCCHIO_WITH_HPP_FCL

    Scalar complementarity, proximal_metric, primal_feasibility, dual_feasibility;
    bool abs_prec_reached = false, rel_prec_reached = false;
    Scalar x_previous_norm_inf = x.template lpNorm<Eigen::Infinity>();

    if (stat_record)
    {
      stats.reserve(this->max_it);
      stats.reset();
    }

    for (; it <= this->max_it; ++it)
    {
      x_previous = x;
      complementarity = Scalar(0);
      dual_feasibility = Scalar(0);
      primal_feasibility = Scalar(0);
      Eigen::DenseIndex row_id = 0;
      for (size_t constraint_id = 0; constraint_id < nc; ++constraint_id)
      {
        const ConstraintModel & cmodel = constraint_models[constraint_id];
        const Eigen::DenseIndex constraint_set_size = cmodel.activeSize();

        auto G_block = G.block(row_id, row_id, constraint_set_size, constraint_set_size);
        auto force = x.segment(row_id, constraint_set_size);

        auto velocity = y.segment(row_id, constraint_set_size);
        auto dual_to_pos = y_to_pos.segment(row_id, constraint_set_size);

        auto time_scaling_acc_to_constraints_segment =
          time_scaling_acc_to_constraints.segment(row_id, constraint_set_size);
        auto time_scaling_constraints_to_pos_segment =
          time_scaling_constraints_to_pos.segment(row_id, constraint_set_size);

        // Update dual variable
        velocity.noalias() = G.middleRows(row_id, constraint_set_size) * x;
        velocity += gs.segment(row_id, constraint_set_size);

        typedef PGSConstraintProjectionStepVisitor<
          Scalar, decltype(G_block), decltype(force), decltype(velocity), decltype(dual_to_pos),
          decltype(time_scaling_acc_to_constraints_segment),
          decltype(time_scaling_constraints_to_pos_segment)>
          Step;
        Step step(over_relax);
        step.run(
          cmodel, G_block, force, velocity, dual_to_pos, time_scaling_acc_to_constraints_segment,
          time_scaling_constraints_to_pos_segment);
        //        PGSConstraintProjectionStep<ConstraintSet> step(over_relax, set);
        //        step.project(G_block, velocity, force);
        //        step.computeFeasibility(velocity, force);

        // Update problem feasibility
        complementarity = math::max(complementarity, step.complementarity);
        dual_feasibility = math::max(dual_feasibility, step.dual_feasibility);
        primal_feasibility = math::max(primal_feasibility, step.primal_feasibility);

        // Update row id for the next constraint
        row_id += constraint_set_size;
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

      if (stat_record)
      {
        VectorXs tmp = G * x;
        tmp += gs;
        VectorXs rhs(tmp.size());
        if (solve_ncp)
        {
          internal::computeDeSaxeCorrection(constraint_models, tmp, rhs);
          tmp += rhs;
        }

        tmp.array() *=
          time_scaling_acc_to_constraints.array(); // back to constraint formulation level
        rhs = tmp;
        internal::computeDualConeProjection(constraint_models, rhs, rhs);
        tmp -= rhs;
        Scalar dual_feasibility_ncp = tmp.template lpNorm<Eigen::Infinity>();
        stats.dual_feasibility_ncp.push_back(dual_feasibility_ncp);
        stats.it = it;
        stats.primal_feasibility.push_back(primal_feasibility);
        stats.dual_feasibility.push_back(dual_feasibility);
        stats.complementarity.push_back(complementarity);
      }

      if (abs_prec_reached || rel_prec_reached)
        break;

      x_previous_norm_inf = x_norm_inf;
    }

#ifdef PINOCCHIO_WITH_HPP_FCL
    timer.stop();
#endif // PINOCCHIO_WITH_HPP_FCL

    PINOCCHIO_EIGEN_MALLOC_ALLOWED();

    // Express the dual variable in the units of constraints
    y.array() *= time_scaling_acc_to_constraints.array();

    this->absolute_residual = math::max(complementarity, dual_feasibility);
    this->relative_residual = proximal_metric;
    this->it = it;

    if (abs_prec_reached || rel_prec_reached)
      return true;

    return false;
  }

  template<typename _Scalar>
  template<typename VectorLike, typename ConstraintModel, typename ConstraintModelAllocator>
  bool PGSContactSolverTpl<_Scalar>::solve(
    const DelassusOperatorDense & delassus,
    const Eigen::MatrixBase<VectorLike> & g,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const Scalar dt,
    const boost::optional<RefConstVectorXs> x_guess,
    const Scalar over_relax,
    const bool solve_ncp,
    const bool stat_record)

  {
    typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
    typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

    WrappedConstraintModelVector wrapped_constraint_models(
      constraint_models.cbegin(), constraint_models.cend());

    return solve(
      delassus, g, wrapped_constraint_models, dt, x_guess, over_relax, solve_ncp, stat_record);
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_pgs_solver_hxx__
