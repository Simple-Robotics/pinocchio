//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_contact_solver_utils_hpp__
#define __pinocchio_algorithm_contact_solver_utils_hpp__

#include "pinocchio/math/fwd.hpp"
#include "pinocchio/math/comparison-operators.hpp"
#include "pinocchio/algorithm/constraints/constraints.hpp"
#include "pinocchio/algorithm/delassus-operator-base.hpp"
#include "pinocchio/algorithm/constraints/visitors/constraint-model-visitor.hpp"

namespace pinocchio
{

  namespace internal
  {

    template<typename ForceVectorLike, typename ResultVectorLike>
    struct ProjectionVisitor
    : visitors::ConstraintUnaryVisitorBase<ProjectionVisitor<ForceVectorLike, ResultVectorLike>>
    {

      typedef boost::fusion::vector<const ForceVectorLike &, ResultVectorLike &> ArgsType;

      typedef visitors::ConstraintUnaryVisitorBase<
        ProjectionVisitor<ForceVectorLike, ResultVectorLike>>
        Base;

      template<typename ConstraintModel, typename Vector1Like, typename Vector2Like>
      static void algo(
        const ConstraintModelBase<ConstraintModel> & cmodel,
        const Eigen::MatrixBase<Vector1Like> & force,
        const Eigen::MatrixBase<Vector2Like> & result)
      {
        result.const_cast_derived() = cmodel.set().project(force.derived());
        //        assert(set.dual().isInside(result, Scalar(1e-12)));
      }

      using Base::run;

      template<typename ConstraintModel, typename Vector1Like, typename Vector2Like>
      static void run(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const Eigen::MatrixBase<Vector1Like> & force,
        const Eigen::MatrixBase<Vector2Like> & result)
      {
        algo(cmodel.derived(), force.derived(), result.const_cast_derived());
      }

      template<
        typename Scalar,
        int Options,
        template<typename S, int O>
        class ConstraintCollectionTpl,
        typename Vector1Like,
        typename Vector2Like>
      static void run(
        const pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        const Eigen::MatrixBase<Vector1Like> & force,
        const Eigen::MatrixBase<Vector2Like> & result)
      {
        ArgsType args(force.derived(), result.const_cast_derived());
        run(cmodel.derived(), args);
      }
    };

    /// \brief Project a vector x on the vector of cones.
    template<
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename VectorLikeIn,
      typename VectorLikeOut>
    void computeConeProjection(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const Eigen::DenseBase<VectorLikeIn> & x,
      const Eigen::DenseBase<VectorLikeOut> & x_proj_)
    {
      assert(x.size() == x_proj_.size());
      Eigen::DenseIndex index = 0;
      VectorLikeOut & x_proj = x_proj_.const_cast_derived();
      for (const auto & cmodel : constraint_models)
      {
        const auto size = cmodel.size();
        auto res = x_proj.segment(index, size);
        const auto force_segment = x.segment(index, size);

        typedef ProjectionVisitor<decltype(force_segment), decltype(res)> Algo;
        Algo::run(cmodel, force_segment, res);

        index += size;
      }
    }

    template<typename VelocityVectorLike, typename ResultVectorLike>
    struct DualProjectionVisitor
    : visitors::ConstraintUnaryVisitorBase<
        DualProjectionVisitor<VelocityVectorLike, ResultVectorLike>>
    {
      typedef boost::fusion::vector<const VelocityVectorLike &, ResultVectorLike &> ArgsType;

      typedef visitors::ConstraintUnaryVisitorBase<
        DualProjectionVisitor<VelocityVectorLike, ResultVectorLike>>
        Base;

      template<typename ConstraintModel, typename Vector1Like, typename Vector2Like>
      static void algo(
        const ConstraintModelBase<ConstraintModel> & cmodel,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & result)
      {
        algo_step(cmodel.set(), velocity.derived(), result.const_cast_derived());
      }

      template<typename Vector1Like, typename Vector2Like>
      static void algo_step(
        const CoulombFrictionConeTpl<double> & set,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & result)
      {
        result.const_cast_derived() = set.dual().project(velocity);
        //        assert(set.dual().isInside(result, Scalar(1e-12)));
      }

      template<typename ConstraintSet, typename Vector1Like, typename Vector2Like>
      static void algo_step(
        const ConstraintSet & set,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & result)
      {
        PINOCCHIO_UNUSED_VARIABLE(set);
        result.const_cast_derived() = velocity;
      }

      using Base::run;

      template<typename ConstraintModel, typename Vector1Like, typename Vector2Like>
      static void run(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & result)
      {
        algo(cmodel.derived(), velocity.derived(), result.const_cast_derived());
      }

      template<
        typename Scalar,
        int Options,
        template<typename S, int O>
        class ConstraintCollectionTpl,
        typename Vector1Like,
        typename Vector2Like>
      static void run(
        const pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & result)
      {
        ArgsType args(velocity.derived(), result.const_cast_derived());
        run(cmodel.derived(), args);
      }
    }; // struct DualProjectionVisitor

    /// \brief Project a vector x on the dual of the cones contained in the vector of cones.
    template<
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename VectorLikeIn,
      typename VectorLikeOut>
    void computeDualConeProjection(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const Eigen::DenseBase<VectorLikeIn> & x,
      const Eigen::DenseBase<VectorLikeOut> & x_proj_)
    {
      assert(x.size() == x_proj_.size());
      VectorLikeOut & x_proj = x_proj_.const_cast_derived();
      Eigen::DenseIndex index = 0;
      for (const auto & cmodel : constraint_models)
      {
        const auto size = cmodel.size();
        const auto velocity_segment = x.segment(index, size);
        auto res_segment = x_proj.segment(index, size);
        typedef DualProjectionVisitor<decltype(velocity_segment), decltype(res_segment)> Algo;
        Algo::run(cmodel, velocity_segment, res_segment);
        index += size;
      }
    }

    template<typename Scalar, typename VelocityVectorLike, typename ForceVectorLike>
    struct ComplementarityVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ComplementarityVisitor<Scalar, VelocityVectorLike, ForceVectorLike>,
        Scalar>
    {
      typedef boost::fusion::vector<const VelocityVectorLike &, const ForceVectorLike &> ArgsType;

      typedef visitors::ConstraintUnaryVisitorBase<
        ComplementarityVisitor<Scalar, VelocityVectorLike, ForceVectorLike>,
        Scalar>
        Base;
      using Base::run;

      template<typename ConstraintModel, typename Vector1Like, typename Vector2Like>
      static Scalar algo(
        const ConstraintModelBase<ConstraintModel> & cmodel,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & force)
      {
        return algo_step(cmodel.set(), velocity.derived(), force.derived());
      }

      template<typename Vector1Like, typename Vector2Like>
      static Scalar algo_step(
        const CoulombFrictionConeTpl<Scalar> & set,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & force)
      {
        return set.computeConicComplementarity(velocity, force);
      }

      template<typename Vector1Like, typename Vector2Like>
      static Scalar algo_step(
        const UnboundedSetTpl<Scalar> & set,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & force)
      {
        PINOCCHIO_UNUSED_VARIABLE(set);
        PINOCCHIO_UNUSED_VARIABLE(velocity);
        PINOCCHIO_UNUSED_VARIABLE(force);
        return Scalar(0);
      }

      template<typename Vector1Like, typename Vector2Like>
      static Scalar algo_step(
        const BoxSetTpl<Scalar> & set,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & force)
      {
        Scalar complementarity = Scalar(0);

        const auto & lb = set.lb();
        const auto & ub = set.ub();
        for (Eigen::DenseIndex row_id = 0; row_id < set.size(); ++row_id)
        {
          const Scalar velocity_positive_part = math::max(Scalar(0), velocity[row_id]);
          const Scalar velocity_negative_part = velocity_positive_part - velocity[row_id];

          Scalar row_complementarity = velocity_positive_part * (force[row_id] - lb[row_id]);
          row_complementarity =
            math::max(row_complementarity, velocity_negative_part * (ub[row_id] - force[row_id]));
          complementarity = math::max(complementarity, row_complementarity);
        }

        return complementarity;
      }

      template<typename ConstraintModel, typename Vector1Like, typename Vector2Like>
      static Scalar run(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & force)
      {
        return algo(cmodel.derived(), velocity.derived(), force.derived());
      }

      template<
        int Options,
        template<typename S, int O>
        class ConstraintCollectionTpl,
        typename Vector1Like,
        typename Vector2Like>
      static Scalar run(
        const pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & force)
      {
        ArgsType args(velocity.derived(), force.derived());
        return run(cmodel.derived(), args);
      }
    }; // struct ComplementarityVisitor

    template<
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename VectorLikeVelocity,
      typename VectorLikeForce>
    typename ConstraintModel::Scalar computeConicComplementarity(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const Eigen::DenseBase<VectorLikeVelocity> & velocities,
      const Eigen::DenseBase<VectorLikeForce> & forces)
    {
      typedef typename ConstraintModel::Scalar Scalar;
      assert(velocities.size() == forces.size());
      Eigen::DenseIndex index = 0;
      Scalar complementarity = Scalar(0);
      for (const auto & cmodel : constraint_models)
      {
        const auto size = cmodel.size();

        const auto velocity_segment = velocities.segment(index, size);
        const auto force_segment = forces.segment(index, size);

        typedef ComplementarityVisitor<Scalar, decltype(velocity_segment), decltype(force_segment)>
          Algo;
        const Scalar constraint_complementarity =
          Algo::run(cmodel, velocity_segment, force_segment);

        complementarity = math::max(complementarity, constraint_complementarity);
        index += size;
      }

      return complementarity;
    }

    template<typename VelocityVectorLike, typename ResultVectorLike>
    struct DeSaxeCorrectionVisitor
    : visitors::ConstraintUnaryVisitorBase<
        DeSaxeCorrectionVisitor<VelocityVectorLike, ResultVectorLike>>
    {

      typedef boost::fusion::vector<const VelocityVectorLike &, ResultVectorLike &> ArgsType;

      typedef visitors::ConstraintUnaryVisitorBase<
        DeSaxeCorrectionVisitor<VelocityVectorLike, ResultVectorLike>>
        Base;
      using Base::run;

      template<typename ConstraintModel, typename Vector1Like, typename Vector2Like>
      static void algo(
        const ConstraintModelBase<ConstraintModel> & cmodel,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & result)
      {
        return algo_step(cmodel.set(), velocity.derived(), result.const_cast_derived());
      }

      template<typename Vector1Like, typename Vector2Like>
      static void algo_step(
        const CoulombFrictionConeTpl<double> & set,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & result)
      {
        result.const_cast_derived() = set.computeNormalCorrection(velocity);
      }

      template<typename ConstraintSet, typename Vector1Like, typename Vector2Like>
      static void algo_step(
        const ConstraintSet & set,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & result)
      {
        PINOCCHIO_UNUSED_VARIABLE(set);
        PINOCCHIO_UNUSED_VARIABLE(velocity);
        result.const_cast_derived().setZero();
      }

      template<typename ConstraintModel, typename Vector1Like, typename Vector2Like>
      static void run(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & result)
      {
        algo(cmodel.derived(), velocity.derived(), result.const_cast_derived());
      }

      template<
        typename Scalar,
        int Options,
        template<typename S, int O>
        class ConstraintCollectionTpl,
        typename Vector1Like,
        typename Vector2Like>
      static void run(
        const pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        const Eigen::MatrixBase<Vector1Like> & velocity,
        const Eigen::MatrixBase<Vector2Like> & result)
      {
        ArgsType args(velocity.derived(), result.const_cast_derived());
        run(cmodel.derived(), args);
      }
    }; // struct DeSaxeCorrectionVisitor

    template<
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename VectorLikeIn,
      typename VectorLikeOut>
    void computeComplementarityShift(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const Eigen::DenseBase<VectorLikeIn> & velocities,
      const Eigen::DenseBase<VectorLikeOut> & shift_)
    {
      assert(velocities.size() == shift_.size());
      Eigen::DenseIndex index = 0;
      VectorLikeOut & shift = shift_.const_cast_derived();
      for (const auto & cmodel : constraint_models)
      {
        const auto size = cmodel.size();

        const auto velocity_segment = velocities.segment(index, size);
        auto result_segment = shift.segment(index, size);
        typedef DeSaxeCorrectionVisitor<decltype(velocity_segment), decltype(result_segment)> Step;

        Step::run(cmodel, velocity_segment, result_segment);

        index += size;
      }
    }

    template<typename ConstraintSet, typename ConstraintAllocator, typename VectorLikeIn>
    typename ConstraintSet::Scalar computePrimalFeasibility(
      const std::vector<ConstraintSet, ConstraintAllocator> & sets,
      const Eigen::DenseBase<VectorLikeIn> & forces)
    {
      typedef typename ConstraintSet::Vector3 Vector3;
      typedef typename ConstraintSet::Scalar Scalar;

      Eigen::DenseIndex index = 0;
      Scalar norm(0);
      for (const auto & set : sets)
      {
        const auto size = set.size();
        Vector3 df_projected =
          set.project(forces.segment(index, size)) - forces.segment(index, size);
        norm = math::max(norm, df_projected.norm());
        index += size;
      }

      return norm;
    }

    template<
      typename ConstraintSet,
      typename ConstraintAllocator,
      typename ForceVector,
      typename VelocityVector>
    typename ConstraintSet::Scalar computeReprojectionError(
      const std::vector<ConstraintSet, ConstraintAllocator> & sets,
      const Eigen::DenseBase<ForceVector> & forces,
      const Eigen::DenseBase<VelocityVector> & velocities)
    {
      typedef typename ConstraintSet::Vector3 Vector3;
      typedef typename ConstraintSet::Scalar Scalar;

      Eigen::DenseIndex index = 0;
      Scalar norm(0);
      for (const auto & set : sets)
      {
        const auto size = set.size();
        const Vector3 df_projected =
          forces.segment(index, size)
          - set.project(forces.segment(index, size) - velocities.segment(index, size));
        norm = math::max(norm, df_projected.norm());
        index += size;
      }

      return norm;
    }

    // template<typename Scalar, typename ConstraintAllocator, typename VectorLikeIn>
    // Scalar computeDualFeasibility(DelassusOperatorBase<DelassusDerived> & delassus,
    //                               const Eigen::MatrixBase<VectorLike> & g,
    //                               const
    //                               std::vector<CoulombFrictionConeTpl<Scalar>,ConstraintAllocator>
    //                               & cones, const Eigen::DenseBase<VectorLikeIn> & forces)
    //{
    //   typedef CoulombFrictionConeTpl<Scalar> Cone;
    //   typedef typename Cone::Vector3 Vector3;
    //
    //   Eigen::DenseIndex index = 0;
    //   Scalar norm = 0;
    //   for(const auto & cone: cones)
    //   {
    //     const Vector3 df_projected = cone.project(forces.template segment<3>(index)) -
    //     forces.template segment<3>(index); norm = math::max(complementarity,
    //     df_projected.norm()); index += 3;
    //   }
    //
    //   return norm;
    // }

  } // namespace internal

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_contact_solver_utils_hpp__
