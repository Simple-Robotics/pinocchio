//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_contact_solver_utils_hpp__
#define __pinocchio_algorithm_contact_solver_utils_hpp__

#include "pinocchio/math/fwd.hpp"
#include "pinocchio/math/comparison-operators.hpp"
#include "pinocchio/algorithm/constraints/coulomb-friction-cone.hpp"
#include "pinocchio/algorithm/delassus-operator-base.hpp"

namespace pinocchio
{

  namespace internal
  {

    /// \brief Project a vector x on the vector of cones.
    template<
      typename ConstraintSet,
      typename ConstraintAllocator,
      typename VectorLikeIn,
      typename VectorLikeOut>
    void computeConeProjection(
      const std::vector<ConstraintSet, ConstraintAllocator> & sets,
      const Eigen::DenseBase<VectorLikeIn> & x,
      const Eigen::DenseBase<VectorLikeOut> & x_proj_)
    {
      assert(x.size() == x_proj_.size());
      Eigen::DenseIndex index = 0;
      VectorLikeOut & x_proj = x_proj_.const_cast_derived();
      for (const auto & set : sets)
      {
        const auto size = set.size();
        x_proj.segment(index, size) = set.project(x.segment(index, size));
        index += size;
      }
    }

    struct DualProjectionVisitor
    {
      template<typename VelocityVectorLike, typename ResultVectorLike>
      static void algo(
        const CoulombFrictionConeTpl<double> & set,
        const Eigen::MatrixBase<VelocityVectorLike> & velocity,
        const Eigen::MatrixBase<ResultVectorLike> & result)
      {
        result.const_cast_derived() = set.dual().project(velocity);
        //        assert(set.dual().isInside(result, Scalar(1e-12)));
      }

      template<typename ConstraintSet, typename VelocityVectorLike, typename ResultVectorLike>
      static void algo(
        const ConstraintSet & set,
        const Eigen::MatrixBase<VelocityVectorLike> & velocity,
        const Eigen::MatrixBase<ResultVectorLike> & result)
      {
        PINOCCHIO_UNUSED_VARIABLE(set);
        result.const_cast_derived() = velocity;
      }
    };

    /// \brief Project a vector x on the dual of the cones contained in the vector of cones.
    template<
      typename ConstraintSet,
      typename ConstraintAllocator,
      typename VectorLikeIn,
      typename VectorLikeOut>
    void computeDualConeProjection(
      const std::vector<ConstraintSet, ConstraintAllocator> & sets,
      const Eigen::DenseBase<VectorLikeIn> & x,
      const Eigen::DenseBase<VectorLikeOut> & x_proj_)
    {
      assert(x.size() == x_proj_.size());
      VectorLikeOut & x_proj = x_proj_.const_cast_derived();
      Eigen::DenseIndex index = 0;
      for (const auto & set : sets)
      {
        const auto size = set.size();
        DualProjectionVisitor::algo(set, x.segment(index, size), x_proj.segment(index, size));
        index += size;
      }
    }

    struct ComplementarityVisitor
    {
      template<typename Scalar, typename VelocityVectorLike, typename ForceVectorLike>
      static Scalar algo(
        const CoulombFrictionConeTpl<Scalar> & set,
        const Eigen::MatrixBase<VelocityVectorLike> & velocity,
        const Eigen::MatrixBase<ForceVectorLike> & force)
      {
        return set.computeConicComplementarity(velocity, force);
      }

      template<typename Scalar, typename VelocityVectorLike, typename ForceVectorLike>
      static Scalar algo(
        const UnboundedSetTpl<Scalar> & set,
        const Eigen::MatrixBase<VelocityVectorLike> & velocity,
        const Eigen::MatrixBase<ForceVectorLike> & force)
      {
        PINOCCHIO_UNUSED_VARIABLE(set);
        PINOCCHIO_UNUSED_VARIABLE(velocity);
        PINOCCHIO_UNUSED_VARIABLE(force);
        return Scalar(0);
      }

      template<typename Scalar, typename VelocityVectorLike, typename ForceVectorLike>
      static Scalar algo(
        const BoxSetTpl<Scalar> & set,
        const Eigen::MatrixBase<VelocityVectorLike> & velocity,
        const Eigen::MatrixBase<ForceVectorLike> & force)
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
    };

    template<
      typename ConstraintSet,
      typename ConstraintAllocator,
      typename VectorLikeVelocity,
      typename VectorLikeForce>
    typename ConstraintSet::Scalar computeConicComplementarity(
      const std::vector<ConstraintSet, ConstraintAllocator> & sets,
      const Eigen::DenseBase<VectorLikeVelocity> & velocities,
      const Eigen::DenseBase<VectorLikeForce> & forces)
    {
      typedef typename ConstraintSet::Scalar Scalar;
      assert(velocities.size() == forces.size());
      Eigen::DenseIndex index = 0;
      Scalar complementarity = Scalar(0);
      for (const auto & set : sets)
      {
        const auto size = set.size();
        const Scalar set_complementarity = ComplementarityVisitor::algo(
          set, velocities.segment(index, size), forces.segment(index, size));
        complementarity = math::max(complementarity, set_complementarity);
        index += size;
      }

      return complementarity;
    }

    struct DeSaxeCorrectionVisitor
    {
      template<typename VelocityVectorLike, typename ResultVectorLike>
      static void algo(
        const CoulombFrictionConeTpl<double> & set,
        const Eigen::MatrixBase<VelocityVectorLike> & velocity,
        const Eigen::MatrixBase<ResultVectorLike> & result)
      {
        result.const_cast_derived() = set.computeNormalCorrection(velocity);
      }

      template<typename ConstraintSet, typename VelocityVectorLike, typename ResultVectorLike>
      static void algo(
        const ConstraintSet & set,
        const Eigen::MatrixBase<VelocityVectorLike> & velocity,
        const Eigen::MatrixBase<ResultVectorLike> & result)
      {
        PINOCCHIO_UNUSED_VARIABLE(set);
        PINOCCHIO_UNUSED_VARIABLE(velocity);
        result.const_cast_derived().setZero();
      }
    };

    template<
      typename ConstraintSet,
      typename ConstraintAllocator,
      typename VectorLikeIn,
      typename VectorLikeOut>
    void computeComplementarityShift(
      const std::vector<ConstraintSet, ConstraintAllocator> & sets,
      const Eigen::DenseBase<VectorLikeIn> & velocities,
      const Eigen::DenseBase<VectorLikeOut> & shift_)
    {
      assert(velocities.size() == shift_.size());
      Eigen::DenseIndex index = 0;
      VectorLikeOut & shift = shift_.const_cast_derived();
      for (const auto & set : sets)
      {
        const auto size = set.size();
        DeSaxeCorrectionVisitor::algo(
          set, velocities.segment(index, size), shift.segment(index, size));
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
