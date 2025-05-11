//
// Copyright (c) 2023-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_jointwise_constraint_base_hpp__
#define __pinocchio_algorithm_constraints_jointwise_constraint_base_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"

namespace pinocchio
{

  template<typename Derived>
  struct JointWiseConstraintModelBase : ConstraintModelBase<Derived>
  {
    typedef typename traits<Derived>::Scalar Scalar;
    enum
    {
      Options = traits<Derived>::Options
    };
    typedef ConstraintModelBase<Derived> Base;

    using Base::derived;
    using typename Base::BooleanVector;
    using typename Base::EigenIndexVector;

    typedef typename traits<Derived>::ConstraintData ConstraintData;
    typedef typename traits<Derived>::ConstraintSet ConstraintSet;

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    template<int Options, template<typename, int> class JointCollectionTpl>
    JointWiseConstraintModelBase(const ModelTpl<Scalar, Options, JointCollectionTpl> & model)
    : Base(model)
    {
    }

    template<typename OtherDerived>
    bool operator==(const JointWiseConstraintModelBase<OtherDerived> & other) const
    {
      return base() == other.base();
    }

    template<typename OtherDerived>
    bool operator!=(const JointWiseConstraintModelBase<OtherDerived> & other) const
    {
      return !(*this == other);
    }

    /// \brief Map the constraint forces (aka constraint Lagrange multipliers) to the joint torques
    /// associated to each independant constraint. This operation corresponds to the mapping of the
    /// constraint multipliers on the joint torque.
    ///
    /// \param[in] model The model of the rigid body system.
    /// \param[in] data The data associated with model.
    /// \param[in] cdata The constraint data associated with the constraint model.
    /// \param[in] constraint_forces Input constraint forces (Lagrange multipliers) associated with
    /// the constraint.
    /// \param[out] joint_torques Output joint torques associated with the model.
    ///
    /// \note The results will be added to the joint_torques ouput argument.
    template<
      template<typename, int> class JointCollectionTpl,
      typename ConstraintForcesLike,
      typename JointTorquesLike>
    void mapConstraintForceToJointTorques(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ConstraintForcesLike> & constraint_forces,
      const Eigen::MatrixBase<JointTorquesLike> & joint_torques) const
    {
      derived().mapConstraintForceToJointTorques(
        model, data, cdata, constraint_forces, joint_torques.const_cast_derived());
    }

    ///
    /// \copydoc Base::mapConstraintForceToJointSpace(const ModelTpl<Scalar, Options,
    /// JointCollectionTpl> &, const DataTpl<Scalar, Options, JointCollectionTpl> , const
    /// ConstraintData &, const Eigen::MatrixBase<ConstraintForceLike> &,
    /// std::vector<ForceTpl<Scalar, Options>, ForceAllocator> &, const
    /// Eigen::MatrixBase<JointTorquesLike> &,ReferenceFrameTag<rf>)
    ///
    template<
      template<typename, int> class JointCollectionTpl,
      typename ConstraintForceLike,
      typename ForceAllocator,
      typename JointTorquesLike,
      ReferenceFrame rf>
    void mapConstraintForceToJointSpace(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ConstraintForceLike> & constraint_forces,
      std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces,
      const Eigen::MatrixBase<JointTorquesLike> & joint_torques,
      ReferenceFrameTag<rf> reference_frame) const
    {
      PINOCCHIO_UNUSED_VARIABLE(joint_forces);
      PINOCCHIO_UNUSED_VARIABLE(reference_frame);
      mapConstraintForceToJointTorques(model, data, cdata, constraint_forces, joint_torques);
    }

    /// \brief Map the joint motions to the constraint motions.
    /// This operation corresponds to the dual mapping wrt mapConstraintForcesToJointTorques.
    ///
    /// \param[in] model The model of the rigid body system.
    /// \param[in] data The data associated with model.
    /// \param[in] cdata The constraint data associated with the constraint model.
    /// \param[in] joint_generalized_velocity Input joint motions associated with the model.
    /// \param[out] constraint_motions Output constraint motions.
    ///
    template<
      template<typename, int> class JointCollectionTpl,
      typename JointMotionsLike,
      typename ConstraintMotionsLike>
    void mapJointMotionsToConstraintMotion(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<JointMotionsLike> & joint_generalized_velocity,
      const Eigen::MatrixBase<ConstraintMotionsLike> & constraint_motions) const
    {
      derived().mapJointMotionsToConstraintMotion(
        model, data, cdata, joint_generalized_velocity, constraint_motions.const_cast_derived());
    }

    ///
    ///\copydoc Base::mapJointSpaceToConstraintMotion(const ModelTpl<Scalar, Options,
    /// JointCollectionTpl> &, const DataTpl<Scalar, Options, JointCollectionTpl> , const
    /// ConstraintData &,
    /// std::vector<MotionTpl<Scalar, Options>, MotionAllocator> &, const
    /// Eigen::MatrixBase<JointMotionsLike> &, const Eigen::MatrixBase<VectorLike>
    /// &,ReferenceFrameTag<rf>)
    ///
    template<
      template<typename, int> class JointCollectionTpl,
      typename MotionAllocator,
      typename JointMotionsLike,
      typename VectorLike,
      ReferenceFrame rf>
    void mapJointSpaceToConstraintMotion(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_motions,
      const Eigen::MatrixBase<JointMotionsLike> & joint_generalized_velocity,
      const Eigen::MatrixBase<VectorLike> & constraint_motions,
      ReferenceFrameTag<rf> reference_frame) const
    {
      PINOCCHIO_UNUSED_VARIABLE(joint_motions);
      PINOCCHIO_UNUSED_VARIABLE(reference_frame);
      mapJointMotionsToConstraintMotion(
        model, data, cdata, joint_generalized_velocity, constraint_motions.const_cast_derived());
    }

  protected:
    /// \brief Default constructor
    JointWiseConstraintModelBase()
    {
    }
  }; // struct JointWiseConstraintModelBase
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_jointwise_constraint_base_hpp__
