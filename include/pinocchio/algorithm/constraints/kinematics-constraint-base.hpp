//
// Copyright (c) 2023-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_kinematics_constraint_base_hpp__
#define __pinocchio_algorithm_constraints_kinematics_constraint_base_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"

namespace pinocchio
{

  template<typename Derived>
  struct KinematicsConstraintModelBase : ConstraintModelBase<Derived>
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

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    template<int Options, template<typename, int> class JointCollectionTpl>
    KinematicsConstraintModelBase(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const JointIndex joint2_id)
    : Base(model)
    , joint1_id(joint1_id)
    , joint2_id(joint2_id)
    {
    }

    /// \brief Index of the first joint in the model tree
    JointIndex joint1_id;

    /// \brief Index of the second joint in the model tree
    JointIndex joint2_id;

    template<typename OtherDerived>
    bool operator==(const KinematicsConstraintModelBase<OtherDerived> & other) const
    {
      return base() == other.base() && joint1_id == other.joint1_id && joint2_id == other.joint2_id;
    }

    template<typename OtherDerived>
    bool operator!=(const KinematicsConstraintModelBase<OtherDerived> & other) const
    {
      return !(*this == other);
    }

    ///
    /// \brief Map the constraint forces (aka constraint Lagrange multipliers) to the forces
    /// supported by the joints expressed in the input reference_frame.
    ///
    /// \param[in] model The model of the rigid body system.
    /// \param[in] data The data associated with model.
    /// \param[in] cdata The constraint data associated with the constraint model.
    /// \param[in] constraint_forces Input constraint forces (Lagrange multipliers) associated with
    /// the constraint.
    /// \param[out] joint_forces Output joint forces associated with each joint of the model.
    /// \param[in] reference_frame Input reference frame in which the forces are expressed.
    ///
    /// \note The results will be added to the joint_torques ouput argument.
    ///
    template<
      template<typename, int> class JointCollectionTpl,
      typename ForceLike,
      typename ForceAllocator,
      ReferenceFrame rf>
    void mapConstraintForceToJointForces(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ForceLike> & constraint_forces,
      std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces,
      ReferenceFrameTag<rf> reference_frame) const
    {
      derived().mapConstraintForceToJointForces(
        model, data, cdata, constraint_forces, joint_forces.const_cast_derived(), reference_frame);
    }

    ///
    /// \brief Map the constraint forces (aka constraint Lagrange multipliers) to the forces
    /// supported by the joints expressed in the LOCAL frame of the joints.
    ///
    /// \param[in] model The model of the rigid body system.
    /// \param[in] data The data associated with model.
    /// \param[in] cdata The constraint data associated with the constraint model.
    /// \param[in] constraint_forces Input constraint forces (Lagrange multipliers) associated with
    /// the constraint.
    /// \param[out] joint_forces Output joint forces associated with each joint of the model.
    ///
    /// \note The results will be added to the joint_torques ouput argument.
    ///
    template<
      template<typename, int> class JointCollectionTpl,
      typename ForceLike,
      typename ForceAllocator>
    void mapConstraintForceToJointForces(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ForceLike> & constraint_forces,
      std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces) const
    {
      mapConstraintForceToJointForces(
        model, data, cdata, constraint_forces, joint_forces.const_cast_derived(), LocalFrameTag());
    }

    ///
    /// \brief Map the joint motions to the constraint motions. The joint motions are expressed in
    /// the frame given by the input argument reference_frame.
    ///
    /// \param[in] model The model of the rigid body system.
    /// \param[in] data The data associated with model.
    /// \param[in] cdata The constraint data associated with the constraint model.
    /// \param[in] joint_motions Input joint motions associated with each joint of the model.
    /// the constraint.
    /// \param[out] constraint_motions Output contraint motions.
    /// \param[in] reference_frame Input reference frame in which the joint motion quantities are
    /// expressed.
    ///
    template<
      template<typename, int> class JointCollectionTpl,
      typename MotionAllocator,
      typename VectorLike,
      ReferenceFrame rf>
    void mapJointMotionsToConstraintMotion(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_motions,
      const Eigen::MatrixBase<VectorLike> & constraint_motions,
      ReferenceFrameTag<rf> reference_frame) const
    {
      derived().mapJointMotionsToConstraintMotion(
        model, data, cdata, joint_motions, constraint_motions.const_cast_derived(),
        reference_frame);
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
      PINOCCHIO_UNUSED_VARIABLE(joint_torques);
      mapConstraintForceToJointForces(
        model, data, cdata, constraint_forces, joint_forces, reference_frame);
    }

    ///
    /// \brief Map the joint motions to the constraint motions expressed.
    ///
    /// \param[in] model The model of the rigid body system.
    /// \param[in] data The data associated with model.
    /// \param[in] cdata The constraint data associated with the constraint model.
    /// \param[in] joint_motions Input joint motions associated with each joint of the model.
    /// the constraint.
    /// \param[out] constraint_motions Output contraint motions.
    ///
    template<
      template<typename, int> class JointCollectionTpl,
      typename MotionAllocator,
      typename VectorLike>
    void mapJointMotionsToConstraintMotion(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_motions,
      const Eigen::MatrixBase<VectorLike> & constraint_motions) const
    {
      derived().mapJointMotionsToConstraintMotion(
        model, data, cdata, joint_motions, constraint_motions.const_cast_derived(),
        LocalFrameTag());
    }

  protected:
    /// \brief Default constructor
    KinematicsConstraintModelBase()
    {
    }
  }; // struct KinematicsConstraintModelBase
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_kinematics_constraint_base_hpp__
