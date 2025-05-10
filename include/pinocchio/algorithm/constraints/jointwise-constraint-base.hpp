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
    void mapConstraintForcesToJointTorques(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ConstraintForcesLike> & constraint_forces,
      const Eigen::MatrixBase<JointTorquesLike> & joint_torques) const
    {
      derived().mapConstraintForcesToJointTorques(
        model, data, cdata, constraint_forces, joint_torques.const_cast_derived());
    }

    /// \brief Map the joint motions to the constraint motions.
    /// This operation corresponds to the dual mapping wrt mapConstraintForcesToJointTorques.
    ///
    /// \param[in] model The model of the rigid body system.
    /// \param[in] data The data associated with model.
    /// \param[in] cdata The constraint data associated with the constraint model.
    /// \param[in] joint_motions Input joint motions associated with the model.
    /// \param[out] constraint_motions Output constraint motions.
    ///
    template<
      template<typename, int> class JointCollectionTpl,
      typename JointMotionsLike,
      typename ConstraintMotionsLike>
    void mapJointMotionsToConstraintMotions(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<JointMotionsLike> & joint_motions,
      const Eigen::MatrixBase<ConstraintMotionsLike> & constraint_motions) const
    {
      derived().mapJointMotionsToConstraintMotions(
        model, data, cdata, joint_motions, constraint_motions.const_cast_derived());
    }

  protected:
    /// \brief Default constructor
    JointWiseConstraintModelBase()
    {
    }
  }; // struct JointWiseConstraintModelBase
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_jointwise_constraint_base_hpp__
