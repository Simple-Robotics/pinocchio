//
// Copyright (c) 2019-2024 INRIA CNRS
//

#ifndef __pinocchio_algorithm_constraints_point_constraint_data_hpp__
#define __pinocchio_algorithm_constraints_point_constraint_data_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"

namespace pinocchio
{

  ///
  ///  \brief Data structure associated with PointConstraint
  ///
  template<typename Derived>
  struct PointConstraintDataBase : ConstraintDataBase<Derived>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef typename traits<Derived>::Scalar Scalar;
    enum
    {
      Options = traits<Derived>::Options
    };

    typedef typename traits<Derived>::ConstraintModel ConstraintModel;
    typedef typename traits<Derived>::ConstraintData ConstraintData;
    typedef ConstraintDataBase<Derived> Base;

    typedef SE3Tpl<Scalar, Options> SE3;
    typedef MotionTpl<Scalar, Options> Motion;
    typedef ForceTpl<Scalar, Options> Force;
    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;
    typedef Eigen::Matrix<Scalar, 6, 6, Options> Matrix6;
    typedef PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(Matrix6) VectorOfMatrix6;
    typedef Eigen::Matrix<Scalar, 6, Eigen::Dynamic, Options> Matrix6x;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> MatrixX;

    // data

    /// \brief Resulting contact forces
    Vector3 constraint_force;

    /// \brief Placement of the constraint frame 1 with respect to the WORLD frame
    SE3 oMc1;

    /// \brief Placement of the constraint frame 2 with respect to the WORLD frame
    SE3 oMc2;

    /// \brief Relative displacement between the two frames
    SE3 c1Mc2;

    /// \brief Constraint position error
    Vector3 constraint_position_error;

    /// \brief Constraint velocity error
    Vector3 constraint_velocity_error;

    /// \brief Constraint acceleration error
    Vector3 constraint_acceleration_error;

    /// \brief Constraint acceleration biais
    Vector3 constraint_acceleration_biais_term;

    //    VectorOfMatrix6 extended_motion_propagators_joint1;
    //    VectorOfMatrix6 lambdas_joint1;
    //    VectorOfMatrix6 extended_motion_propagators_joint2;

    //    Matrix6x dv1_dq, da1_dq, da1_dv, da1_da;
    //    Matrix6x dv2_dq, da2_dq, da2_dv, da2_da;
    //    MatrixX dvc_dq, dac_dq, dac_dv, dac_da;

    /// \brief Default constructor
    PointConstraintDataBase()
    {
    }

    explicit PointConstraintDataBase(const ConstraintModel & constraint_model)
    : constraint_force(Vector3::Zero())
    , oMc1(SE3::Identity())
    , oMc2(SE3::Identity())
    , c1Mc2(SE3::Identity())
    , constraint_position_error(Vector3::Zero())
    , constraint_velocity_error(Vector3::Zero())
    , constraint_acceleration_error(Vector3::Zero())
    , constraint_acceleration_biais_term(Vector3::Zero())
    //    , extended_motion_propagators_joint1(constraint_model.depth_joint1, Matrix6::Zero())
    //    , lambdas_joint1(constraint_model.depth_joint1, Matrix6::Zero())
    //    , extended_motion_propagators_joint2(constraint_model.depth_joint2, Matrix6::Zero())
    //    , dv1_dq(6, constraint_model.nv)
    //    , da1_dq(6, constraint_model.nv)
    //    , da1_dv(6, constraint_model.nv)
    //    , da1_da(6, constraint_model.nv)
    //    , dv2_dq(6, constraint_model.nv)
    //    , da2_dq(6, constraint_model.nv)
    //    , da2_dv(6, constraint_model.nv)
    //    , da2_da(6, constraint_model.nv)
    //    , dvc_dq(constraint_model.size(), constraint_model.nv)
    //    , dac_dq(constraint_model.size(), constraint_model.nv)
    //    , dac_dv(constraint_model.size(), constraint_model.nv)
    //    , dac_da(constraint_model.size(), constraint_model.nv)
    {
      PINOCCHIO_UNUSED_VARIABLE(constraint_model);
    }

    bool operator==(const PointConstraintDataBase & other) const
    {
      return constraint_force == other.constraint_force && oMc1 == other.oMc1 && oMc2 == other.oMc2
             && c1Mc2 == other.c1Mc2 && constraint_position_error == other.constraint_position_error
             && constraint_velocity_error == other.constraint_velocity_error
             && constraint_acceleration_error == other.constraint_acceleration_error
             && constraint_acceleration_biais_term == other.constraint_acceleration_biais_term
        //      && extended_motion_propagators_joint1 == other.extended_motion_propagators_joint1
        //      && lambdas_joint1 == other.lambdas_joint1
        //      && extended_motion_propagators_joint2 == other.extended_motion_propagators_joint2
        //
        //      && dv1_dq == other.dv1_dq && da1_dq == other.da1_dq && da1_dv == other.da1_dv
        //      && da1_da == other.da1_da
        //        //
        //      && dv2_dq == other.dv2_dq && da2_dq == other.da2_dq && da2_dv == other.da2_dv
        //      && da2_da == other.da2_da
        //        //
        //      && dvc_dq == other.dvc_dq && dac_dq == other.dac_dq && dac_dv == other.dac_dv
        //      && dac_da == other.dac_da
        ;
    }

    bool operator!=(const PointConstraintDataBase & other) const
    {
      return !(*this == other);
    }

    using Base::derived;

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }
  };

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_point_constraint_data_hpp__
