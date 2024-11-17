//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__
#define __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__

#include "pinocchio/math/fwd.hpp"

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"

namespace pinocchio
{

  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, JointLimitConstraintModelTpl<Scalar, Options>>
  {
    typedef JointLimitConstraintModelTpl<NewScalar, Options> type;
  };

  template<typename _Scalar, int _Options>
  struct traits<JointLimitConstraintModelTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef JointLimitConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef JointLimitConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef BoxSetTpl<Scalar, Options> ConstraintSet;
  };

  template<typename _Scalar, int _Options>
  struct traits<JointLimitConstraintDataTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef JointLimitConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef JointLimitConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef BoxSetTpl<Scalar, Options> ConstraintSet;
  };

  template<typename _Scalar, int _Options>
  struct JointLimitConstraintModelTpl
  : ConstraintModelBase<JointLimitConstraintModelTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef ConstraintModelBase<JointLimitConstraintModelTpl> Base;
    typedef std::vector<JointIndex> JointIndexVector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;

    using typename Base::BooleanVector;
    using typename Base::EigenIndexVector;

    typedef std::vector<BooleanVector> VectorOfBooleanVector;
    typedef std::vector<EigenIndexVector> VectofOfEigenIndexVector;

    typedef JointLimitConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef BoxSetTpl<Scalar, Options> ConstraintSet;
    typedef BoxSetTpl<Scalar, Options> BoxSet;

    typedef JointModelRevoluteTpl<Scalar, Options, 0> JointModelRX;
    typedef JointModelRevoluteTpl<Scalar, Options, 1> JointModelRY;
    typedef JointModelRevoluteTpl<Scalar, Options, 2> JointModelRZ;
    typedef JointModelRevoluteUnalignedTpl<Scalar, Options> JointModelRevoluteUnaligned;

    typedef JointModelPrismaticTpl<Scalar, Options, 0> JointModelPX;
    typedef JointModelPrismaticTpl<Scalar, Options, 1> JointModelPY;
    typedef JointModelPrismaticTpl<Scalar, Options, 2> JointModelPZ;
    typedef JointModelPrismaticUnalignedTpl<Scalar, Options> JointModelPrismaticUnaligned;

    typedef boost::mpl::vector<
      JointModelRX,
      JointModelRY,
      JointModelRZ,
      JointModelRevoluteUnaligned,
      JointModelPX,
      JointModelPY,
      JointModelPZ,
      JointModelPrismaticUnaligned>
      ValidJointTypes;

    template<template<typename, int> class JointCollectionTpl>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & active_joints)
    {
      init(model, active_joints, model.lowerPositionLimit, model.upperPositionLimit);
    }

    template<
      template<typename, int>
      class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & active_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub)
    {
      init(model, active_joints, lb, ub);
    }

    ConstraintData createData() const
    {
      return ConstraintData(*this);
    }

    int size() const
    {
      return int(row_active_indexes.size());
    }
    template<template<typename, int> class JointCollectionTpl>
    void calc(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const;

    template<template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
    void jacobian(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & _jacobian_matrix) const;

    /// \brief Returns the sparsity associated with a given row
    const BooleanVector & getRowSparsityPattern(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < size());
      return row_sparsity_pattern[size_t(row_id)];
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    const EigenIndexVector & getRowActiveIndexes(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < size());
      return row_active_indexes[size_t(row_id)];
    }

    /// \brief Returns the vector of configuration vector index for active lower bounds
    const EigenIndexVector & getActiveLowerBoundConstraints() const
    {
      return active_lower_bound_constraints;
    }

    /// \brief Returns the vector of tangent vector index for active lower bounds
    const EigenIndexVector & getActiveLowerBoundConstraintsTangent() const
    {
      return active_lower_bound_constraints_tangent;
    }

    /// \brief Returns the vector of configuration vector index for active upper bounds
    const EigenIndexVector & getActiveUpperBoundConstraints() const
    {
      return active_upper_bound_constraints;
    }

    /// \brief Returns the vector of tangent vector index for active upper bounds
    const EigenIndexVector & getActiveUpperBoundConstraintsTangent() const
    {
      return active_upper_bound_constraints_tangent;
    }

    const ConstraintSet & set() const
    {
      return m_set;
    }

    ConstraintSet & set()
    {
      return m_set;
    }

  protected:
    template<
      template<typename, int>
      class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration>
    void init(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & active_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub);

    /// \brief Check whether the active joints are bound to the joint types contained in
    /// SupportedJointTypes.
    template<template<typename, int> class JointCollectionTpl>
    static int check_active_joints(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & active_joints);

    /// \brief Selected dof in the configuration vector
    EigenIndexVector active_configuration_components;
    /// \brief Active lower bound constraints
    EigenIndexVector active_lower_bound_constraints, active_lower_bound_constraints_tangent;
    /// \brief Active upper bound constraints
    EigenIndexVector active_upper_bound_constraints, active_upper_bound_constraints_tangent;
    /// \brief Lower and upper limit values for active constraints
    BoxSet active_configuration_limits;

    VectorOfBooleanVector row_sparsity_pattern;
    VectofOfEigenIndexVector row_active_indexes;

    ConstraintSet m_set;
  };

  template<typename _Scalar, int _Options>
  struct JointLimitConstraintDataTpl
  : ConstraintDataBase<JointLimitConstraintDataTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef ConstraintModelBase<JointLimitConstraintDataTpl> Base;
    typedef std::vector<JointIndex> JointIndexVector;

    typedef JointLimitConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef typename ConstraintModel::VectorXs VectorXs;

    explicit JointLimitConstraintDataTpl(const ConstraintModel & constraint_model)
    : constraint_residual(constraint_model.size())
    {
    }

    bool operator==(const JointLimitConstraintDataTpl & other) const
    {
      if (this == &other)
        return true;
      return this->constraint_residual == other.constraint_residual;
    }

    /// \brief Residual of the constraint
    VectorXs constraint_residual;
  };
} // namespace pinocchio

#include "pinocchio/algorithm/constraints/joint-limit-constraint.hxx"

#endif // ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__
