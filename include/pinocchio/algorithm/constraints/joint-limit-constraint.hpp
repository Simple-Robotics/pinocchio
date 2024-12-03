//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__
#define __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__

#include "pinocchio/math/fwd.hpp"

#include <boost/mpl/vector.hpp>

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/joint-limit-constraint-cone.hpp"
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
    typedef JointLimitConstraintConeTpl<Scalar> ConstraintSet;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef VectorXs VectorConstraintSize;

    typedef VectorXs ComplianceVectorType;
    typedef ComplianceVectorType & ComplianceVectorTypeRef;
    typedef const ComplianceVectorType & ComplianceVectorTypeConstRef;

    template<typename InputMatrix>
    struct JacobianMatrixProductReturnType
    {
      typedef typename InputMatrix::Scalar Scalar;
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(InputMatrix) InputMatrixPlain;
      typedef Eigen::
        Matrix<Scalar, Eigen::Dynamic, InputMatrix::ColsAtCompileTime, InputMatrixPlain::Options>
          type;
    };

    template<typename InputMatrix>
    struct JacobianTransposeMatrixProductReturnType
    {
      typedef typename InputMatrix::Scalar Scalar;
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(InputMatrix) InputMatrixPlain;
      typedef Eigen::Matrix<
        Scalar,
        Eigen::Dynamic,
        InputMatrixPlain::ColsAtCompileTime,
        InputMatrixPlain::Options>
        type;
    };
  };

  template<typename _Scalar, int _Options>
  struct traits<JointLimitConstraintDataTpl<_Scalar, _Options>>
  : traits<JointLimitConstraintModelTpl<_Scalar, _Options>>
  {
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
    typedef VectorXs VectorConstraintSize;
    typedef VectorXs ComplianceVectorType;

    using typename Base::BooleanVector;
    using typename Base::EigenIndexVector;

    typedef std::vector<BooleanVector> VectorOfBooleanVector;
    typedef std::vector<EigenIndexVector> VectofOfEigenIndexVector;

    typedef JointLimitConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef JointLimitConstraintConeTpl<Scalar> ConstraintSet;
    typedef BoxSetTpl<Scalar, Options> BoxSet;

    //    typedef JointModelRevoluteTpl<Scalar, Options, 0> JointModelRX;
    //    typedef JointModelRevoluteTpl<Scalar, Options, 1> JointModelRY;
    //    typedef JointModelRevoluteTpl<Scalar, Options, 2> JointModelRZ;
    //    typedef JointModelRevoluteUnalignedTpl<Scalar, Options> JointModelRevoluteUnaligned;
    //
    //    typedef JointModelPrismaticTpl<Scalar, Options, 0> JointModelPX;
    //    typedef JointModelPrismaticTpl<Scalar, Options, 1> JointModelPY;
    //    typedef JointModelPrismaticTpl<Scalar, Options, 2> JointModelPZ;
    //    typedef JointModelPrismaticUnalignedTpl<Scalar, Options> JointModelPrismaticUnaligned;

    //    typedef boost::mpl::vector<
    //      JointModelRX,
    //      JointModelRY,
    //      JointModelRZ,
    //      JointModelRevoluteUnaligned,
    //      JointModelPX,
    //      JointModelPY,
    //      JointModelPZ,
    //      JointModelPrismaticUnaligned>
    //      ValidJointTypes;

    template<template<typename, int> class JointCollectionTpl>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & active_joints)
    {
      init(model, active_joints, model.lowerPositionLimit, model.upperPositionLimit);
    }

    template<
      template<typename, int> class JointCollectionTpl,
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

    /// \brief Cast operator
    template<typename NewScalar>
    typename CastType<NewScalar, JointLimitConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, JointLimitConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::template cast<NewScalar>(res);

      res.active_configuration_components = active_configuration_components;
      res.active_lower_bound_constraints = active_lower_bound_constraints;
      res.active_lower_bound_constraints_tangent = active_lower_bound_constraints_tangent;
      res.active_upper_bound_constraints = active_upper_bound_constraints;
      res.active_upper_bound_constraints_tangent = active_upper_bound_constraints_tangent;
      res.active_configuration_limits = active_configuration_limits;
      res.row_active_indexes = row_active_indexes;
      res.row_sparsity_pattern = row_sparsity_pattern;
      res.active_configuration_components = active_configuration_components;

      res.m_set = m_set.template cast<NewScalar>();
      res.m_compliance = m_compliance.template cast<NewScalar>();
      return res;
    }

    ConstraintData createData() const
    {
      return ConstraintData(*this);
    }

    int size() const
    {
      return int(row_active_indexes.size());
    }

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
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

    /// \brief Returns the compliance internally stored in the constraint model
    const ComplianceVectorType & compliance() const
    {
      return m_compliance;
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ComplianceVectorType & compliance()
    {
      return m_compliance;
    }

    ///
    ///  \brief Comparison operator
    ///
    /// \param[in] other Other FrictionalJointConstraintModelTpl to compare with.
    ///
    /// \returns true if the two *this is equal to other (type, joint1_id and placement attributs
    /// must be the same).
    ///
    bool operator==(const JointLimitConstraintModelTpl & other) const
    {
      return base() == other.base()
             && active_configuration_components == other.active_configuration_components
             && active_lower_bound_constraints == other.active_lower_bound_constraints
             && active_lower_bound_constraints_tangent
                  == other.active_lower_bound_constraints_tangent
             && active_upper_bound_constraints == other.active_upper_bound_constraints
             && active_upper_bound_constraints_tangent
                  == other.active_upper_bound_constraints_tangent
             && active_configuration_limits == other.active_configuration_limits
             && row_active_indexes == other.row_active_indexes
             && row_sparsity_pattern == other.row_sparsity_pattern && m_set == other.m_set
             && m_compliance == other.m_compliance;
    }

  protected:
    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration>
    void init(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & active_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub);

    //    /// \brief Check whether the active joints are bound to the joint types contained in
    //    /// SupportedJointTypes.
    //    template<template<typename, int> class JointCollectionTpl>
    //    static int check_active_joints(
    //      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    //      const JointIndexVector & active_joints);

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
    ComplianceVectorType m_compliance;
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
