//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_frictional_joint_constraint_hpp__
#define __pinocchio_algorithm_constraints_frictional_joint_constraint_hpp__

#include "pinocchio/math/fwd.hpp"

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/box-set.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"

namespace pinocchio
{

  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, FrictionalJointConstraintModelTpl<Scalar, Options>>
  {
    typedef FrictionalJointConstraintModelTpl<NewScalar, Options> type;
  };

  template<typename _Scalar, int _Options>
  struct traits<FrictionalJointConstraintModelTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef FrictionalJointConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef FrictionalJointConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef BoxSetTpl<Scalar, Options> ConstraintSet;

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
  struct traits<FrictionalJointConstraintDataTpl<_Scalar, _Options>>
  : traits<FrictionalJointConstraintModelTpl<_Scalar, _Options>>
  {
  };

  template<typename _Scalar, int _Options>
  struct FrictionalJointConstraintModelTpl
  : ConstraintModelBase<FrictionalJointConstraintModelTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef ConstraintModelBase<FrictionalJointConstraintModelTpl> Base;
    typedef std::vector<JointIndex> JointIndexVector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef VectorXs VectorConstraintSize;
    typedef VectorXs ComplianceVectorType;

    using typename Base::BooleanVector;
    using typename Base::EigenIndexVector;

    typedef std::vector<BooleanVector> VectorOfBooleanVector;
    typedef std::vector<EigenIndexVector> VectofOfEigenIndexVector;

    typedef FrictionalJointConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef BoxSetTpl<Scalar, Options> ConstraintSet;

    template<template<typename, int> class JointCollectionTpl>
    FrictionalJointConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & active_joints)
    {
      init(model, active_joints);
      m_set = ConstraintSet(size());
      m_compliance = ComplianceVectorType::Zero(size());
    }

    /// \brief Cast operator
    template<typename NewScalar>
    typename CastType<NewScalar, FrictionalJointConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, FrictionalJointConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::template cast<NewScalar>(res);

      res.active_dofs = active_dofs;
      res.row_sparsity_pattern = row_sparsity_pattern;
      res.row_active_indexes = row_active_indexes;

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
      return int(active_dofs.size());
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
      ConstraintData & cdata) const
    {
      PINOCCHIO_UNUSED_VARIABLE(model);
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);
    }

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

    /// \brief Returns the vector of active rows
    const EigenIndexVector & getActiveDofs() const
    {
      return active_dofs;
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
    bool operator==(const FrictionalJointConstraintModelTpl & other) const
    {
      return base() == other.base() && active_dofs == other.active_dofs
             && row_sparsity_pattern == other.row_sparsity_pattern
             && row_active_indexes == other.row_active_indexes && m_set == other.m_set
             && m_compliance == other.m_compliance;
    }

  protected:
    template<template<typename, int> class JointCollectionTpl>
    void init(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & active_joints);

    EigenIndexVector active_dofs;
    VectorOfBooleanVector row_sparsity_pattern;
    VectofOfEigenIndexVector row_active_indexes;

    ConstraintSet m_set;
    ComplianceVectorType m_compliance;
  };

  template<typename _Scalar, int _Options>
  struct FrictionalJointConstraintDataTpl
  : ConstraintDataBase<FrictionalJointConstraintDataTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef ConstraintModelBase<FrictionalJointConstraintDataTpl> Base;
    typedef std::vector<JointIndex> JointIndexVector;

    typedef FrictionalJointConstraintModelTpl<Scalar, Options> ConstraintModel;

    explicit FrictionalJointConstraintDataTpl(const ConstraintModel & /*constraint_model*/)
    {
    }

    bool operator==(const FrictionalJointConstraintDataTpl & /*other*/) const
    {
      return true;
    }
  };
} // namespace pinocchio

#include "pinocchio/algorithm/constraints/frictional-joint-constraint.hxx"

#endif // ifndef __pinocchio_algorithm_constraints_frictional_joint_constraint_hpp__
