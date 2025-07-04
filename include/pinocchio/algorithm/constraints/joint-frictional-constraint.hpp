//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_frictional_joint_constraint_hpp__
#define __pinocchio_algorithm_constraints_frictional_joint_constraint_hpp__

#include "pinocchio/math/fwd.hpp"

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/box-set.hpp"
#include "pinocchio/algorithm/constraints/jointwise-constraint-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-common-parameters.hpp"

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

    static constexpr ConstraintFormulationLevel constraint_formulation_level =
      ConstraintFormulationLevel::VELOCITY_LEVEL;
    static constexpr bool has_baumgarte_corrector = false;
    static constexpr bool has_baumgarte_corrector_vector = false;

    typedef FrictionalJointConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef FrictionalJointConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef BoxSetTpl<Scalar, Options> ConstraintSet;

    typedef ConstraintModel Model;
    typedef ConstraintData Data;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> JacobianMatrixType;
    typedef VectorXs VectorConstraintSize;

    typedef VectorXs ComplianceVectorType;
    typedef ComplianceVectorType & ComplianceVectorTypeRef;
    typedef const ComplianceVectorType & ComplianceVectorTypeConstRef;

    typedef ComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef ComplianceVectorTypeConstRef ActiveComplianceVectorTypeConstRef;

    typedef Eigen::Matrix<Scalar, 0, 0> BaumgarteVectorType; // empty vector
    typedef BaumgarteCorrectorVectorParametersTpl<BaumgarteVectorType>
      BaumgarteCorrectorVectorParameters;
    typedef BaumgarteCorrectorVectorParameters & BaumgarteCorrectorVectorParametersRef;
    typedef const BaumgarteCorrectorVectorParameters & BaumgarteCorrectorVectorParametersConstRef;

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
  : JointWiseConstraintModelBase<FrictionalJointConstraintModelTpl<_Scalar, _Options>>
  , ConstraintModelCommonParameters<FrictionalJointConstraintModelTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef FrictionalJointConstraintModelTpl Self;
    typedef JointWiseConstraintModelBase<Self> Base;
    typedef ConstraintModelCommonParameters<Self> BaseCommonParameters;
    typedef ConstraintModelBase<FrictionalJointConstraintModelTpl> RootBase;

    template<typename NewScalar, int NewOptions>
    friend struct FrictionalJointConstraintModelTpl;

    static const ConstraintFormulationLevel constraint_formulation_level =
      traits<FrictionalJointConstraintModelTpl>::constraint_formulation_level;
    typedef typename traits<Self>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
    typedef typename traits<Self>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;
    typedef typename traits<Self>::ActiveComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef
      typename traits<Self>::ActiveComplianceVectorTypeConstRef ActiveComplianceVectorTypeConstRef;
    typedef typename traits<Self>::ComplianceVectorType ComplianceVectorType;

    typedef FrictionalJointConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef BoxSetTpl<Scalar, Options> ConstraintSet;

    using RootBase::jacobian;
    using typename Base::BooleanVector;
    using typename Base::EigenIndexVector;

    typedef std::vector<BooleanVector> VectorOfBooleanVector;
    typedef std::vector<EigenIndexVector> VectofOfEigenIndexVector;
    typedef std::vector<JointIndex> JointIndexVector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef VectorXs VectorConstraintSize;

    FrictionalJointConstraintModelTpl()
    {
    }

    template<template<typename, int> class JointCollectionTpl>
    FrictionalJointConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & active_joints)
    : active_joints(active_joints)
    {
      init(model, active_joints);
    }

    /// \brief Cast operator
    template<typename NewScalar>
    typename CastType<NewScalar, FrictionalJointConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, FrictionalJointConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::cast(res);
      BaseCommonParameters::template cast<NewScalar>(res);

      res.active_joints = active_joints;
      res.active_dofs = active_dofs;
      res.row_sparsity_pattern = row_sparsity_pattern;
      res.row_active_indexes = row_active_indexes;

      res.m_set = m_set.template cast<NewScalar>();
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

    int activeSize() const
    {
      return size();
    }

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    BaseCommonParameters & base_common_parameters()
    {
      return static_cast<BaseCommonParameters &>(*this);
    }
    const BaseCommonParameters & base_common_parameters() const
    {
      return static_cast<const BaseCommonParameters &>(*this);
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

    template<typename InputMatrix, template<typename, int> class JointCollectionTpl>
    typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
    jacobianMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      typedef typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
        ReturnType;
      ReturnType res(size(), mat.cols());
      jacobianMatrixProduct(model, data, cdata, mat.derived(), res);
      return res;
    }

    template<
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op = SETTO>
    void jacobianMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot = SetTo()) const;

    template<typename InputMatrix, template<typename, int> class JointCollectionTpl>
    typename traits<Self>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
    jacobianTransposeMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      typedef
        typename traits<Self>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
          ReturnType;
      ReturnType res(model.nv, mat.cols());
      jacobianTransposeMatrixProduct(model, data, cdata, mat.derived(), res);
      return res;
    }

    template<
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op = SETTO>
    void jacobianTransposeMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot = SetTo()) const;

    /// \brief Returns the sparsity associated with a given row
    const BooleanVector & getRowActivableSparsityPattern(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < size());
      return row_sparsity_pattern[size_t(row_id)];
    }

    /// \brief Returns the sparsity associated with a given row
    const BooleanVector & getRowActiveSparsityPattern(const Eigen::DenseIndex row_id) const
    {
      return getRowActivableSparsityPattern(row_id);
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    const EigenIndexVector & getRowActivableIndexes(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < size());
      return row_active_indexes[size_t(row_id)];
    }
    // row_active_indexes[size_t(row_id)]

    /// \brief Returns the vector of the active indexes associated with a given row
    const EigenIndexVector & getRowActiveIndexes(const Eigen::DenseIndex row_id) const
    {
      return getRowActivableIndexes(row_id);
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ActiveComplianceVectorTypeConstRef getActiveCompliance_impl() const
    {
      return this->compliance();
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ActiveComplianceVectorTypeRef getActiveCompliance_impl()
    {
      return this->compliance();
    }

    /// \brief Returns the vector of active joints
    const JointIndexVector & getActiveJoints() const
    {
      return active_joints;
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

    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorNLike,
      ReferenceFrame rf>
    void appendCouplingConstraintInertias(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
      const ReferenceFrameTag<rf> reference_frame) const;

    /// \copydoc Base::mapConstraintForcesToJointTorques
    template<
      template<typename, int> class JointCollectionTpl,
      typename ConstraintForcesLike,
      typename JointTorquesLike>
    void mapConstraintForceToJointTorques(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ConstraintForcesLike> & constraint_forces,
      const Eigen::MatrixBase<JointTorquesLike> & joint_torques) const;

    /// \copydoc Base::mapJointMotionsToConstraintMotions
    template<
      template<typename, int> class JointCollectionTpl,
      typename JointMotionsLike,
      typename ConstraintMotionsLike>
    void mapJointMotionsToConstraintMotion(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<JointMotionsLike> & joint_motions,
      const Eigen::MatrixBase<ConstraintMotionsLike> & constraint_motions) const;

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
      if (this == &other)
        return true;
      return base() == other.base() && base_common_parameters() == other.base_common_parameters()
             && active_dofs == other.active_dofs
             && row_sparsity_pattern == other.row_sparsity_pattern
             && row_active_indexes == other.row_active_indexes && m_set == other.m_set;
    }

    static std::string classname()
    {
      return std::string("FrictionalJointConstraintModel");
    }
    std::string shortname() const
    {
      return classname();
    }

    bool operator!=(const FrictionalJointConstraintModelTpl & other) const
    {
      return !(*this == other);
    }

  protected:
    template<template<typename, int> class JointCollectionTpl>
    void init(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & active_joints);

    JointIndexVector active_joints;
    EigenIndexVector active_dofs;
    VectorOfBooleanVector row_sparsity_pattern;
    VectofOfEigenIndexVector row_active_indexes;

    ConstraintSet m_set;
    using BaseCommonParameters::m_compliance;
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
    typedef ConstraintDataBase<FrictionalJointConstraintDataTpl> Base;
    typedef std::vector<JointIndex> JointIndexVector;

    typedef FrictionalJointConstraintModelTpl<Scalar, Options> ConstraintModel;

    FrictionalJointConstraintDataTpl()
    {
    }

    explicit FrictionalJointConstraintDataTpl(const ConstraintModel & /*constraint_model*/)
    {
    }

    bool operator==(const FrictionalJointConstraintDataTpl & /*other*/) const
    {
      return true;
    }

    bool operator!=(const FrictionalJointConstraintDataTpl & other) const
    {
      return !(*this == other);
    }

    static std::string classname()
    {
      return std::string("FrictionalJointConstraintData");
    }
    std::string shortname() const
    {
      return classname();
    }
  };
} // namespace pinocchio

#include "pinocchio/algorithm/constraints/joint-frictional-constraint.hxx"

#endif // ifndef __pinocchio_algorithm_constraints_frictional_joint_constraint_hpp__
