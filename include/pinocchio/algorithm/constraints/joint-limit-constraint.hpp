//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__
#define __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__

#include "pinocchio/math/fwd.hpp"

#include <boost/mpl/vector.hpp>

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/joint-limit-constraint-cone.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-common-parameters.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-vector-parameters.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hpp"
#include "pinocchio/container/storage.hpp"

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

    static constexpr ConstraintFormulationLevel constraint_formulation_level =
      ConstraintFormulationLevel::POSITION_LEVEL;

    typedef JointLimitConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef JointLimitConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef JointLimitConstraintConeTpl<Scalar> ConstraintSet;

    typedef ConstraintModel Model;
    typedef ConstraintData Data;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef VectorXs VectorConstraintSize;

    typedef VectorXs ComplianceVectorType;
    typedef ComplianceVectorType & ComplianceVectorTypeRef;
    typedef const ComplianceVectorType & ComplianceVectorTypeConstRef;

    typedef EigenStorageTpl<VectorXs> EigenStorageVector;
    typedef typename EigenStorageVector::RefMapType ActiveComplianceVectorTypeRef;
    typedef typename EigenStorageVector::ConstRefMapType ActiveComplianceVectorTypeConstRef;
    // typedef Eigen::Ref<ComplianceVectorType>  ActiveComplianceVectorTypeRef;
    // typedef Eigen::Ref<const ComplianceVectorType> ActiveComplianceVectorTypeConstRef;

    static constexpr bool has_baumgarte_corrector = true;
    typedef VectorXs BaumgarteVectorType;
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
  struct traits<JointLimitConstraintDataTpl<_Scalar, _Options>>
  : traits<JointLimitConstraintModelTpl<_Scalar, _Options>>
  {
  };

  template<typename _Scalar, int _Options>
  struct JointLimitConstraintModelTpl
  : ConstraintModelBase<JointLimitConstraintModelTpl<_Scalar, _Options>>
  , ConstraintModelCommonParameters<JointLimitConstraintModelTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    typedef JointLimitConstraintModelTpl Self;
    enum
    {
      Options = _Options
    };

    typedef ConstraintModelBase<JointLimitConstraintModelTpl> Base;
    typedef ConstraintModelCommonParameters<JointLimitConstraintModelTpl> BaseCommonParameters;

    typedef std::vector<JointIndex> JointIndexVector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef typename traits<Self>::EigenStorageVector EigenStorageVector;
    typedef VectorXs VectorConstraintSize;
    typedef VectorXs MarginVectorType;
    typedef typename traits<Self>::ComplianceVectorType ComplianceVectorType;
    typedef typename traits<Self>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
    typedef typename traits<Self>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;
    typedef typename traits<Self>::ActiveComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef
      typename traits<Self>::ActiveComplianceVectorTypeConstRef ActiveComplianceVectorTypeConstRef;
    // typedef
    //   typename traits<Self>::BaumgarteCorrectorVectorParameters
    //   BaumgarteCorrectorVectorParameters;
    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    static const ConstraintFormulationLevel constraint_formulation_level =
      traits<JointLimitConstraintModelTpl>::constraint_formulation_level;

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

    JointLimitConstraintModelTpl()
    : active_compliance_storage(size(), 1)
    , active_compliance(active_compliance_storage.map())
    {
    }

    JointLimitConstraintModelTpl(const JointLimitConstraintModelTpl & other)
    : active_compliance(active_compliance_storage.map())
    {
      *this = other;
    }

    template<template<typename, int> class JointCollectionTpl>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & activable_joints)
    : active_compliance_storage(size(), 1)
    , active_compliance(active_compliance_storage.map())
    {
      init(model, activable_joints, model.lowerPositionLimit, model.upperPositionLimit);
    }

    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & activable_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub)
    : active_compliance_storage(size(), 1)
    , active_compliance(active_compliance_storage.map())
    {
      init(model, activable_joints, lb, ub);
    }

    template<typename NewScalar, int NewOptions>
    friend struct JointLimitConstraintModelTpl;

    /// \brief Cast operator
    template<typename NewScalar>
    typename CastType<NewScalar, JointLimitConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, JointLimitConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::cast(res);
      BaseCommonParameters::template cast<NewScalar>(res);

      res.activable_configuration_components = activable_configuration_components;
      res.activable_lower_bound_constraints = activable_lower_bound_constraints;
      res.activable_lower_bound_constraints_tangent = activable_lower_bound_constraints_tangent;
      res.activable_upper_bound_constraints = activable_upper_bound_constraints;
      res.activable_upper_bound_constraints_tangent = activable_upper_bound_constraints_tangent;
      res.activable_configuration_limits =
        activable_configuration_limits.template cast<NewScalar>();
      res.row_activable_indexes = row_activable_indexes;
      res.row_activable_sparsity_pattern = row_activable_sparsity_pattern;

      res.active_set_indexes = active_set_indexes;
      res.active_lower_bound_constraints = active_lower_bound_constraints;
      res.active_lower_bound_constraints_tangent = active_lower_bound_constraints_tangent;
      res.active_upper_bound_constraints = active_upper_bound_constraints;
      res.active_upper_bound_constraints_tangent = active_upper_bound_constraints_tangent;
      res.active_compliance_storage = active_compliance_storage.template cast<NewScalar>();
      res.active_compliance = res.active_compliance_storage.map();
      res.m_set = m_set.template cast<NewScalar>();
      return res;
    }

    ConstraintData createData() const
    {
      return ConstraintData(*this);
    }

    int size() const
    {
      return int(row_activable_indexes.size());
    }

    int activeSize() const
    {
      return int(m_set.size());
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

    /// \brief Resize the constraint by computing the current active set.
    template<template<typename, int> class JointCollectionTpl>
    void resize(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata);

    /// \brief Compute the constraint residual for the active constraints. It assumes that the
    /// active set has been computed by calling `resize` beforehand.
    template<template<typename, int> class JointCollectionTpl>
    void calc(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const;

    /// \brief Returns the vector of the active indexes associated with a given row
    /// This vector is computed when calling the calc method.
    const std::vector<std::size_t> & getActiveSetIndexes() const
    {
      return active_set_indexes;
    }

    /// \brief Returns the sparsity associated with a given row
    const BooleanVector & getRowActivableSparsityPattern(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < size());
      return row_activable_sparsity_pattern[std::size_t(row_id)];
    }

    /// \brief Returns the sparsity associated with a given row
    /// This vector is computed when calling the calc method.
    const BooleanVector & getRowActiveSparsityPattern(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < activeSize());
      return row_activable_sparsity_pattern[active_set_indexes[std::size_t(row_id)]];
    }

    /// \brief Returns the vector of the activable indexes associated with a given row
    const EigenIndexVector & getRowActivableIndexes(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < size());
      return row_activable_indexes[size_t(row_id)];
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    /// This vector is computed when calling the calc method.
    const EigenIndexVector & getRowActiveIndexes(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < activeSize());
      return row_activable_indexes[active_set_indexes[std::size_t(row_id)]];
    }

    /// \brief Returns the vector of configuration vector index for activable lower bounds
    const EigenIndexVector & getActivableLowerBoundConstraints() const
    {
      return activable_lower_bound_constraints;
    }

    /// \brief Returns the vector of configuration vector index for active lower bounds
    /// This vector is computed when calling the calc method.
    const EigenIndexVector & getActiveLowerBoundConstraints() const
    {
      return active_lower_bound_constraints;
    }

    /// \brief Returns the vector of tangent vector index for activable lower bounds
    const EigenIndexVector & getActivableLowerBoundConstraintsTangent() const
    {
      return activable_lower_bound_constraints_tangent;
    }

    /// \brief Returns the vector of tangent vector index for active lower bounds
    /// This vector is computed when calling the calc method.
    const EigenIndexVector & getActiveLowerBoundConstraintsTangent() const
    {
      return active_lower_bound_constraints_tangent;
    }

    /// \brief Returns the vector of configuration vector index for activable upper bounds
    const EigenIndexVector & getActivableUpperBoundConstraints() const
    {
      return activable_upper_bound_constraints;
    }

    /// \brief Returns the vector of configuration vector index for active upper bounds
    /// This vector is computed when calling the calc method.
    const EigenIndexVector & getActiveUpperBoundConstraints() const
    {
      return active_upper_bound_constraints;
    }

    /// \brief Returns the vector of tangent vector index for activable upper bounds
    const EigenIndexVector & getActivableUpperBoundConstraintsTangent() const
    {
      return activable_upper_bound_constraints_tangent;
    }

    /// \brief Returns the vector of tangent vector index for active upper bounds
    /// This vector is computed when calling the calc method.
    const EigenIndexVector & getActiveUpperBoundConstraintsTangent() const
    {
      return active_upper_bound_constraints_tangent;
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ActiveComplianceVectorTypeConstRef getActiveCompliance_impl() const
    {
      return active_compliance;
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ActiveComplianceVectorTypeRef getActiveCompliance_impl()
    {
      return active_compliance;
    }

    const ConstraintSet & set() const
    {
      return m_set;
    }

    ConstraintSet & set()
    {
      return m_set;
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
      return base() == other.base() && base_common_parameters() == other.base_common_parameters()
             && activable_configuration_components == other.activable_configuration_components
             && activable_lower_bound_constraints == other.activable_lower_bound_constraints
             && activable_lower_bound_constraints_tangent
                  == other.activable_lower_bound_constraints_tangent
             && activable_upper_bound_constraints == other.activable_upper_bound_constraints
             && activable_upper_bound_constraints_tangent
                  == other.activable_upper_bound_constraints_tangent
             && activable_configuration_limits == other.activable_configuration_limits
             && row_activable_indexes == other.row_activable_indexes
             && row_activable_sparsity_pattern == other.row_activable_sparsity_pattern
             && active_lower_bound_constraints == other.active_lower_bound_constraints
             && active_lower_bound_constraints_tangent
                  == other.active_lower_bound_constraints_tangent
             && active_upper_bound_constraints == other.active_upper_bound_constraints
             && active_upper_bound_constraints_tangent
                  == other.active_upper_bound_constraints_tangent
             && active_set_indexes == other.active_set_indexes
             && active_compliance == other.active_compliance && m_compliance == other.m_compliance
             && m_set == other.m_set;
    }

    bool operator!=(const JointLimitConstraintModelTpl & other) const
    {
      return !(*this == other);
    }

    JointLimitConstraintModelTpl & operator=(const JointLimitConstraintModelTpl & other)
    {
      if (this != &other)
      {
        base_common_parameters() = other.base_common_parameters();
        activable_configuration_components = other.activable_configuration_components;
        activable_lower_bound_constraints = other.activable_lower_bound_constraints;
        activable_lower_bound_constraints_tangent = other.activable_lower_bound_constraints_tangent;
        activable_upper_bound_constraints = other.activable_upper_bound_constraints;
        activable_upper_bound_constraints_tangent = other.activable_upper_bound_constraints_tangent;
        activable_configuration_limits = other.activable_configuration_limits;
        row_activable_indexes = other.row_activable_indexes;
        row_activable_sparsity_pattern = other.row_activable_sparsity_pattern;
        active_lower_bound_constraints = other.active_lower_bound_constraints;
        active_lower_bound_constraints_tangent = other.active_lower_bound_constraints_tangent;
        active_upper_bound_constraints = other.active_upper_bound_constraints;
        active_upper_bound_constraints_tangent = other.active_upper_bound_constraints_tangent;
        active_set_indexes = other.active_set_indexes;
        m_set = other.m_set;
      }
      return *this;
    }

    // Jacobian operations

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

  public:
    static std::string classname()
    {
      return std::string("JointLimitConstraintModel");
    }
    std::string shortname() const
    {
      return classname();
    }

  protected:
    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration>
    void init(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & activable_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub);

    //    /// \brief Check whether the active joints are bound to the joint types contained in
    //    /// SupportedJointTypes.
    //    template<template<typename, int> class JointCollectionTpl>
    //    static int check_activable_joints(
    //      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    //      const JointIndexVector & activable_joints);

    /// \brief Selected dof in the configuration vector
    EigenIndexVector activable_configuration_components;
    /// \brief Active lower bound constraints
    EigenIndexVector activable_lower_bound_constraints, activable_lower_bound_constraints_tangent;
    /// \brief Active upper bound constraints
    EigenIndexVector activable_upper_bound_constraints, activable_upper_bound_constraints_tangent;
    /// \brief Lower and upper limit values for active constraints
    // TODO: this should be removed (not used anymore)
    BoxSet activable_configuration_limits;

    /// \brief Vector containing the indexes of the constraints in the active set.
    std::vector<std::size_t> active_set_indexes;

    /// \brief Active lower bound constraints that are active in the current configuration
    /// These vectors are computed during the call to the calc method.
    EigenIndexVector active_lower_bound_constraints, active_lower_bound_constraints_tangent;
    /// \brief Active upper bound constraints that are active in the current configuration
    /// These vectors are computed during the call to the calc method.
    EigenIndexVector active_upper_bound_constraints, active_upper_bound_constraints_tangent;

    VectorOfBooleanVector row_activable_sparsity_pattern;
    VectofOfEigenIndexVector row_activable_indexes;

    ConstraintSet m_set;
    using BaseCommonParameters::m_baumgarte_parameters;
    using BaseCommonParameters::m_compliance;

    /// \brief Compliance of the active constraints
    EigenStorageVector active_compliance_storage;
    typename EigenStorageVector::RefMapType active_compliance;
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
    typedef ConstraintDataBase<JointLimitConstraintDataTpl> Base;
    typedef std::vector<JointIndex> JointIndexVector;

    typedef JointLimitConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef typename ConstraintModel::VectorXs VectorXs;

    typedef typename ConstraintModel::EigenStorageVector EigenStorageVector;
    typedef typename ConstraintModel::Base::BooleanVector BooleanVector;
    typedef typename ConstraintModel::Base::EigenIndexVector EigenIndexVector;

    typedef std::vector<BooleanVector> VectorOfBooleanVector;
    typedef std::vector<EigenIndexVector> VectofOfEigenIndexVector;

    JointLimitConstraintDataTpl()
    : constraint_residual(constraint_residual_storage.map())
    {
    }

    JointLimitConstraintDataTpl(const JointLimitConstraintDataTpl & other)
    : constraint_residual(constraint_residual_storage.map())
    {
      *this = other;
    }

    explicit JointLimitConstraintDataTpl(const ConstraintModel & constraint_model)
    : activable_constraint_residual(constraint_model.size())
    , constraint_residual_storage(constraint_model.size(), 1)
    , constraint_residual(constraint_residual_storage.map())
    {
      constraint_residual_storage.resize(0);
    }

    bool operator==(const JointLimitConstraintDataTpl & other) const
    {
      if (this == &other)
        return true;
      return this->constraint_residual == other.constraint_residual;
    }

    bool operator!=(const JointLimitConstraintDataTpl & other) const
    {
      return !(*this == other);
    }

    JointLimitConstraintDataTpl & operator=(const JointLimitConstraintDataTpl & other)
    {
      if (this != &other)
      {
        constraint_residual_storage = other.constraint_residual_storage;
        activable_constraint_residual = other.activable_constraint_residual;
      }
      return *this;
    }

    /// \brief Residual of all the activable constraints
    VectorXs activable_constraint_residual;

    /// \brief Residual of the active constraints
    EigenStorageVector constraint_residual_storage;
    typename EigenStorageVector::RefMapType constraint_residual;

    static std::string classname()
    {
      return std::string("JointLimitConstraintData");
    }
    std::string shortname() const
    {
      return classname();
    }
  };
} // namespace pinocchio

#include "pinocchio/algorithm/constraints/joint-limit-constraint.hxx"

#endif // ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__
