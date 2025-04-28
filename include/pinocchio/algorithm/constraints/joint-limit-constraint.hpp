//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__
#define __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__

#include "pinocchio/math/fwd.hpp"

#include <boost/mpl/vector.hpp>

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/joint-limit-constraint-cone.hpp"
#include "pinocchio/algorithm/constraints/jointwise-constraint-base.hpp"
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
    static constexpr bool has_baumgarte_corrector = true;
    static constexpr bool has_baumgarte_corrector_vector = false;

    typedef JointLimitConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef JointLimitConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef JointLimitConstraintConeTpl<Scalar> ConstraintSet;

    typedef ConstraintModel Model;
    typedef ConstraintData Data;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> JacobianMatrixType;
    typedef VectorXs VectorConstraintSize;

    typedef VectorXs ComplianceVectorType;
    typedef ComplianceVectorType & ComplianceVectorTypeRef;
    typedef const ComplianceVectorType & ComplianceVectorTypeConstRef;

    typedef EigenStorageTpl<VectorXs> EigenStorageVector;
    typedef typename EigenStorageVector::RefMapType ActiveComplianceVectorTypeRef;
    typedef typename EigenStorageVector::ConstRefMapType ActiveComplianceVectorTypeConstRef;

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
  : JointWiseConstraintModelBase<JointLimitConstraintModelTpl<_Scalar, _Options>>
  , ConstraintModelCommonParameters<JointLimitConstraintModelTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef JointLimitConstraintModelTpl Self;
    typedef JointWiseConstraintModelBase<Self> Base;
    typedef ConstraintModelBase<Self> RootBase;

    typedef ConstraintModelCommonParameters<JointLimitConstraintModelTpl> BaseCommonParameters;

    template<typename NewScalar, int NewOptions>
    friend struct JointLimitConstraintModelTpl;

    static const ConstraintFormulationLevel constraint_formulation_level =
      traits<JointLimitConstraintModelTpl>::constraint_formulation_level;
    typedef typename traits<Self>::ComplianceVectorType ComplianceVectorType;
    typedef typename traits<Self>::EigenStorageVector EigenStorageVector;
    typedef typename traits<Self>::ActiveComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef
      typename traits<Self>::ActiveComplianceVectorTypeConstRef ActiveComplianceVectorTypeConstRef;
    typedef
      typename traits<Self>::BaumgarteCorrectorVectorParameters BaumgarteCorrectorVectorParameters;
    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    typedef JointLimitConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef JointLimitConstraintConeTpl<Scalar> ConstraintSet;
    typedef BoxSetTpl<Scalar, Options> BoxSet;

    using RootBase::jacobian;
    using typename Base::BooleanVector;
    using typename Base::EigenIndexVector;

    typedef std::vector<BooleanVector> VectorOfBooleanVector;
    typedef std::vector<EigenIndexVector> VectofOfEigenIndexVector;
    typedef std::vector<size_t> VectorOfSize;
    typedef std::vector<JointIndex> JointIndexVector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef VectorXs VectorConstraintSize;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, MAX_JOINT_NV, Eigen::RowMajor>
      CompactTangentMap_t;

    JointLimitConstraintModelTpl()
    : active_compliance(active_compliance_storage.map())
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
      const JointIndexVector & _activable_joints)
    : active_compliance_storage(0, 1)
    , active_compliance(active_compliance_storage.map())
    {
      init(model, _activable_joints, model.lowerPositionLimit, model.upperPositionLimit);
    }

    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & _activable_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub)
    : active_compliance_storage(0, 1)
    , active_compliance(active_compliance_storage.map())
    {
      init(model, _activable_joints, lb, ub);
    }

    /// \brief Cast operator
    template<typename NewScalar>
    typename CastType<NewScalar, JointLimitConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, JointLimitConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::cast(res);
      BaseCommonParameters::template cast<NewScalar>(res);

      res.activable_joints = activable_joints;
      res.nq_reduce = nq_reduce;
      res.lower_activable_size = lower_activable_size;
      res.lower_active_size = lower_active_size;
      res.row_sparsity_pattern = row_sparsity_pattern;
      res.row_indexes = row_indexes;
      res.bound_position_limit = bound_position_limit.template cast<NewScalar>();
      res.bound_position_margin = bound_position_margin.template cast<NewScalar>();
      res.activable_idx_qs = activable_idx_qs;
      res.active_set_indexes = active_set_indexes;
      res.activable_idx_rows = activable_idx_rows;
      res.activable_idx_qs_reduce = activable_idx_qs_reduce;
      res.activable_nvs = activable_nvs;
      res.activable_idx_vs = activable_idx_vs;
      res.active_idx_rows = active_idx_rows;
      res.active_idx_qs_reduce = active_idx_qs_reduce;
      res.active_nvs = active_nvs;
      res.active_idx_vs = active_idx_vs;
      res.m_set = m_set.template cast<NewScalar>();
      res.active_compliance_storage = active_compliance_storage.template cast<NewScalar>();
      res.active_compliance = res.active_compliance_storage.map();

      return res;
    }

    ConstraintData createData() const
    {
      return ConstraintData(*this);
    }

    int size() const
    {
      return int(activable_idx_row.size());
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

    /// \brief Returns the sparsity associated with a given row
    const BooleanVector & getRowActivableSparsityPattern(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < size());
      return row_sparsity_pattern[activable_idx_rows[static_cast<size_t>(row_id)]];
    }

    /// \brief Returns the sparsity associated with a given row
    /// This vector is computed when calling the calc method.
    const BooleanVector & getRowActiveSparsityPattern(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < activeSize());
      return row_sparsity_pattern[active_idx_rows[static_cast<size_t>(row_id)]];
    }

    /// \brief Returns the vector of the activable indexes associated with a given row
    const EigenIndexVector & getRowActivableIndexes(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < size());
      return row_indexes[activable_idx_rows[static_cast<size_t>(row_id)]];
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    /// This vector is computed when calling the calc method.
    const EigenIndexVector & getRowActiveIndexes(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < activeSize());
      return row_indexes[active_idx_rows[static_cast<size_t>(row_id)]];
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
      
    /// \brief Returns the vector of the active indexes associated with a given row
    /// This vector is computed when calling the calc method.
    const VectorOfSize & getActiveSetIndexes() const
    {
      return active_set_indexes;
    }

    /// Specialized accessor

    const JointIndexVector & getActivableJoints() const
    {
      return activable_joints;
    }
    int getNqReduce() const
    {
      return nq_reduce;
    }
    int lowerSize() const
    {
      return lower_activable_size;
    }
    int lowerActiveSize() const
    {
      return lower_active_size;
    }
    int upperSize() const
    {
      return size() - lowerSize();
    }
    int upperActiveSize() const
    {
      return activeSize() - lowerActiveSize();
    }

    const VectorXs & getBoundPositionLimit() const
    {
      return bound_position_limit;
    }
    const VectorXs & getBoundPositionMargin() const
    {
      return bound_position_margin;
    }

    const EigenIndexVector getActivableIdxQs() const
    {
      return activable_idx_qs;
    }

    const VectorOfSize & getActiveSetIndexes() const
    {
      return active_set_indexes;
    }

    const EigenIndexVector getActivableIdxQsReduce() const
    {
      return activable_idx_qs_reduce;
    }
    const EigenIndexVector getActiveIdxQsReduce() const
    {
      return active_idx_qs_reduce;
    }
    const EigenIndexVector getActivableNvs() const
    {
      return activable_nvs;
    }
    const EigenIndexVector getActiveNvs() const
    {
      return active_nvs;
    }
    const EigenIndexVector getActivableIdxVs() const
    {
      return activable_idx_vs;
    }
    const EigenIndexVector getActiveIdxVs() const
    {
      return active_idx_vs;
    }

    // row_sparsity_pattern, row_indexes, activable_idx_rows, active_idx_rows are
    // not exposed as they only privately allow getRowActiv[e/able]SparsityPattern and
    // getRowActiv[e/able]Indexes

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
             && activable_joints == other.activable_joints && nq_reduce == other.nq_reduce
             && lower_activable_size == other.lower_activable_size
             && lower_active_size == other.lower_active_size
             && row_sparsity_pattern == other.row_sparsity_pattern
             && row_indexes == other.row_indexes
             && bound_position_limit == other.bound_position_limit
             && bound_position_margin == other.bound_position_margin
             && activable_idx_qs == other.activable_idx_qs
             && active_set_indexes == other.active_set_indexes
             && activable_idx_rows == other.activable_idx_rows
             && activable_idx_qs_reduce == other.activable_idx_qs_reduce
             && activable_nvs == other.activable_nvs && activable_idx_vs == other.activable_idx_vs
             && active_idx_rows == other.active_idx_rows
             && active_idx_qs_reduce == other.active_idx_qs_reduce && active_nvs == other.active_nvs
             && active_idx_vs == other.active_idx_vs && m_set == other.m_set
             && active_compliance_storage == other.active_compliance_storage;
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
        activable_joints = other.activable_joints;
        nq_reduce = other.nq_reduce;
        lower_activable_size = other.lower_activable_size;
        lower_active_size = other.lower_active_size;
        row_sparsity_pattern = other.row_sparsity_pattern;
        row_indexes = other.row_indexes;
        bound_position_limit = other.bound_position_limit;
        bound_position_margin = other.bound_position_margin;
        activable_idx_qs = other.activable_idx_qs;
        active_set_indexes = other.active_set_indexes;
        activable_idx_rows = other.activable_idx_rows;
        activable_idx_qs_reduce = other.activable_idx_qs_reduce;
        activable_nvs = other.activable_nvs;
        activable_idx_vs = other.activable_idx_vs;
        active_idx_rows = other.active_idx_rows;
        active_idx_qs_reduce = other.active_idx_qs_reduce;
        active_nvs = other.active_nvs;
        active_idx_vs = other.active_idx_vs;
        m_set = other.m_set;
        active_compliance_storage = other.active_compliance_storage;
        active_compliance = active_compliance_storage.map();
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
      const JointIndexVector & _activable_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub);

    /// @brief List of joints that are concerned by the constraint. size nja
    JointIndexVector activable_joints;
    /// @brief nq size given the considered joints
    /// nq_reduce = SUM(j in activable_joints) j.nq
    int nq_reduce;
    /// @brief number of lower bound limite activable and active
    int lower_activable_size, lower_active_size;

    /// @brief Sparsity pattern for each considered joint. size nja
    VectorOfBooleanVector row_sparsity_pattern;
    VectofOfEigenIndexVector row_indexes;

    /// @brief Limit value of lower and upper bound in the constraint (size size()=lsize+usize)
    VectorXs bound_position_limit;
    /// @brief Margin value of lower and upper bound in the constraint (size size()=lsize+usize)
    VectorXs bound_position_margin;

    /// @brief give for each activable constraint the qs in the configuration vector
    EigenIndexVector activable_idx_qs;

    /// \brief Vector containing the indexes of the constraints in the active set.
    /// the size of the vector is activeSize()
    /// each element have value < size()
    VectorOfSize active_set_indexes;

    /// @brief give for each active/activable constraint the row_id of sparsity pattern
    VectorOfSize activable_idx_rows, active_idx_rows;
    /// @brief give for each active/activable constraint  of sparsity pattern
    EigenIndexVector activable_idx_qs_reduce, active_idx_qs_reduce;
    /// @brief For each dof, the associated nv and idx_v to exploit tangent map sparsity
    EigenIndexVector activable_nvs, active_nvs;
    EigenIndexVector activable_idx_vs, active_idx_vs;

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
    typedef typename ConstraintModel::CompactTangentMap_t CompactTangentMap_t;
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
    , compact_tangent_map(CompactTangentMap_t::Zero(constraint_model.getNqReduce(), MAX_JOINT_NV))
    {
      constraint_residual_storage.resize(0);
    }

    bool operator==(const JointLimitConstraintDataTpl & other) const
    {
      if (this == &other)
        return true;
      return (
        this->constraint_residual == other.constraint_residual
        && this->compact_tangent_map == other.compact_tangent_map)
    }

    bool operator!=(const JointLimitConstraintDataTpl & other) const
    {
      return !(*this == other);
    }

    JointLimitConstraintDataTpl & operator=(const JointLimitConstraintDataTpl & other)
    {
      if (this != &other)
      {
        activable_constraint_residual = other.activable_constraint_residual;
        constraint_residual_storage = other.constraint_residual_storage;
        compact_tangent_map = other.compact_tangent_map;
      }
      return *this;
    }

    /// \brief Residual of all the activable constraints
    VectorXs activable_constraint_residual;

    /// \brief Residual of the active constraints
    EigenStorageVector constraint_residual_storage;
    typename EigenStorageVector::RefMapType constraint_residual;

    /// @brief Compact storage of the tangent map
    CompactTangentMap_t compact_tangent_map;

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
