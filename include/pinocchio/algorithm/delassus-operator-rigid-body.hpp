//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__
#define __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/delassus-operator-base.hpp"
#include "pinocchio/utils/reference.hpp"

#include "pinocchio/algorithm/constraints/constraint-collection-default.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-generic.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-generic.hpp"

namespace pinocchio
{

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    template<typename T> class Holder = std::reference_wrapper>
  struct DelassusOperatorRigidBodySystemsTpl;

  template<
    typename _Scalar,
    int _Options,
    template<typename, int> class JointCollectionTpl,
    class _ConstraintModel,
    template<typename T> class Holder>
  struct traits<DelassusOperatorRigidBodySystemsTpl<
    _Scalar,
    _Options,
    JointCollectionTpl,
    _ConstraintModel,
    Holder>>
  {
    typedef _Scalar Scalar;

    enum
    {
      Options = _Options,
      RowsAtCompileTime = Eigen::Dynamic
    };

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> DenseMatrix;
    typedef DenseMatrix Matrix;

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::Data Data;

    typedef _ConstraintModel ConstraintModel;
    typedef typename ConstraintModel::ConstraintData ConstraintData;

    typedef PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(ConstraintModel) ConstraintModelVector;
    typedef PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(ConstraintData) ConstraintDataVector;

    typedef const Vector & getDampingReturnType;
  };

  template<
    typename _Scalar,
    int _Options,
    template<typename, int> class _JointCollectionTpl,
    class _ConstraintModel,
    template<typename T> class Holder>
  struct DelassusOperatorRigidBodySystemsTpl
  : DelassusOperatorBase<DelassusOperatorRigidBodySystemsTpl<
      _Scalar,
      _Options,
      _JointCollectionTpl,
      _ConstraintModel,
      Holder>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef DelassusOperatorRigidBodySystemsTpl Self;
    typedef DelassusOperatorBase<Self> Base;

    typedef typename traits<Self>::Scalar Scalar;
    enum
    {
      Options = traits<Self>::Options
    };

    typedef typename traits<Self>::Vector Vector;
    typedef typename traits<Self>::DenseMatrix DenseMatrix;

    typedef typename traits<Self>::Model Model;
    typedef Holder<const Model> ModelHolder;
    typedef typename traits<Self>::Data Data;
    typedef Holder<Data> DataHolder;

    typedef typename Data::Force Force;
    typedef typename Data::VectorXs VectorXs;
    typedef PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(Force) ForceVector;

    typedef typename traits<Self>::ConstraintModel ConstraintModel;
    typedef typename traits<Self>::ConstraintModelVector ConstraintModelVector;
    typedef Holder<const ConstraintModelVector> ConstraintModelVectorHolder;

    typedef typename traits<Self>::ConstraintData ConstraintData;
    typedef typename traits<Self>::ConstraintDataVector ConstraintDataVector;
    typedef Holder<ConstraintDataVector> ConstraintDataVectorHolder;

    DelassusOperatorRigidBodySystemsTpl(
      const ModelHolder & model_ref,
      const DataHolder & data_ref,
      const ConstraintModelVectorHolder & constraint_models_ref,
      const ConstraintDataVectorHolder & constraint_datas_ref,
      const Scalar min_damping_value = Scalar(1e-8))
    : Base()
    , m_size(evalConstraintSize(helper::get_ref(constraint_models_ref)))
    , m_model_ref(model_ref)
    , m_data_ref(data_ref)
    , m_constraint_models_ref(constraint_models_ref)
    , m_constraint_datas_ref(constraint_datas_ref)
    , m_custom_data(helper::get_ref(model_ref))
    , m_dirty(true)
    , m_damping(Vector::Zero(m_size))
    , m_compliance(Vector::Zero(m_size))
    , m_sum_compliance_damping(Vector::Zero(m_size))
    , m_sum_compliance_damping_inverse(Vector::Zero(m_size))
    {
      assert(model().check(data()) && "data is not consistent with model.");
      PINOCCHIO_CHECK_ARGUMENT_SIZE(
        constraint_models().size(), constraint_datas().size(),
        "The sizes of contact vector models and contact vector data are not the same.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        min_damping_value >= Scalar(0) && "The damping value should be positive.");

      updateDamping(min_damping_value);
    }

    ///
    /// \brief Update the intermediate computations according to a new configuration vector entry
    ///
    /// \param[in] q Configuration vector
    ///
    template<typename ConfigVectorType>
    void compute(const Eigen::MatrixBase<ConfigVectorType> & q);

    ///
    /// \brief Update the intermediate computations before calling solveInPlace or operator*
    ///
    void compute();

    const Model & model() const
    {
      return helper::get_ref(m_model_ref);
    }

    Data & data()
    {
      return helper::get_ref(m_data_ref);
    }
    const Data & data() const
    {
      return helper::get_ref(m_data_ref);
    }

    const ConstraintModelVector & constraint_models() const
    {
      return helper::get_ref(m_constraint_models_ref);
    }

    const ConstraintDataVector & constraint_datas() const
    {
      return helper::get_ref(m_constraint_datas_ref);
    }
    ConstraintDataVector & constraint_datas()
    {
      return helper::get_ref(m_constraint_datas_ref);
    }

    Eigen::DenseIndex size() const
    {
      return m_size;
    }
    Eigen::DenseIndex rows() const
    {
      return m_size;
    }
    Eigen::DenseIndex cols() const
    {
      return m_size;
    }

    void update(
      const ConstraintModelVectorHolder & constraint_models_ref,
      const ConstraintDataVectorHolder & constraint_datas_ref)
    {
      if (
        helper::get_pointer(m_constraint_models_ref) == helper::get_pointer(constraint_models_ref)
        && helper::get_pointer(m_constraint_datas_ref) == helper::get_pointer(constraint_datas_ref))
        return;
      m_constraint_models_ref = constraint_models_ref;
      m_constraint_datas_ref = constraint_datas_ref;
      m_dirty = true;
    }

    template<typename MatrixIn, typename MatrixOut>
    void applyOnTheRight(
      const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & res) const;

    template<typename VectorLike>
    void updateDamping(const Eigen::MatrixBase<VectorLike> & vec)
    {
      m_damping = vec;
      updateSumComplianceDamping();
    }

    void updateDamping(const Scalar & mu)
    {
      updateDamping(Vector::Constant(size(), mu));
    }

    const Vector & getDamping() const
    {
      return m_damping;
    }

    template<typename VectorLike>
    void updateCompliance(const Eigen::MatrixBase<VectorLike> & compliance_vector)
    {
      m_compliance = compliance_vector;
      updateSumComplianceDamping();
    }

    void updateCompliance(const Scalar & compliance_value)
    {
      updateCompliance(Vector::Constant(size(), compliance_value));
    }

    const Vector & getCompliance() const
    {
      return m_compliance;
    }

    template<typename MatrixLike>
    void solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat) const;

    struct CustomData
    {
      typedef typename Data::Motion Motion;
      typedef typename Data::Force Force;

      typedef typename PINOCCHIO_ALIGNED_STD_VECTOR(Motion) MotionVector;
      typedef typename PINOCCHIO_ALIGNED_STD_VECTOR(Force) ForceVector;

      CustomData(const Model & model)
      : a(size_t(model.njoints), Motion::Zero())
      , oa_augmented(size_t(model.njoints), Motion::Zero())
      , u(model.nv)
      , ddq(model.nv)
      , f(size_t(model.njoints))
      , of_augmented(size_t(model.njoints))
      {
      }

      MotionVector a, oa_augmented;
      VectorXs u, ddq;
      ForceVector f, of_augmented;
    };

    const CustomData & getCustomData() const
    {
      return m_custom_data;
    }

    CustomData & getCustomData()
    {
      return m_custom_data;
    }

    struct AugmentedMassMatrixOperator
    {
      AugmentedMassMatrixOperator(const DelassusOperatorRigidBodySystemsTpl & delassus_operator)
      : m_self(delassus_operator)
      {
      }

      template<typename MatrixLike>
      void solveInPlace(
        const Eigen::MatrixBase<MatrixLike> & mat, bool reset_joint_force_vector = true) const;

    protected:
      const DelassusOperatorRigidBodySystemsTpl & m_self;
    };

    AugmentedMassMatrixOperator getAugmentedMassMatrixOperator() const
    {
      return AugmentedMassMatrixOperator(*this);
    }

  protected:
    static Eigen::DenseIndex evalConstraintSize(const ConstraintModelVector & constraint_models)
    {
      Eigen::DenseIndex size = 0;
      for (const auto & cm : constraint_models)
        size += cm.size();

      return size;
    }

    inline void compute_conclude()
    {
      m_dirty = false;
    }

    void updateSumComplianceDamping()
    {
      m_sum_compliance_damping = m_damping + m_compliance;
      m_sum_compliance_damping_inverse = m_sum_compliance_damping.cwiseInverse();
      m_dirty = true;
    }

    // Holders
    Eigen::DenseIndex m_size;
    ModelHolder m_model_ref;
    DataHolder m_data_ref;
    ConstraintModelVectorHolder m_constraint_models_ref;
    ConstraintDataVectorHolder m_constraint_datas_ref;

    mutable CustomData m_custom_data;
    bool m_dirty;
    Vector m_damping, m_compliance;
    Vector m_sum_compliance_damping, m_sum_compliance_damping_inverse;
  };

} // namespace pinocchio

#include "pinocchio/algorithm/delassus-operator-rigid-body.hxx"

#endif // ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__
