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
#include "pinocchio/utils/template-template-parameter.hpp"

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
    typedef typename helper::remove_ref<ConstraintModel>::type InnerConstraintModel;

    typedef typename helper::remove_ref<ConstraintModel>::type::ConstraintData InnerConstraintData;
    typedef typename std::conditional<
      helper::is_type_holder<ConstraintModel>::value,
      typename internal::extract_template_template_parameter<ConstraintModel>::template type<
        InnerConstraintData>,
      InnerConstraintData>::type ConstraintData;

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
    typedef typename traits<Self>::InnerConstraintModel InnerConstraintModel;
    typedef typename traits<Self>::ConstraintModelVector ConstraintModelVector;
    typedef Holder<const ConstraintModelVector> ConstraintModelVectorHolder;

    typedef typename traits<Self>::ConstraintData ConstraintData;
    typedef typename traits<Self>::InnerConstraintData InnerConstraintData;
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
        "The sizes of contact vector models and contact vector datas are not the same.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        min_damping_value >= Scalar(0) && "The damping value should be positive.");

      updateDamping(min_damping_value);
      update(constraint_models_ref, constraint_datas_ref);
    }

    /// \brief Update the constraint model and data vectors.
    ///
    /// \param[in] constraint_models_ref Vector of constraint models
    /// \param[in] constraint_datas_ref Vector of constraint datas
    ///
    void update(
      const ConstraintModelVectorHolder & constraint_models_ref,
      const ConstraintDataVectorHolder & constraint_datas_ref);

    ///
    /// \brief Update the intermediate computations according to a new configuration vector entry
    ///
    /// \param[in] q Configuration vector
    ///
    template<typename ConfigVectorType>
    void compute(
      const Eigen::MatrixBase<ConfigVectorType> & q,
      bool apply_on_the_right = true,
      bool solve_in_place = true);

    DenseMatrix matrix(bool enforce_symmetry = false) const
    {
      DenseMatrix res(this->size(), this->size());

      typedef Eigen::Map<VectorXs> MapVectorXs;
      MapVectorXs x = MapVectorXs(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, this->size(), 1));

      for (Eigen::DenseIndex i = 0; i < this->size(); ++i)
      {
        x = VectorXs::Unit(this->size(), i);
        this->applyOnTheRight(x, res.col(i));
      }
      if (enforce_symmetry)
      {
        res = 0.5 * (res + res.transpose());
      }
      return res;
    }

  protected:
    void compute_or_update_decomposition(bool apply_on_the_right, bool solve_in_place);

    ///
    /// \brief Update the internal factorization because the damping or compliance vectors have been
    /// modified
    ///
    void updateDecomposition()
    {
      compute_or_update_decomposition(false, true);
    }

  public:
    ///
    /// \brief Update the intermediate computations before calling solveInPlace or operator*
    ///
    /// \param[in] apply_on_the_right If true, this will update the quantities related to the
    /// applyOnTheRight method
    /// \param[in] solve_in_place If true, this will update the quantities related to the
    /// solveInPlace method
    ///
    /// \remarks By activating or deactivating apply_on_the_right and solve_in_place, this enables
    /// to lower the quantities updated to the minimum, helping to save time overall.
    void compute(bool apply_on_the_right = true, bool solve_in_place = true)
    {
      const ConstraintModelVector & constraint_models_ref = constraint_models();
      ConstraintDataVector & constraint_datas_ref = constraint_datas();

      for (size_t ee_id = 0; ee_id < constraint_models_ref.size(); ++ee_id)
      {
        const auto & cmodel = helper::get_ref<ConstraintModel>(constraint_models_ref[ee_id]);
        auto & cdata = helper::get_ref<ConstraintData>(constraint_datas_ref[ee_id]);
        cmodel.calc(model(), data(), cdata);
      }

      compute_or_update_decomposition(apply_on_the_right, solve_in_place);
    }

  public:
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

    bool isDirty() const
    {
      return m_dirty;
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

    /// \brief solveInPlace operation returning the results of the inverse of the Delassus operator
    /// times the input matrix mat
    ///
    /// \param[in,out] mat Input/output argument containing the right hand side and the result of
    /// the operation
    ///
    /// \warning The parameter is only marked 'const' to make the C++ compiler accept a temporary
    /// expression here. This function will const_cast it, so constness isn't honored here.
    ///
    /// \remarks Even if the method is marked 'const', it will update the internal decomposition if
    /// the Delassus operator is dirty after an update of the damping or compliance values.
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
      for (const ConstraintModel & cm : constraint_models)
      {
        const InnerConstraintModel & cmodel = helper::get_ref<ConstraintModel>(cm);
        size += cmodel.size();
      }

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

    DelassusOperatorRigidBodySystemsTpl & self_const_cast() const
    {
      return const_cast<DelassusOperatorRigidBodySystemsTpl &>(*this);
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
