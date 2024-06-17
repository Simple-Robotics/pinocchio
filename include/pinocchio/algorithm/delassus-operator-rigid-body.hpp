//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__
#define __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/delassus-operator-base.hpp"
#include "pinocchio/utils/reference.hpp"

namespace pinocchio {

  template<typename _Scalar, int _Options, template<typename,int> class JointCollectionTpl, template<typename T> class Holder = std::reference_wrapper>
  struct DelassusOperatorRigidBodyTpl;
  
  template<typename _Scalar, int _Options, template<typename,int> class JointCollectionTpl, template<typename T> class Holder>
  struct traits<DelassusOperatorRigidBodyTpl<_Scalar,_Options,JointCollectionTpl,Holder> >
  {
    typedef _Scalar Scalar;
    
    enum { Options = _Options, RowsAtCompileTime = Eigen::Dynamic };
    
    typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1,Options> Vector;
    typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Options> DenseMatrix;
    typedef DenseMatrix Matrix;

    typedef ModelTpl<Scalar,Options,JointCollectionTpl> Model;
    typedef typename Model::Data Data;

    typedef RigidConstraintModelTpl<Scalar,Options> ConstraintModel;
    typedef typename ConstraintModel::ConstraintData ConstraintData;
    
    typedef PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(ConstraintModel) ConstraintModelVector;
    typedef PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(ConstraintData) ConstraintDataVector;
  };
  
  template<typename _Scalar, int _Options, template<typename,int> class JointCollectionTpl, template<typename T> class Holder>
  struct DelassusOperatorRigidBodyTpl
  : DelassusOperatorBase< DelassusOperatorRigidBodyTpl<_Scalar,_Options,JointCollectionTpl,Holder> >
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    typedef DelassusOperatorRigidBodyTpl Self;
    typedef DelassusOperatorBase<Self> Base;
    
    typedef typename traits<Self>::Scalar Scalar;
    enum { Options = traits<Self>::Options };
    
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

    DelassusOperatorRigidBodyTpl(const ModelHolder &model_ref,
                                 const DataHolder &data_ref,
                                 const ConstraintModelVectorHolder &constraint_models_ref,
                                 const ConstraintDataVectorHolder &constraint_datas_ref,
                                 const Scalar min_damping_value = Scalar(1e-8))
    : Base()
    , m_size(evalConstraintSize(helper::get_ref(constraint_models_ref)))
    , m_model_ref(model_ref)
    , m_data_ref(data_ref)
    , m_constraint_models_ref(constraint_models_ref)
    , m_constraint_datas_ref(constraint_datas_ref)
    , m_custom_data(helper::get_ref(model_ref), helper::get_ref(data_ref),evalConstraintSize(helper::get_ref(constraint_models_ref)))
    , m_dirty(true)
    , m_damping(Vector::Constant(m_size,min_damping_value))
    , m_damping_inverse(Vector::Constant(m_size,Scalar(1)/min_damping_value))
    {
      assert(model().check(data()) && "data is not consistent with model.");
      PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_models().size(), constraint_datas().size(),
                                    "The sizes of contact vector models and contact vector data are not the same.");
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
    
    const Model & model() const { return helper::get_ref(m_model_ref); }
    
    Data & data() { return helper::get_ref(m_data_ref); }
    const Data & data() const { return helper::get_ref(m_data_ref); }
    
    const ConstraintModelVector & constraint_models() const { return helper::get_ref(m_constraint_models_ref); }
    
    const ConstraintDataVector & constraint_datas() const { return helper::get_ref(m_constraint_datas_ref); }
    ConstraintDataVector & constraint_datas() { return helper::get_ref(m_constraint_datas_ref); }
    
    Eigen::DenseIndex size() const { return m_size; }
    Eigen::DenseIndex rows() const { return m_size; }
    Eigen::DenseIndex cols() const { return m_size; }
    
    void update(const ConstraintModelVectorHolder &constraint_models_ref,
                const ConstraintDataVectorHolder &constraint_datas_ref)
    {
      if(helper::get_pointer(m_constraint_models_ref) == helper::get_pointer(constraint_models_ref)
         && helper::get_pointer(m_constraint_datas_ref) == helper::get_pointer(constraint_datas_ref))
        return;
      m_constraint_models_ref = constraint_models_ref;
      m_constraint_datas_ref = constraint_datas_ref;
      m_dirty = true;
    }
    
    template<typename MatrixIn, typename MatrixOut>
    void applyOnTheRight(const Eigen::MatrixBase<MatrixIn> & x,
                         const Eigen::MatrixBase<MatrixOut> & res) const;
    
    template<typename VectorLike>
    void updateDamping(const Eigen::MatrixBase<VectorLike> & vec)
    {
      m_damping = vec;
      m_damping_inverse = m_damping.cwiseInverse();
//      mat_tmp = delassus_matrix;
//      mat_tmp += vec.asDiagonal();
//      llt.compute(mat_tmp);
    }
    
    void updateDamping(const Scalar & mu)
    {
      updateDamping(Vector::Constant(size(),mu));
    }
    
    const Vector & getDamping() const { return m_damping; }
    
    template<typename MatrixLike>
    void solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat) const;

    struct CustomData
    {
      typedef typename Data::SE3 SE3;
      typedef typename Data::Inertia Inertia;
      typedef typename Data::Motion Motion;
      typedef typename Data::Matrix6 Matrix6;
      typedef typename Data::Force Force;
      
      typedef typename PINOCCHIO_ALIGNED_STD_VECTOR(SE3) SE3Vector;
      typedef typename PINOCCHIO_ALIGNED_STD_VECTOR(Motion) MotionVector;
      typedef typename PINOCCHIO_ALIGNED_STD_VECTOR(Matrix6) Matrix6Vector;

      explicit CustomData(const Model & model,
                          const Data & data,
                          const Eigen::DenseIndex size)
      : liMi(size_t(model.njoints),SE3::Identity())
      , oMi(size_t(model.njoints),SE3::Identity())
      , a(size_t(model.njoints),Motion::Zero())
      , a_augmented(size_t(model.njoints),Motion::Zero())
      , Yaba(size_t(model.njoints),Inertia::Zero())
      , Yaba_augmented(size_t(model.njoints),Inertia::Zero())
      , joints(data.joints)
      , joints_augmented(data.joints)
      , u(model.nv)
      , ddq(model.nv)
      , f(size_t(model.njoints))
      , tmp_vec(size)
      {}
      
      SE3Vector liMi, oMi;
      MotionVector a, a_augmented;
      Matrix6Vector Yaba, Yaba_augmented;
      
      typename Data::JointDataVector joints;
      typename Data::JointDataVector joints_augmented;
      VectorXs u, ddq;
      ForceVector f;
      Vector tmp_vec;
    };
    
    const CustomData & getCustomData() const { return m_custom_data; }
    
  protected:
    
    static Eigen::DenseIndex evalConstraintSize(const ConstraintModelVector & constraint_models)
    {
      Eigen::DenseIndex size = 0;
      for(const auto & cm: constraint_models)
        size += cm.size();
      
      return size;
    }
    
    inline void compute_conclude()
    {
      m_dirty = false;
    }
    
    // Holders
    Eigen::DenseIndex m_size;
    ModelHolder m_model_ref;
    DataHolder m_data_ref;
    ConstraintModelVectorHolder m_constraint_models_ref;
    ConstraintDataVectorHolder m_constraint_datas_ref;
    
    mutable CustomData m_custom_data;
    bool m_dirty;
    Vector m_damping, m_damping_inverse;
  };
  
} // namespace pinocchio

#include "pinocchio/algorithm/delassus-operator-rigid-body.hxx"

#endif // ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__
