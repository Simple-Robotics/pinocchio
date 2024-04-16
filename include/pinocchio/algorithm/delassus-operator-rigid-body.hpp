//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__
#define __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__

#include <functional>
#include <memory>

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
    typedef DelassusOperatorRigidBodyTpl Self;
    typedef DelassusOperatorBase<Self> Base;
    
    typedef typename traits<Self>::Scalar Scalar;
    enum { Options = traits<Self>::Options };
    
    typedef typename traits<Self>::Model Model;
    typedef Holder<const Model> ModelHolder;
    typedef typename traits<Self>::Data Data;
    typedef Holder<Data> DataHolder;
    
    typedef typename traits<Self>::ConstraintModel ConstraintModel;
    typedef typename traits<Self>::ConstraintModelVector ConstraintModelVector;
    typedef Holder<const ConstraintModelVector> ConstraintModelVectorHolder;
    
    typedef typename traits<Self>::ConstraintData ConstraintData;
    typedef typename traits<Self>::ConstraintDataVector ConstraintDataVector;
    typedef Holder<ConstraintDataVector> ConstraintDataVectorHolder;

    DelassusOperatorRigidBodyTpl(const ModelHolder &model_ref,
                                        const DataHolder &data_ref,
                                        const ConstraintModelVectorHolder &constraint_models_ref,
                                        const ConstraintDataVectorHolder &constraint_datas_ref)
    : Base()
    , m_size(evalConstraintSize(get_ref(constraint_models_ref)))
    , m_model_ref(model_ref)
    , m_data_ref(data_ref)
    , m_constraint_models_ref(constraint_models_ref)
    , m_constraint_datas_ref(constraint_datas_ref)
    {
      assert(model().check(data()) && "data is not consistent with model.");
      PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_models().size(), constraint_datas().size(),
                                    "The sizes of contact vector models and contact vector data are not the same.");
    }
    
    /// \brief Update the intermediate computations according to a new configuration vector entry
    ///
    /// \param[in] q Configuration vector
    ///
    template<typename ConfigVectorType>
    void calc(const Eigen::MatrixBase<ConfigVectorType> & q);
    
    const Model & model() const { return get_ref(m_model_ref); }
    
    Data & data() { return get_ref(m_data_ref); }
    const Data & data() const { return get_ref(m_data_ref); }
    
    const ConstraintModelVector & constraint_models() const { return get_ref(m_constraint_models_ref); }
    
    const ConstraintDataVector & constraint_datas() const { return get_ref(m_constraint_datas_ref); }
    ConstraintDataVector & constraint_datas() { return get_ref(m_constraint_datas_ref); }
    
    Eigen::DenseIndex size() const { return m_size; }

  protected:
    
    static Eigen::DenseIndex evalConstraintSize(const ConstraintModelVector & constraint_models)
    {
      Eigen::DenseIndex size = 0;
      for(const auto & cm: constraint_models)
        size += cm.size();
      
      return size;
    }
    
    // Holders
    Eigen::DenseIndex m_size;
    ModelHolder m_model_ref;
    DataHolder m_data_ref;
    ConstraintModelVectorHolder m_constraint_models_ref;
    ConstraintDataVectorHolder m_constraint_datas_ref;
  };
  
} // namespace pinocchio

#include "pinocchio/algorithm/delassus-operator-rigid-body.hxx"

#endif // ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__
