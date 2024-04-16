//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hxx__
#define __pinocchio_algorithm_delassus_operator_linear_complexity_hxx__

#include "pinocchio/algorithm/check.hpp"
#include "pinocchio/algorithm/contact-info.hpp"
#include "pinocchio/algorithm/model.hpp"
#include "pinocchio/multibody/fwd.hpp"

namespace pinocchio {
  
  template<typename DelassusOperator, typename ConfigVectorType>
  struct DelassusOperatorRigidBodyTplCalcForwardPass
  : public fusion::JointUnaryVisitorBase< DelassusOperatorRigidBodyTplCalcForwardPass<DelassusOperator,ConfigVectorType> >
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::Data Data;
    
    typedef boost::fusion::vector<const Model &,
    Data &,
    const ConfigVectorType &
    > ArgsType;
    
    template<typename JointModel>
    static void algo(const pinocchio::JointModelBase<JointModel> & jmodel,
                     pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const Model & model,
                     Data & data,
                     const Eigen::MatrixBase<ConfigVectorType> & q)
    {
      typedef typename Model::JointIndex JointIndex;
      
      const JointIndex i = jmodel.id();
      jmodel.calc(jdata.derived(),q.derived());
      
      const JointIndex parent = model.parents[i];
      data.liMi[i] = model.jointPlacements[i] * jdata.M();
      if(parent > 0)
        data.oMi[i] = data.oMi[parent] * data.liMi[i];
      else
        data.oMi[i] = data.liMi[i];
      
      jmodel.jointCols(data.J) = data.oMi[i].act(jdata.S());
      
      data.oYcrb[i] = data.oMi[i].act(model.inertias[i]);
      data.oYaba[i] = data.oYcrb[i].matrix();
    }
    
  };
  
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, template<typename T> class Holder>
  template<typename ConfigVectorType>
  void DelassusOperatorRigidBodyTpl<Scalar,Options,JointCollectionTpl,Holder>::calc(const Eigen::MatrixBase<ConfigVectorType> & q)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(q.size(), model().nq,
                                  "The joint configuration vector is not of right size");
    typedef DelassusOperatorRigidBodyTplCalcForwardPass<DelassusOperatorRigidBodyTpl,ConfigVectorType> Pass1;
    for(JointIndex i=1; i<(JointIndex)model().njoints; ++i)
    {
      typename Pass1::ArgsType args(this->model,this->data,q.derived());
      Pass1::run(this->model().joints[i],this->data().joints[i],args);
    }
  }
  
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hxx__
