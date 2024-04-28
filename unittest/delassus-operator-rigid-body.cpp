//
// Copyright (c) 2024 INRIA
//

#include <iostream>
#include <memory>

#include "pinocchio/algorithm/cholesky.hpp"
#include "pinocchio/algorithm/contact-info.hpp"
#include "pinocchio/algorithm/contact-dynamics.hpp"
#include "pinocchio/algorithm/contact-jacobian.hpp"
#include "pinocchio/algorithm/contact-cholesky.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/context.hpp"
#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/algorithm/delassus-operator-rigid-body.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;
using namespace Eigen;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(default_constructor_shared_ptr)
{
  typedef DelassusOperatorRigidBodyTpl<double,0,JointCollectionDefaultTpl,std::shared_ptr> DelassusOperatorRigidBodySharedPtr;
  typedef typename DelassusOperatorRigidBodySharedPtr::ConstraintModelVector ConstraintModelVector;
  typedef typename DelassusOperatorRigidBodySharedPtr::ConstraintDataVector ConstraintDataVector;

  std::shared_ptr<Model> model_shared_ptr = std::make_shared<Model>();
  Model * model_ptr = model_shared_ptr.get();
  Model & model = *model_ptr;
  buildModels::humanoidRandom(model,true);
  
  std::shared_ptr<Data> data_shared_ptr = std::make_shared<Data>(model);
  Data * data_ptr = data_shared_ptr.get();
//  Data & data = *data_ptr;
  
  std::shared_ptr<ConstraintModelVector> constraint_models_shared_ptr = std::make_shared<ConstraintModelVector>();
  ConstraintModelVector * constraint_models_ptr = constraint_models_shared_ptr.get();
  std::shared_ptr<ConstraintDataVector> constraint_datas_shared_ptr = std::make_shared<ConstraintDataVector>();
  ConstraintDataVector * constraint_datas_ptr = constraint_datas_shared_ptr.get();

  DelassusOperatorRigidBodySharedPtr delassus_operator(model_shared_ptr,
                                                       data_shared_ptr,
                                                       constraint_models_shared_ptr,
                                                       constraint_datas_shared_ptr);
  
  BOOST_CHECK(delassus_operator.size() == 0);
  BOOST_CHECK(&delassus_operator.model() == model_ptr);
  BOOST_CHECK(&delassus_operator.data() == data_ptr);
  BOOST_CHECK(&delassus_operator.constraint_models() == constraint_models_ptr);
  BOOST_CHECK(&delassus_operator.constraint_datas() == constraint_datas_ptr);
}

BOOST_AUTO_TEST_CASE(default_constructor_reference_wrapper)
{
  typedef DelassusOperatorRigidBodyTpl<double,0,JointCollectionDefaultTpl,std::reference_wrapper> DelassusOperatorRigidBodyReferenceWrapper;
  typedef typename DelassusOperatorRigidBodyReferenceWrapper::ConstraintModelVector ConstraintModelVector;
  typedef typename DelassusOperatorRigidBodyReferenceWrapper::ConstraintDataVector ConstraintDataVector;
  
  Model model;
  std::reference_wrapper<Model> model_ref = model;
  buildModels::humanoidRandom(model,true);
  
  Data data(model);
  std::reference_wrapper<Data> data_ref = data;
  
  ConstraintModelVector constraint_models;
  std::reference_wrapper<ConstraintModelVector> constraint_models_ref = constraint_models;
  ConstraintDataVector constraint_datas;
  std::reference_wrapper<ConstraintDataVector> constraint_datas_ref = constraint_datas;

  DelassusOperatorRigidBodyReferenceWrapper delassus_operator(model_ref, data_ref,
                                                              constraint_models_ref,
                                                              constraint_datas_ref);
  
  BOOST_CHECK(delassus_operator.size() == 0);
  BOOST_CHECK(&delassus_operator.model() == &model);
  BOOST_CHECK(&delassus_operator.data() == &data);
  BOOST_CHECK(&delassus_operator.constraint_models() == &constraint_models);
  BOOST_CHECK(&delassus_operator.constraint_datas() == &constraint_datas);
}

BOOST_AUTO_TEST_CASE(test_compute)
{
  typedef DelassusOperatorRigidBodyTpl<double,0,JointCollectionDefaultTpl,std::reference_wrapper> DelassusOperatorRigidBodyReferenceWrapper;
  typedef DelassusOperatorRigidBodyReferenceWrapper::CustomData CustomData;
  typedef typename DelassusOperatorRigidBodyReferenceWrapper::ConstraintModelVector ConstraintModelVector;
  typedef typename DelassusOperatorRigidBodyReferenceWrapper::ConstraintDataVector ConstraintDataVector;
  
  Model model;
  std::reference_wrapper<Model> model_ref = model;
  buildModels::humanoidRandom(model,true);
  model.gravity.setZero();
  const Eigen::VectorXd q_neutral = neutral(model);
  
  const std::string RF = "rleg6_joint";
  const std::string LF = "lleg6_joint";
  
  Data data(model), data_gt(model), data_aba(model);
  std::reference_wrapper<Data> data_ref = data;
  
  ConstraintModelVector constraint_models;
  ConstraintDataVector constraint_datas;
  const RigidConstraintModel cm_RF_LOCAL(CONTACT_3D,model,model.getJointId(RF),SE3::Random(),LOCAL);
  constraint_models.push_back(cm_RF_LOCAL);
  constraint_datas.push_back(cm_RF_LOCAL.createData());
  const RigidConstraintModel cm_LF_LOCAL(CONTACT_3D,model,model.getJointId(LF),SE3::Random(),LOCAL);
  constraint_models.push_back(cm_LF_LOCAL);
  constraint_datas.push_back(cm_LF_LOCAL.createData());
  
  ConstraintDataVector constraint_datas_gt = constraint_datas;
  
  std::reference_wrapper<ConstraintModelVector> constraint_models_ref = constraint_models;
  std::reference_wrapper<ConstraintDataVector> constraint_datas_ref = constraint_datas;
  
  DelassusOperatorRigidBodyReferenceWrapper delassus_operator(model_ref, data_ref,
                                                              constraint_models_ref,
                                                              constraint_datas_ref);
  // Eval J Minv Jt
  auto Minv_gt = computeMinverse(model, data_gt, q_neutral);
  make_symmetric(Minv_gt);
  BOOST_CHECK(Minv_gt.isApprox(Minv_gt.transpose()));
  
  Eigen::MatrixXd constraints_jacobian_gt(delassus_operator.size(),model.nv);
  constraints_jacobian_gt.setZero();
  evalConstraints(model,data_gt,constraint_models,constraint_datas_gt);
  getConstraintsJacobian(model,data_gt,constraint_models,constraint_datas_gt,constraints_jacobian_gt);
  
  const Eigen::MatrixXd delassus_dense_gt = constraints_jacobian_gt * Minv_gt * constraints_jacobian_gt.transpose();
  
  // Test Jacobian transpose operator
  {
    const Eigen::VectorXd rhs = Eigen::VectorXd::Random(delassus_operator.size());
    
    computeJointJacobians(model,data,q_neutral);
    for(Model::JointIndex joint_id = 1; joint_id < Model::JointIndex(model.njoints); ++joint_id)
    {
      BOOST_CHECK(data.joints[joint_id].S().isApprox(data_gt.joints[joint_id].S()));
      BOOST_CHECK(data.liMi[joint_id].isApprox(data_gt.liMi[joint_id]));
      BOOST_CHECK(data.oMi[joint_id].isApprox(data_gt.oMi[joint_id]));
    }
    
    const Eigen::VectorXd Jt_rhs_gt = constraints_jacobian_gt.transpose() * rhs;
    Eigen::VectorXd Jt_rhs(model.nv);
    
    evalConstraints(model,data,constraint_models,constraint_datas);
    evalConstraintJacobianTransposeProduct(model,data,constraint_models,constraint_datas,rhs,Jt_rhs);
    
    BOOST_CHECK(Jt_rhs.isApprox(Jt_rhs_gt));
    
  };
  
  // Test operator *
  {
    const Eigen::VectorXd rhs = Eigen::VectorXd::Random(delassus_operator.size());
    Eigen::VectorXd res(delassus_operator.size());
    
    delassus_operator.compute(q_neutral);
    delassus_operator.applyOnTheRight(rhs,res);
    
    // Eval Jt*rhs vs internal computations. This test is useful to check intermediate computation.
//    Eigen::VectorXd Jt_rhs(model.nv);
//    evalConstraintJacobianTransposeProduct(model,data,constraint_models,constraint_datas,rhs,Jt_rhs);
//    BOOST_CHECK(delassus_operator.getCustomData().u.isApprox(Jt_rhs));
//    
//    std::cout << "delassus_operator.getCustomData().u: " << delassus_operator.getCustomData().u.transpose() << std::endl;
//    std::cout << "Jt_rhs: " << Jt_rhs.transpose() << std::endl;
//    const Eigen::VectorXd Jt_rhs_gt = constraints_jacobian_gt.transpose() * rhs;
//    BOOST_CHECK(delassus_operator.getCustomData().u.isApprox(Jt_rhs_gt));

    pinocchio::container::aligned_vector<Data::Force> joint_forces_gt(size_t(model.njoints),Data::Force::Zero());
    mapConstraintForceToJointForces(model,data_gt,constraint_models,constraint_datas_gt,rhs,joint_forces_gt);
    minimal::aba(model, data_aba, q_neutral, Eigen::VectorXd::Zero(model.nv), Eigen::VectorXd::Zero(model.nv), joint_forces_gt);
    
    for(Model::JointIndex joint_id = 1; joint_id < Model::JointIndex(model.njoints); ++joint_id)
    {
      const CustomData & custom_data = delassus_operator.getCustomData();
      BOOST_CHECK(custom_data.joints[joint_id].S().isApprox(data_aba.joints[joint_id].S()));
      BOOST_CHECK(custom_data.liMi[joint_id].isApprox(data_aba.liMi[joint_id]));
//      BOOST_CHECK(custom_data.oMi[joint_id].isApprox(data_aba.oMi[joint_id])); // minimal::ABA does not compute this quantity
      BOOST_CHECK(custom_data.Yaba[joint_id].isApprox(data_aba.Yaba[joint_id]));
    }
    BOOST_CHECK(delassus_operator.getCustomData().u.isApprox(data_aba.u));

//    std::cout << "delassus_operator.getCustomData().u: " << delassus_operator.getCustomData().u.transpose() << std::endl;
//    std::cout << "data_aba.u: " << data_aba.u.transpose() << std::endl;
    
    const Eigen::VectorXd Jt_rhs_gt = constraints_jacobian_gt.transpose() * rhs;
    const Eigen::VectorXd Minv_Jt_rhs_gt = Minv_gt * Jt_rhs_gt;
    
//    std::cout << "Minv_Jt_rhs: " << delassus_operator.getCustomData().ddq.transpose() << std::endl;
//    std::cout << "Minv_Jt_rhs_gt: " << Minv_Jt_rhs_gt.transpose() << std::endl;
    
    BOOST_CHECK(delassus_operator.getCustomData().ddq.isApprox(Minv_Jt_rhs_gt));
    
    const auto res_gt = (delassus_dense_gt * rhs).eval();
    BOOST_CHECK(res.isApprox(res_gt));
    
//    std::cout << "res: " << res.transpose() << std::endl;
//    std::cout << "res_gt: " << res_gt.transpose() << std::endl;
    
    // Multiple call and operator *
    {
      for(int i = 0; i < 100; ++i)
      {
        Eigen::VectorXd res(delassus_operator.size());
        delassus_operator.applyOnTheRight(rhs,res);
        BOOST_CHECK(res.isApprox(res_gt));
        
        const Eigen::VectorXd res2 = delassus_operator * rhs;
        BOOST_CHECK(res2 == res); // Should be exactly the same
        BOOST_CHECK(res2.isApprox(res_gt));
      }
    }
  }
  
  // Update damping
  {
    const Eigen::VectorXd rhs = Eigen::VectorXd::Random(delassus_operator.size());
    const double mu = 1 ;
    delassus_operator.updateDamping(mu);
    BOOST_CHECK(delassus_operator.getDamping().isApproxToConstant(mu));
    
    Eigen::VectorXd res_damped(delassus_operator.size());
    delassus_operator.applyOnTheRight(rhs,res_damped);
    const auto res_gt_damped = ((delassus_dense_gt + mu * Eigen::MatrixXd::Identity(delassus_operator.size(),delassus_operator.size())) * rhs).eval();
    BOOST_CHECK(res_damped.isApprox(res_gt_damped));
  }
}

BOOST_AUTO_TEST_SUITE_END()
