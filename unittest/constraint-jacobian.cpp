//
// Copyright (c) 2024 INRIA
//

#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/contact-info.hpp"
#include "pinocchio/algorithm/contact-jacobian.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/parsers/sample-models.hpp"

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;
using namespace Eigen;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(constraint_jacobian_operations)
{
  
  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model,true);
  Data data(model), data_ref(model);
  
  model.lowerPositionLimit.head<3>().fill(-1.);
  model.upperPositionLimit.head<3>().fill( 1.);
  VectorXd q = randomConfiguration(model);
  computeJointJacobians(model,data,q);
  computeJointJacobians(model,data_ref,q);

  const std::string RF = "rleg6_joint";
  const std::string LF = "lleg6_joint";
  
  // 3D - LOCAL
  {
    RigidConstraintModel cm_RF_LOCAL(CONTACT_3D,model,model.getJointId(RF),SE3::Random(),LOCAL);
    RigidConstraintData cd_RF_LOCAL(cm_RF_LOCAL);
    RigidConstraintModel cm_LF_LOCAL(CONTACT_3D,model,model.getJointId(LF),SE3::Random(),LOCAL);
    RigidConstraintData cd_LF_LOCAL(cm_LF_LOCAL);
    
    const std::vector<RigidConstraintModel> constraints_models{cm_RF_LOCAL,cm_LF_LOCAL};
    std::vector<RigidConstraintData> constraints_datas{cd_RF_LOCAL,cd_LF_LOCAL};
    std::vector<RigidConstraintData> constraints_datas_ref{cd_RF_LOCAL,cd_LF_LOCAL};
    
    const Eigen::DenseIndex m = Eigen::DenseIndex(getTotalConstraintSize(constraints_models));
    
    Eigen::VectorXd res(model.nv);
    const Eigen::VectorXd rhs = Eigen::VectorXd::Random(m);
    
    evalConstraints(model,data,constraints_models,constraints_datas);
    evalConstraintJacobianTransposeProduct(model,data,constraints_models,constraints_datas,rhs,res);
    
    // Check Jacobian
    {
      Eigen::VectorXd res_ref = Eigen::VectorXd::Zero(model.nv);
      Data::MatrixXs J_RF_LOCAL_sparse(3,model.nv); J_RF_LOCAL_sparse.setZero(); // TODO: change input type when all the API would be refactorized with CRTP on contact constraints
      getConstraintJacobian(model,data,cm_RF_LOCAL,cd_RF_LOCAL,J_RF_LOCAL_sparse);
      res_ref += J_RF_LOCAL_sparse.transpose() * rhs.segment<3>(0);
      
      Data::MatrixXs J_LF_LOCAL_sparse(3,model.nv); J_LF_LOCAL_sparse.setZero(); // TODO: change input type when all the API would be refactorized with CRTP on contact constraints
      getConstraintJacobian(model,data,cm_LF_LOCAL,cd_LF_LOCAL,J_LF_LOCAL_sparse);
      res_ref += J_LF_LOCAL_sparse.transpose() * rhs.segment<3>(3);
      
      BOOST_CHECK(res.isApprox(res_ref));
    }
    
    // Alternative way to compute the Jacobians
    {
      Eigen::MatrixXd J_ref(6,model.nv); J_ref.setZero();
      getConstraintsJacobian(model, data_ref, constraints_models, constraints_datas_ref, J_ref);
      const Eigen::VectorXd res_ref = J_ref.transpose() * rhs;
      BOOST_CHECK(res.isApprox(res_ref));
    }
  }
  
}

BOOST_AUTO_TEST_SUITE_END()
