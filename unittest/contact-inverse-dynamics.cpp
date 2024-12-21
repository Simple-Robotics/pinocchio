//
// Copyright (c) 2024 INRIA
//

#include <iostream>
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/centroidal.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/contact-info.hpp"
#include "pinocchio/algorithm/compute-all-terms.hpp"
#include "pinocchio/algorithm/constrained-dynamics.hpp"
#include "pinocchio/algorithm/contact-dynamics.hpp"
#include "pinocchio/algorithm/contact-inverse-dynamics.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/multibody/sample-models.hpp"
#include "pinocchio/utils/timer.hpp"
#include "pinocchio/spatial/classic-acceleration.hpp"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(test_contact_inverse_dynamics_3D)
{
  using namespace Eigen;
  using namespace pinocchio;

  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);
  pinocchio::Data data(model), data_ref(model);

  model.lowerPositionLimit.head<3>().fill(-1.);
  model.upperPositionLimit.head<3>().fill(1.);
  VectorXd q = randomConfiguration(model);

  VectorXd v = VectorXd::Random(model.nv);
  VectorXd tau = VectorXd::Random(model.nv);

  const std::string RF = "rleg6_joint";
  //  const Model::JointIndex RF_id = model.getJointId(RF);
  const std::string LF = "lleg6_joint";
  //  const Model::JointIndex LF_id = model.getJointId(LF);

  // Contact models and data
  typedef PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(FrictionalPointConstraintModel)
    FrictionalConstraintModelVector;
  typedef PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(FrictionalPointConstraintData)
    FrictionalConstraintDataVector;

  FrictionalConstraintModelVector contact_models;
  FrictionalConstraintDataVector contact_datas;
  FrictionalPointConstraintModel ci_RF(model, model.getJointId(RF));
  ci_RF.set().mu = 0.4;
  contact_models.push_back(ci_RF);
  contact_datas.push_back(FrictionalPointConstraintData(ci_RF));
  FrictionalPointConstraintModel ci_LF(model, model.getJointId(LF));
  ci_LF.set().mu = 0.4;
  contact_models.push_back(ci_LF);
  contact_datas.push_back(FrictionalPointConstraintData(ci_LF));

  FrictionalConstraintDataVector contact_datas_ref(contact_datas);

  Eigen::DenseIndex constraint_dim = 0;
  for (size_t k = 0; k < contact_models.size(); ++k)
    constraint_dim += contact_models[k].size();

  Eigen::MatrixXd J_ref(constraint_dim, model.nv);
  J_ref.setZero();

  double dt = 1e-3;
  Eigen::VectorXd R = Eigen::VectorXd::Zero(constraint_dim);
  Eigen::VectorXd constraint_correction = Eigen::VectorXd::Zero(constraint_dim);
  boost::optional<Eigen::VectorXd> lambda_guess = boost::none;
  Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  ProximalSettings prox_settings(1e-12, 1e-6, 1);
  initConstraintDynamics(model, data_ref, contact_models);
  contactInverseDynamics(
    model, data_ref, q, v, a, dt, contact_models, contact_datas, R, constraint_correction,
    prox_settings, lambda_guess);
}

BOOST_AUTO_TEST_SUITE_END()
