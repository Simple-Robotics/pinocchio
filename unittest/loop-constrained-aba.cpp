//
// Copyright (c) 2024-2025 INRIA
//

#include <iostream>

#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/contact-info.hpp"
#include "pinocchio/algorithm/proximal.hpp"
#include "pinocchio/algorithm/constrained-dynamics.hpp"
#include "pinocchio/algorithm/contact-dynamics.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/multibody/sample-models.hpp"
#include "pinocchio/spatial/classic-acceleration.hpp"
#include "pinocchio/spatial/explog.hpp"
#include "pinocchio/algorithm/loop-constrained-aba.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

using namespace pinocchio;
using namespace Eigen;

void build_trident_model(Model & model)
{

  Inertia I1 = Inertia::Random();
  Inertia I2 = Inertia::Random();
  model.gravity.linear() = model.gravity981;

  JointModel branch_connector = JointModelRX();

  JointIndex joint0 = model.addJoint(0, JointModelFreeFlyer(), SE3::Random(), "joint0");
  model.appendBodyToJoint(joint0, I1, SE3::Random());

  int num_children = 10;

  JointIndex parent_id;
  std::string joint_name;
  JointIndex mid_chain_joint;

  JointIndex joint1 = model.addJoint(joint0, branch_connector, SE3::Random(), "joint1");
  parent_id = joint1;
  model.appendBodyToJoint(joint1, I1, SE3::Random());
  for (int k = 1; k <= num_children; ++k)
  {
    joint_name = "joint1" + std::to_string(k);
    parent_id = model.addJoint(parent_id, JointModelRX(), SE3::Random(), joint_name);
    model.appendBodyToJoint(parent_id, Inertia::Random(), SE3::Random());

    if (k == 3)
      mid_chain_joint = parent_id;
  }

  JointIndex joint2 = model.addJoint(mid_chain_joint, branch_connector, SE3::Random(), "joint2");
  parent_id = joint2;
  model.appendBodyToJoint(joint2, I2, SE3::Random());
  for (int k = 1; k <= num_children; ++k)
  {
    joint_name = "joint2" + std::to_string(k);
    parent_id = model.addJoint(parent_id, JointModelRX(), SE3::Random(), joint_name);
    model.appendBodyToJoint(parent_id, Inertia::Random(), SE3::Random());

    if (k == 3)
      mid_chain_joint = parent_id;
  }

  JointIndex joint3 = model.addJoint(mid_chain_joint, branch_connector, SE3::Random(), "joint3");
  parent_id = joint3;
  model.appendBodyToJoint(joint3, Inertia::Random(), SE3::Random());
  for (int k = 1; k <= num_children; ++k)
  {
    joint_name = "joint3" + std::to_string(k);
    parent_id = model.addJoint(parent_id, JointModelRX(), SE3::Random(), joint_name);
    model.appendBodyToJoint(parent_id, Inertia::Random(), SE3::Random());
  }

  model.lowerPositionLimit.fill(-5.);
  model.upperPositionLimit.fill(5.);
}

BOOST_AUTO_TEST_CASE(test_6D_unconstrained)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  const double mu0 = 1e-3;
  ProximalSettings prox_settings(1e-12, mu0, 1);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);
  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_6D_descendants)
{
  Model model;

  build_trident_model(model);
  // model.gravity.setZero();
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint12"), model.getJointId("joint17"),
    LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();

  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  const double mu0 = 1e-5;
  ProximalSettings prox_settings(1e-14, mu0, 100);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);

  std::cout << "data_ref.ddq: " << data_ref.ddq.transpose() << std::endl;
  std::cout << "data.ddq: " << data.ddq.transpose() << std::endl;
  std::cout << "|| data_ref.ddq - data.ddq ||: " << (data_ref.ddq - data.ddq).norm() << std::endl;
  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-8));
}

BOOST_AUTO_TEST_CASE(test_6D_descendants_reversed)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint17"), model.getJointId("joint12"),
    LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();
  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  const double mu0 = 1e-1;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);

  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_12D_descendants_redundant_reversed)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint12"), model.getJointId("joint17"),
    LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();
  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  RigidConstraintModel rcm2 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint17"), model.getJointId("joint12"),
    LOCAL_WORLD_ALIGNED);
  rcm2.joint1_placement.setRandom();
  rcm2.joint2_placement.setRandom();
  contact_models.push_back(rcm2);
  contact_datas.push_back(RigidConstraintData(rcm2));

  const double mu0 = 1e-3;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);
  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_6D_different_branches)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint210"), model.getJointId("joint38"),
    LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();

  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  const double mu0 = 1e-3;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);

  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_12D_coupled_loop_common_link)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = pinocchio::randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint210"), model.getJointId("joint38"),
    LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();
  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  RigidConstraintModel rcm2 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint38"), model.getJointId("joint18"),
    LOCAL_WORLD_ALIGNED);
  rcm2.joint1_placement.setRandom();
  rcm2.joint2_placement.setRandom();
  contact_models.push_back(rcm2);
  contact_datas.push_back(RigidConstraintData(rcm2));

  const double mu0 = 1e-1;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);

  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_24D_coupling_with_double_ground)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = pinocchio::randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint1"), model.getJointId("joint14"),
    LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();
  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  RigidConstraintModel rcm2 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint1"), model.getJointId("joint24"),
    LOCAL_WORLD_ALIGNED);
  rcm2.joint1_placement.setRandom();
  rcm2.joint2_placement.setRandom();
  contact_models.push_back(rcm2);
  contact_datas.push_back(RigidConstraintData(rcm2));

  RigidConstraintModel rcm3 =
    RigidConstraintModel(CONTACT_6D, model, model.getJointId("joint24"), LOCAL_WORLD_ALIGNED);
  rcm3.joint1_placement.setRandom();
  rcm3.joint2_placement.setRandom();
  contact_models.push_back(rcm3);
  contact_datas.push_back(RigidConstraintData(rcm3));

  RigidConstraintModel rcm4 =
    RigidConstraintModel(CONTACT_6D, model, model.getJointId("joint14"), LOCAL_WORLD_ALIGNED);
  rcm4.joint1_placement.setRandom();
  rcm4.joint2_placement.setRandom();
  contact_models.push_back(rcm4);
  contact_datas.push_back(RigidConstraintData(rcm4));

  const double mu0 = 1e-1;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);

  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_6D_consecutive_links)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint14"), model.getJointId("joint15"),
    LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();
  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  const double mu0 = 1e-3;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);

  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_12D_coupled_on_a_chain)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint12"), model.getJointId("joint19"),
    LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();
  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  RigidConstraintModel rcm2 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint14"), model.getJointId("joint18"),
    LOCAL_WORLD_ALIGNED);
  rcm2.joint1_placement.setRandom();
  rcm2.joint2_placement.setRandom();
  contact_models.push_back(rcm2);
  contact_datas.push_back(RigidConstraintData(rcm2));

  const double mu0 = 1e-3;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);

  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_12D_cross_coupled_on_a_chain)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint12"), model.getJointId("joint19"),
    LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();
  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  RigidConstraintModel rcm2 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint14"), model.getJointId("joint18"),
    LOCAL_WORLD_ALIGNED);
  rcm2.joint1_placement.setRandom();
  rcm2.joint2_placement.setRandom();
  contact_models.push_back(rcm2);
  contact_datas.push_back(RigidConstraintData(rcm2));

  const double mu0 = 1e-3;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);

  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_24D_cross_coupling)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint11"), model.getJointId("joint39"),
    LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();
  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  RigidConstraintModel rcm2 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint21"), model.getJointId("joint38"),
    LOCAL_WORLD_ALIGNED);
  rcm2.joint1_placement.setRandom();
  rcm2.joint2_placement.setRandom();
  contact_models.push_back(rcm2);
  contact_datas.push_back(RigidConstraintData(rcm2));

  RigidConstraintModel rcm3 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint19"), model.getJointId("joint29"),
    LOCAL_WORLD_ALIGNED);
  rcm3.joint1_placement.setRandom();
  rcm3.joint2_placement.setRandom();
  contact_models.push_back(rcm3);
  contact_datas.push_back(RigidConstraintData(rcm3));

  RigidConstraintModel rcm4 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint29"), model.getJointId("joint39"),
    LOCAL_WORLD_ALIGNED);
  rcm4.joint1_placement.setRandom();
  rcm4.joint2_placement.setRandom();
  contact_models.push_back(rcm4);
  contact_datas.push_back(RigidConstraintData(rcm4));

  const double mu0 = 1e-3;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);

  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_6D_cons_baumgarte)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 =
    RigidConstraintModel(CONTACT_6D, model, model.getJointId("joint11"), LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();
  rcm1.corrector.Kd.setIdentity();
  rcm1.corrector.Kp.setIdentity();
  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  RigidConstraintModel rcm2 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint21"), model.getJointId("joint38"),
    LOCAL_WORLD_ALIGNED);
  rcm2.joint1_placement.setRandom();
  rcm2.joint2_placement.setRandom();
  rcm2.corrector.Kd.setIdentity();
  rcm2.corrector.Kp.setIdentity();
  contact_models.push_back(rcm2);
  contact_datas.push_back(RigidConstraintData(rcm2));

  RigidConstraintModel rcm3 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint19"), model.getJointId("joint29"),
    LOCAL_WORLD_ALIGNED);
  rcm3.joint1_placement.setRandom();
  rcm3.joint2_placement.setRandom();
  rcm3.corrector.Kd.setIdentity();
  rcm3.corrector.Kp.setIdentity();
  contact_models.push_back(rcm3);
  contact_datas.push_back(RigidConstraintData(rcm3));

  RigidConstraintModel rcm4 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint29"), model.getJointId("joint39"),
    LOCAL_WORLD_ALIGNED);
  rcm4.joint1_placement.setRandom();
  rcm4.joint2_placement.setRandom();
  rcm4.corrector.Kd.setIdentity();
  rcm4.corrector.Kp.setIdentity();
  contact_models.push_back(rcm4);
  contact_datas.push_back(RigidConstraintData(rcm4));

  const double mu0 = 1e-3;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);

  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_3D_cons_baumgarte)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 =
    RigidConstraintModel(CONTACT_3D, model, model.getJointId("joint11"), LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();
  rcm1.corrector.Kd.setIdentity();
  rcm1.corrector.Kp.setIdentity();
  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  RigidConstraintModel rcm2 = RigidConstraintModel(
    CONTACT_3D, model, model.getJointId("joint21"), model.getJointId("joint38"),
    LOCAL_WORLD_ALIGNED);
  rcm2.joint1_placement.setRandom();
  rcm2.joint2_placement.setRandom();
  rcm2.corrector.Kd.setIdentity();
  rcm2.corrector.Kp.setIdentity();
  contact_models.push_back(rcm2);
  contact_datas.push_back(RigidConstraintData(rcm2));

  RigidConstraintModel rcm3 = RigidConstraintModel(
    CONTACT_3D, model, model.getJointId("joint19"), model.getJointId("joint29"),
    LOCAL_WORLD_ALIGNED);
  rcm3.joint1_placement.setRandom();
  rcm3.joint2_placement.setRandom();
  rcm3.corrector.Kd.setIdentity();
  rcm3.corrector.Kp.setIdentity();
  contact_models.push_back(rcm3);
  contact_datas.push_back(RigidConstraintData(rcm3));

  RigidConstraintModel rcm4 = RigidConstraintModel(
    CONTACT_3D, model, model.getJointId("joint29"), model.getJointId("joint39"),
    LOCAL_WORLD_ALIGNED);
  rcm4.joint1_placement.setRandom();
  rcm4.joint2_placement.setRandom();
  rcm4.corrector.Kd.setIdentity();
  rcm4.corrector.Kp.setIdentity();
  contact_models.push_back(rcm4);
  contact_datas.push_back(RigidConstraintData(rcm4));

  const double mu0 = 1e-3;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);
  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_loop_con_and_ground_con)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint12"), model.getJointId("joint17"),
    LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();

  RigidConstraintModel rcm2 =
    RigidConstraintModel(CONTACT_6D, model, model.getJointId("joint14"), LOCAL_WORLD_ALIGNED);
  rcm2.joint1_placement.setRandom();

  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  contact_models.push_back(rcm2);
  contact_datas.push_back(RigidConstraintData(rcm2));

  const double mu0 = 1e-1;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);

  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_loop_con_and_ground_con3D)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint12"), model.getJointId("joint17"),
    LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();

  RigidConstraintModel rcm2 =
    RigidConstraintModel(CONTACT_3D, model, model.getJointId("joint24"), LOCAL_WORLD_ALIGNED);
  rcm2.joint1_placement.setRandom();

  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  contact_models.push_back(rcm2);
  contact_datas.push_back(RigidConstraintData(rcm2));

  const double mu0 = 1e-1;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);

  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_loop_con3D_ground_con3D)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 = RigidConstraintModel(
    CONTACT_3D, model, model.getJointId("joint12"), model.getJointId("joint17"),
    LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();

  RigidConstraintModel rcm2 =
    RigidConstraintModel(CONTACT_3D, model, model.getJointId("joint24"), LOCAL_WORLD_ALIGNED);
  rcm2.joint1_placement.setRandom();

  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  contact_models.push_back(rcm2);
  contact_datas.push_back(RigidConstraintData(rcm2));

  const double mu0 = 1e-1;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);

  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_CASE(test_coupled_3D_6D_loops)
{
  Model model;

  build_trident_model(model);
  Data data(model), data_ref(model);

  // Contact models and data
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_models;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_datas;

  const VectorXd q = randomConfiguration(model);
  const VectorXd v = VectorXd::Random(model.nv);
  const VectorXd tau = VectorXd::Random(model.nv);

  RigidConstraintModel rcm1 = RigidConstraintModel(
    CONTACT_3D, model, model.getJointId("joint12"), model.getJointId("joint27"),
    LOCAL_WORLD_ALIGNED);
  rcm1.joint1_placement.setRandom();
  rcm1.joint2_placement.setRandom();

  RigidConstraintModel rcm2 = RigidConstraintModel(
    CONTACT_6D, model, model.getJointId("joint24"), model.getJointId("joint11"),
    LOCAL_WORLD_ALIGNED);
  rcm2.joint1_placement.setRandom();
  rcm2.joint2_placement.setRandom();

  contact_models.push_back(rcm1);
  contact_datas.push_back(RigidConstraintData(rcm1));

  contact_models.push_back(rcm2);
  contact_datas.push_back(RigidConstraintData(rcm2));

  const double mu0 = 1e-1;
  ProximalSettings prox_settings(1e-12, mu0, 3);

  initConstraintDynamics(model, data_ref, contact_models);
  constraintDynamics(model, data_ref, q, v, tau, contact_models, contact_datas, prox_settings);

  initClCaba(model, data, contact_models);
  closedLoopCaba(model, data, q, v, tau, contact_models, contact_datas, prox_settings);

  BOOST_CHECK(data_ref.ddq.isApprox(data.ddq, 1e-10));
}

BOOST_AUTO_TEST_SUITE_END()
