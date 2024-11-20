//
// Copyright (c) 2024 INRIA
//

#include "pinocchio/algorithm/constraints/constraints.hpp"
#include "pinocchio/algorithm/contact-cholesky.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/contact-jacobian.hpp"
#include "pinocchio/algorithm/admm-solver.hpp"
#include "pinocchio/algorithm/pgs-solver.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/delassus.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;

double mu = 1e-4;

template<typename _ConstraintModel>
struct TestBoxTpl
{
  typedef _ConstraintModel ConstraintModel;

  typedef typename ConstraintModel::ConstraintData ConstraintData;
  typedef typename ConstraintModel::ConstraintSet ConstraintSet;

  TestBoxTpl(const Model & model, const std::vector<ConstraintModel> & constraint_models)
  : model(model)
  , data(model)
  , constraint_models(constraint_models)
  , v_next(Eigen::VectorXd::Zero(model.nv))
  {
    for (const auto & cm : constraint_models)
    {
      constraint_datas.push_back(cm.createData());
      constraint_sets.push_back(cm.set());
    }

    const Eigen::DenseIndex constraint_size = getTotalConstraintSize(constraint_models);
    primal_solution = dual_solution = dual_solution_sparse = Eigen::VectorXd::Zero(constraint_size);
  }

  void operator()(
    const Eigen::VectorXd & q0,
    const Eigen::VectorXd & v0,
    const Eigen::VectorXd & tau0,
    const Force & fext,
    const double dt)
  {
    std::vector<Force> external_forces(size_t(model.njoints), Force::Zero());
    external_forces[1] = fext;

    const Eigen::VectorXd v_free =
      dt * aba(model, data, q0, v0, tau0, external_forces, Convention::WORLD);

    // Cholesky of the Delassus matrix
    crba(model, data, q0, Convention::WORLD);
    ContactCholeskyDecomposition chol(model, constraint_models);
    chol.compute(model, data, constraint_models, constraint_datas, 1e-10);

    const Eigen::MatrixXd delassus_matrix_plain = chol.getDelassusCholeskyExpression().matrix();
    const auto & G = delassus_matrix_plain;
    auto G_expression = chol.getDelassusCholeskyExpression();
    //    std::cout << "G:\n" << delassus_matrix_plain << std::endl;

    Eigen::MatrixXd constraint_jacobian(delassus_matrix_plain.rows(), model.nv);
    constraint_jacobian.setZero();
    getConstraintsJacobian(model, data, constraint_models, constraint_datas, constraint_jacobian);

    const Eigen::VectorXd g = constraint_jacobian * v_free;
    //    std::cout << "g: " << g.transpose() << std::endl;

    ADMMContactSolverTpl<double> admm_solver(int(delassus_matrix_plain.rows()));
    admm_solver.setAbsolutePrecision(1e-10);
    admm_solver.setRelativePrecision(1e-14);

    PGSContactSolver pgs_solver(int(delassus_matrix_plain.rows()));
    pgs_solver.setAbsolutePrecision(1e-10);
    pgs_solver.setRelativePrecision(1e-14);

    const Eigen::VectorXd sol = G_expression.solve(-g);
    has_converged = admm_solver.solve(G_expression, g, constraint_models, primal_solution);
    primal_solution = admm_solver.getPrimalSolution();
    dual_solution = admm_solver.getDualSolution();
    //    std::cout << "constraint_velocity: " << constraint_velocity.transpose() << std::endl;

    const Eigen::VectorXd tau_ext = constraint_jacobian.transpose() * primal_solution / dt;

    v_next =
      v0
      + dt * aba(model, data, q0, v0, (tau0 + tau_ext).eval(), external_forces, Convention::WORLD);
  }

  Model model;
  Data data;
  std::vector<ConstraintModel> constraint_models;
  std::vector<ConstraintData> constraint_datas;
  std::vector<ConstraintSet> constraint_sets;
  Eigen::VectorXd v_next;

  Eigen::VectorXd primal_solution, dual_solution, dual_solution_sparse;
  bool has_converged;
};

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

Eigen::Vector3d computeFtot(const Eigen::VectorXd & contact_forces)
{
  Eigen::Vector3d f_tot = Eigen::Vector3d::Zero();
  for (int k = 0; k < contact_forces.size() / 3; ++k)
  {
    f_tot += contact_forces.segment(3 * k, 3);
  }
  return f_tot;
}

BOOST_AUTO_TEST_CASE(box)
{
  Model model;
  model.addJoint(0, JointModelFreeFlyer(), SE3::Identity(), "free_flyer");

  const int num_tests =
#ifdef NDEBUG
    100000
#else
    100
#endif
    ;

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  const double box_mass = 10;
  const Inertia box_inertia = Inertia::FromBox(box_mass, box_dims[0], box_dims[1], box_dims[2]);

  model.appendBodyToJoint(1, box_inertia);

  BOOST_CHECK(model.check(model.createData()));

  Eigen::VectorXd q0 = neutral(model);
  q0.const_cast_derived()[2] += box_dims[2] / 2;
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd tau0 = Eigen::VectorXd::Zero(model.nv);

  const double dt = 1e-3;

  typedef FrictionalPointConstraintModel ConstraintModel;
  typedef TestBoxTpl<ConstraintModel> TestBox;
  std::vector<ConstraintModel> constraint_models;

  const double friction_value = 0.4;
  {
    const SE3 local_placement_box(
      SE3::Matrix3::Identity(), 0.5 * SE3::Vector3(box_dims[0], box_dims[1], -box_dims[2]));
    SE3::Matrix3 rot = SE3::Matrix3::Identity();
    for (int i = 0; i < 4; ++i)
    {
      const SE3 local_placement(SE3::Matrix3::Identity(), rot * local_placement_box.translation());
      ConstraintModel cm(model, 0, SE3::Identity(), 1, local_placement);
      cm.set() = CoulombFrictionCone(friction_value);
      constraint_models.push_back(cm);
      rot = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitZ()).toRotationMatrix() * rot;
    }
  }

  // Test static motion with zero external force
  {
    const Force fext = Force::Zero();

    TestBox test(model, constraint_models);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.dual_solution.isZero(2e-10));
    const Force::Vector3 f_tot_ref = (-box_mass * Model::gravity981 - fext.linear()) * dt;
    BOOST_CHECK(computeFtot(test.primal_solution).isApprox(f_tot_ref, 1e-8));
    BOOST_CHECK(test.v_next.isZero(2e-10));
  }

  const double f_sliding = friction_value * Model::gravity981.norm() * box_mass;

  // Test static motion with small external force
  for (int k = 0; k < num_tests; ++k)
  {
    const double scaling = 0.9;
    Force fext = Force::Zero();
    fext.linear().head<2>().setRandom().normalize();
    fext.linear() *= scaling * f_sliding;

    TestBox test(model, constraint_models);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.dual_solution.isZero(1e-8));
    const Force::Vector3 f_tot_ref = (-box_mass * Model::gravity981 - fext.linear()) * dt;
    BOOST_CHECK(computeFtot(test.primal_solution).isApprox(f_tot_ref, 1e-6));
    BOOST_CHECK(test.v_next.isZero(1e-8));
  }

  // Test slidding motion
  for (int k = 0; k < num_tests; ++k)
  {
    const double scaling = 1.1;
    Force fext = Force::Zero();
    fext.linear().head<2>().setRandom().normalize();
    fext.linear() *= scaling * f_sliding;

    TestBox test(model, constraint_models);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    const Force::Vector3 f_tot_ref =
      (-box_mass * Model::gravity981 - 1 / scaling * fext.linear()) * dt;
    BOOST_CHECK(computeFtot(test.primal_solution).isApprox(f_tot_ref, 1e-6));
    BOOST_CHECK(
      math::fabs(Motion(test.v_next).linear().norm() - (f_sliding * 0.1 / box_mass * dt)) <= 1e-6);
    BOOST_CHECK(Motion(test.v_next).angular().isZero(1e-6));
  }
}

BOOST_AUTO_TEST_CASE(bilateral_box)
{
  Model model;
  model.addJoint(0, JointModelFreeFlyer(), SE3::Identity(), "free_flyer");

  const int num_tests =
#ifdef NDEBUG
    100000
#else
    100
#endif
    ;

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  const double box_mass = 10;
  const Inertia box_inertia = Inertia::FromBox(box_mass, box_dims[0], box_dims[1], box_dims[2]);

  model.appendBodyToJoint(1, box_inertia);

  BOOST_CHECK(model.check(model.createData()));

  Eigen::VectorXd q0 = neutral(model);
  q0.const_cast_derived()[2] += box_dims[2] / 2;
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd tau0 = Eigen::VectorXd::Zero(model.nv);

  const double dt = 1e-3;

  typedef BilateralPointConstraintModel ConstraintModel;
  typedef TestBoxTpl<ConstraintModel> TestBox;
  std::vector<ConstraintModel> constraint_models;

  {
    const SE3 local_placement_box(
      SE3::Matrix3::Identity(), 0.5 * SE3::Vector3(box_dims[0], box_dims[1], -box_dims[2]));
    SE3::Matrix3 rot = SE3::Matrix3::Identity();
    for (int i = 0; i < 4; ++i)
    {
      const SE3 local_placement(SE3::Matrix3::Identity(), rot * local_placement_box.translation());
      ConstraintModel cm(model, 0, SE3::Identity(), 1, local_placement);
      constraint_models.push_back(cm);
      rot = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitZ()).toRotationMatrix() * rot;
    }
  }

  // Test static motion with zero external force
  {
    const Force fext = Force::Zero();

    TestBox test(model, constraint_models);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.dual_solution.isZero(2e-10));
    BOOST_CHECK(
      computeFtot(test.primal_solution).isApprox(-box_mass * Model::gravity981 * dt, 1e-8));
    BOOST_CHECK(test.v_next.isZero(2e-10));
  }

  for (int k = 0; k < num_tests; ++k)
  {
    Force fext = Force::Zero();
    fext.linear().setRandom();

    TestBox test(model, constraint_models);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.dual_solution.isZero(1e-8));
    const Force::Vector3 f_tot_ref = (-box_mass * Model::gravity981 - fext.linear()) * dt;
    BOOST_CHECK(computeFtot(test.primal_solution).isApprox(f_tot_ref, 1e-6));
    BOOST_CHECK(test.v_next.isZero(1e-8));
  }
}

BOOST_AUTO_TEST_CASE(dry_friction_box)
{
  Model model;
  model.addJoint(0, JointModelFreeFlyer(), SE3::Identity(), "free_flyer");

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  const double box_mass = 10;
  const Inertia box_inertia = Inertia::FromBox(box_mass, box_dims[0], box_dims[1], box_dims[2]);

  model.appendBodyToJoint(1, box_inertia);
  model.gravity.setZero();
  Data data(model);

  Eigen::VectorXd q0 = neutral(model);
  q0.const_cast_derived()[2] += box_dims[2] / 2;
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd tau0 = Eigen::VectorXd::Zero(model.nv);

  const double dt = 1e-3;

  typedef FrictionalJointConstraintModel ConstraintModel;
  typedef ConstraintModel::ConstraintData ConstraintData;
  typedef ConstraintModel::ConstraintSet ConstraintSet;
  std::vector<ConstraintModel> constraint_models;
  std::vector<ConstraintData> constraint_datas;

  ConstraintModel dry_friction_free_flyer(model, ConstraintModel::JointIndexVector(1, 1));
  constraint_models.push_back(dry_friction_free_flyer);

  for (const auto & cm : constraint_models)
    constraint_datas.push_back(cm.createData());

  std::vector<ConstraintSet> constraint_sets;
  constraint_sets.push_back(
    BoxSet(Eigen::VectorXd::Constant(6, -1.), Eigen::VectorXd::Constant(6, +1.)));

  const auto & box_set = constraint_sets[0];
  constraint_models[0].set() = box_set;

  const Eigen::VectorXd v_free = dt * aba(model, data, q0, v0, tau0, Convention::WORLD);

  // Cholesky of the Delassus matrix
  crba(model, data, q0, Convention::WORLD);
  ContactCholeskyDecomposition chol(model, constraint_models);
  chol.compute(model, data, constraint_models, constraint_datas, 1e-10);

  auto G_expression = chol.getDelassusCholeskyExpression();
  const Eigen::MatrixXd delassus_matrix_plain = chol.getDelassusCholeskyExpression().matrix();
  const auto & G = delassus_matrix_plain;
  //    std::cout << "G:\n" << delassus_matrix_plain << std::endl;

  Eigen::MatrixXd constraint_jacobian(dry_friction_free_flyer.size(), model.nv);
  constraint_jacobian.setZero();
  getConstraintsJacobian(model, data, constraint_models, constraint_datas, constraint_jacobian);

  const Eigen::VectorXd g = constraint_jacobian * v_free;

  Eigen::VectorXd dual_solution = Eigen::VectorXd::Zero(g.size());
  Eigen::VectorXd primal_solution = Eigen::VectorXd::Zero(g.size());
  ADMMContactSolver admm_solver(int(delassus_matrix_plain.rows()));
  admm_solver.setAbsolutePrecision(1e-10);
  admm_solver.setRelativePrecision(1e-14);

  const bool has_converged = admm_solver.solve(G_expression, g, constraint_models, dual_solution);
  BOOST_CHECK(has_converged);

  primal_solution = G * dual_solution + g;

  BOOST_CHECK(std::fabs(primal_solution.dot(dual_solution)) <= 1e-8);
  BOOST_CHECK(dual_solution.isZero());

  typedef TestBoxTpl<ConstraintModel> TestBox;
  // Test static motion with zero external force
  {
    TestBox test(model, constraint_models);
    test(q0, v0, tau0, Force::Zero(), dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.primal_solution.isZero(2e-10));
    BOOST_CHECK(test.v_next.isZero(2e-10));
    BOOST_CHECK(box_set.isInside(test.dual_solution));
  }

  for (int i = 0; i < 6; ++i)
  {
    TestBox test(model, constraint_models);
    test(q0, v0, tau0 + 2 * Force::Vector6::Unit(i) / dt, Force::Zero(), dt);

    //    std::cout << "test.dual_solution: " << test.dual_solution.transpose() << std::endl;
    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(!test.primal_solution.isZero(2e-10));
    BOOST_CHECK(!test.v_next.isZero(2e-10));
    BOOST_CHECK(box_set.isInside(test.dual_solution));
    std::cout << "dual_solution[i]: " << dual_solution[i] << std::endl;
    std::cout << "box_set.lb()[i]: " << box_set.lb()[i] << std::endl;
    BOOST_CHECK(std::fabs(test.dual_solution[i] - box_set.lb()[i]) < 1e-8);
  }

  // Sign reversed
  for (int i = 0; i < 6; ++i)
  {
    TestBox test(model, constraint_models);
    test(q0, v0, tau0 - 2 * Force::Vector6::Unit(i) / dt, Force::Zero(), dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(!test.primal_solution.isZero(2e-10));
    BOOST_CHECK(!test.v_next.isZero(2e-10));
    BOOST_CHECK(box_set.isInside(test.dual_solution));
    std::cout << "dual_solution[i]: " << dual_solution[i] << std::endl;
    std::cout << "box_set.ub()[i]: " << box_set.ub()[i] << std::endl;
    BOOST_CHECK(std::fabs(test.dual_solution[i] - box_set.ub()[i]) < 1e-8);
  }
}

BOOST_AUTO_TEST_SUITE_END()
