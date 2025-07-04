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
    const double dt,
    const bool test_warmstart = false)
  {
    std::vector<Force> external_forces(size_t(model.njoints), Force::Zero());
    external_forces[1] = fext;

    const Eigen::VectorXd v_free =
      v0 + dt * aba(model, data, q0, v0, tau0, external_forces, Convention::WORLD);

    // Cholesky of the Delassus matrix
    crba(model, data, q0, Convention::WORLD);
    ContactCholeskyDecomposition chol(model, constraint_models);
    chol.resize(model, constraint_models);
    chol.compute(model, data, constraint_models, constraint_datas, 1e-10);
    // std::cout << "chol.getDamping() :   " << chol.getDamping().transpose() << std::endl;

    const Eigen::MatrixXd delassus_matrix_plain = chol.getDelassusCholeskyExpression().matrix();
    std::cout << "delassus_matrix_plain:    " << delassus_matrix_plain << std::endl;
    auto G_expression = chol.getDelassusCholeskyExpression();

    Eigen::MatrixXd constraint_jacobian(delassus_matrix_plain.rows(), model.nv);
    constraint_jacobian.setZero();
    getConstraintsJacobian(model, data, constraint_models, constraint_datas, constraint_jacobian);

    Eigen::VectorXd time_scaling(delassus_matrix_plain.rows());
    internal::getTimeScalingFromAccelerationToConstraints(constraint_models, dt, time_scaling);
    const Eigen::VectorXd g = (constraint_jacobian * v_free).cwiseProduct(time_scaling / dt);

    ADMMContactSolverTpl<double> admm_solver(int(delassus_matrix_plain.rows()));
    admm_solver.setAbsolutePrecision(1e-13);
    admm_solver.setRelativePrecision(1e-14);
    admm_solver.setLanczosSize((int)g.size());

    PGSContactSolver pgs_solver(int(delassus_matrix_plain.rows()));
    pgs_solver.setAbsolutePrecision(1e-10);
    pgs_solver.setRelativePrecision(1e-14);

    Eigen::VectorXd compliance = Eigen::VectorXd::Zero(g.size());
    G_expression.updateCompliance(compliance);
    Eigen::VectorXd mean_inertia = Eigen::VectorXd::Constant(g.size(), data.M.diagonal().trace());
    mean_inertia /= (double)(g.size());
    mean_inertia = mean_inertia.array().sqrt();
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> preconditioner_vec(mean_inertia);
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution_constref(primal_solution);
    has_converged = admm_solver.solve(
      G_expression, g, constraint_models, dt, preconditioner_vec, primal_solution_constref);
    primal_solution = admm_solver.getPrimalSolution();
    // std::cout << "primal_solution :   " << (primal_solution).transpose() << std::endl;
    // std::cout << "g :   " << (g).transpose() << std::endl;
    // std::cout << "delassus_matrix_plain :   " << (delassus_matrix_plain).transpose() <<
    // std::endl; std::cout << "delassus_matrix_plain * primal_solution :   " <<
    // (delassus_matrix_plain * primal_solution).transpose() << std::endl; std::cout <<
    // "time_scaling.asDiagonal()* delassus_matrix_plain * primal_solution :   " <<
    // (time_scaling.asDiagonal()*delassus_matrix_plain * primal_solution).transpose() << std::endl;
    // std::cout << "time_scaling.asDiagonal()* delassus_matrix_plain * primal_solution + g :   " <<
    // (time_scaling.asDiagonal()*delassus_matrix_plain * primal_solution + g).transpose() <<
    // std::endl;

    if (test_warmstart)
    {
      boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution_warmstart(primal_solution);
      has_converged =
        has_converged
        && admm_solver.solve(
          G_expression, g, constraint_models, dt, preconditioner_vec, primal_solution_warmstart);
      primal_solution = admm_solver.getPrimalSolution();
    }

    dual_solution = admm_solver.getDualSolution();
    n_iter = admm_solver.getIterationCount();
    const Eigen::VectorXd tau_ext = constraint_jacobian.transpose() * primal_solution;

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
  int n_iter;
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

BOOST_AUTO_TEST_CASE(ball)
{
  Model model;
  model.addJoint(0, JointModelFreeFlyer(), SE3::Identity(), "free_flyer");

  const double ball_dim = 1;
  const double ball_mass = 10;
  const Inertia ball_inertia = Inertia::FromSphere(ball_mass, ball_dim);

  model.appendBodyToJoint(1, ball_inertia);

  BOOST_CHECK(model.check(model.createData()));

  Eigen::VectorXd q0 = neutral(model);
  q0.const_cast_derived()[2] += ball_dim / 2;
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd tau0 = Eigen::VectorXd::Zero(model.nv);

  const double dt = 1e-3;

  typedef FrictionalPointConstraintModel ConstraintModel;
  typedef TestBoxTpl<ConstraintModel> TestBox;
  std::vector<ConstraintModel> constraint_models;

  const double friction_value = 0.4;
  {
    const SE3 local_placement_ball(SE3::Matrix3::Identity(), SE3::Vector3(0, 0, -ball_dim));
    ConstraintModel cm(model, 0, SE3::Identity(), 1, local_placement_ball);
    cm.set() = CoulombFrictionCone(friction_value);
    constraint_models.push_back(cm);
  }

  // Test static motion with zero external force
  {
    const Force fext = Force::Zero();

    TestBox test(model, constraint_models);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.dual_solution.isZero(2e-10));
    const Force::Vector3 f_tot_ref = -ball_mass * Model::gravity981 - fext.linear();
    BOOST_CHECK(computeFtot(test.primal_solution).isApprox(f_tot_ref, 1e-8));
    BOOST_CHECK(test.v_next.isZero(2e-10));

    // Test warmstart
    test(q0, v0, tau0, fext, dt, true);
    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.dual_solution.isZero(2e-10));
    BOOST_CHECK(computeFtot(test.primal_solution).isApprox(f_tot_ref, 1e-8));
    BOOST_CHECK(test.v_next.isZero(2e-10));
    BOOST_CHECK(test.n_iter == 0);
  }
}

void buildStackOfCubesModel(
  std::vector<double> masses,
  ::pinocchio::Model & model,
  std::vector<FrictionalPointConstraintModel> & constraint_models)
{
  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  const int n_cubes = (int)masses.size();

  for (int i = 0; i < n_cubes; i++)
  {
    const double box_mass = masses[(std::size_t)i];
    const Inertia box_inertia = Inertia::FromBox(box_mass, box_dims[0], box_dims[1], box_dims[2]);
    JointIndex joint_id =
      model.addJoint(0, JointModelFreeFlyer(), SE3::Identity(), "free_flyer_" + std::to_string(i));
    model.appendBodyToJoint(joint_id, box_inertia);
  }

  const double friction_value = 0.4;
  for (int i = 0; i < n_cubes; i++)
  {
    const SE3 local_placement_box_1(
      SE3::Matrix3::Identity(), 0.5 * SE3::Vector3(box_dims[0], box_dims[1], box_dims[2]));
    const SE3 local_placement_box_2(
      SE3::Matrix3::Identity(), 0.5 * SE3::Vector3(box_dims[0], box_dims[1], -box_dims[2]));
    SE3::Matrix3 rot = SE3::Matrix3::Identity();
    for (int j = 0; j < 4; ++j)
    {
      const SE3 local_placement_1(
        SE3::Matrix3::Identity(), rot * local_placement_box_1.translation());
      const SE3 local_placement_2(
        SE3::Matrix3::Identity(), rot * local_placement_box_2.translation());
      FrictionalPointConstraintModel cm(
        model, (JointIndex)i, local_placement_1, (JointIndex)i + 1, local_placement_2);
      cm.set() = CoulombFrictionCone(friction_value);
      constraint_models.push_back(cm);
      rot = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitZ()).toRotationMatrix() * rot;
    }
  }
}

BOOST_AUTO_TEST_CASE(box)
{
  Model model;
  typedef FrictionalPointConstraintModel ConstraintModel;
  typedef TestBoxTpl<ConstraintModel> TestBox;
  std::vector<ConstraintModel> constraint_models;
  const double box_mass = 1e1;
  const std::vector<double> masses = {box_mass};

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  buildStackOfCubesModel(masses, model, constraint_models);

  const int num_tests =
#ifdef NDEBUG
    100000
#else
    100
#endif
    ;

  BOOST_CHECK(model.check(model.createData()));

  Eigen::VectorXd q0 = neutral(model);
  q0[2] += box_dims[2] / 2;
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd tau0 = Eigen::VectorXd::Zero(model.nv);

  const double dt = 1e-3;

  // Test static motion with zero external force
  {
    const Force fext = Force::Zero();

    TestBox test(model, constraint_models);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.dual_solution.isZero(2e-10));
    const Force::Vector3 f_tot_ref = -box_mass * Model::gravity981 - fext.linear();
    BOOST_CHECK(computeFtot(test.primal_solution).isApprox(f_tot_ref, 1e-8));
    std::cout << "test.primal_solution: " << test.primal_solution.transpose() << std::endl;
    std::cout << "computeFtot(test.primal_solution): "
              << computeFtot(test.primal_solution).transpose() << std::endl;
    std::cout << "f_tot_ref: " << f_tot_ref.transpose() << std::endl;
    BOOST_CHECK(test.v_next.isZero(2e-10));
  }

  const double friction_value = 0.4;
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
    const Force::Vector3 f_tot_ref = -box_mass * Model::gravity981 - fext.linear();
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
    const Force::Vector3 f_tot_ref = -box_mass * Model::gravity981 - 1 / scaling * fext.linear();
    BOOST_CHECK(computeFtot(test.primal_solution).isApprox(f_tot_ref, 1e-6));
    BOOST_CHECK(
      math::fabs(Motion(test.v_next).linear().norm() - (f_sliding * 0.1 / box_mass * dt)) <= 1e-6);
    BOOST_CHECK(Motion(test.v_next).angular().isZero(1e-6));
  }
}

BOOST_AUTO_TEST_CASE(stack_of_boxes)
{
  const int n_cubes = 10;
  const double conditionning = 1e6;
  const double mass_factor = std::pow(conditionning, 1. / (n_cubes - 1));
  std::vector<double> masses;
  double mass_tot = 0;
  for (int i = 0; i < n_cubes; i++)
  {
    const double box_mass = 1e-3 * std::pow(mass_factor, i);
    masses.push_back(box_mass);
    mass_tot += box_mass;
  }

  Model model;
  typedef FrictionalPointConstraintModel ConstraintModel;
  typedef TestBoxTpl<ConstraintModel> TestBox;
  std::vector<ConstraintModel> constraint_models;

  buildStackOfCubesModel(masses, model, constraint_models);
  BOOST_CHECK(model.check(model.createData()));

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();

  Eigen::VectorXd q0 = neutral(model);
  for (int i = 0; i < n_cubes; i++)
  {
    q0[7 * i + 2] = i * box_dims[2];
    q0[7 * i + 2] += box_dims[2] / 2;
  }
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd tau0 = Eigen::VectorXd::Zero(model.nv);

  const double dt = 1e-3;

  // Test static motion with zero external force
  {
    const Force fext = Force::Zero();

    TestBox test(model, constraint_models);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.dual_solution.isZero(2e-10));
    const Force::Vector3 f_tot_ref = -mass_tot * Model::gravity981;
    BOOST_CHECK(computeFtot(test.primal_solution).head(3).isApprox(f_tot_ref, 1e-8));
    BOOST_CHECK(test.v_next.isZero(2e-10));
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
    BOOST_CHECK(computeFtot(test.primal_solution).isApprox(-box_mass * Model::gravity981, 1e-8));
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
    const Force::Vector3 f_tot_ref = -box_mass * Model::gravity981 - fext.linear();
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

  const Eigen::VectorXd v_free = v0 + dt * aba(model, data, q0, v0, tau0, Convention::WORLD);

  // Cholesky of the Delassus matrix
  crba(model, data, q0, Convention::WORLD);
  ContactCholeskyDecomposition chol(model, constraint_models);
  chol.resize(model, constraint_models);
  chol.compute(model, data, constraint_models, constraint_datas, 1e-10);

  auto G_expression = chol.getDelassusCholeskyExpression();
  const Eigen::MatrixXd delassus_matrix_plain = chol.getDelassusCholeskyExpression().matrix();
  const auto & G = delassus_matrix_plain;
  //    std::cout << "G:\n" << delassus_matrix_plain << std::endl;

  Eigen::MatrixXd constraint_jacobian(dry_friction_free_flyer.size(), model.nv);
  constraint_jacobian.setZero();
  getConstraintsJacobian(model, data, constraint_models, constraint_datas, constraint_jacobian);

  Eigen::VectorXd time_scaling(delassus_matrix_plain.rows());
  internal::getTimeScalingFromAccelerationToConstraints(constraint_models, dt, time_scaling);
  const Eigen::VectorXd g = (constraint_jacobian * v_free).cwiseProduct(time_scaling / dt);

  Eigen::VectorXd compliance = Eigen::VectorXd::Zero(g.size());
  G_expression.updateCompliance(compliance);
  boost::optional<Eigen::Ref<const Eigen::VectorXd>> preconditioner_vec(
    Eigen::VectorXd::Ones(g.size()));
  Eigen::VectorXd dual_solution(Eigen::VectorXd::Zero(g.size()));
  boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution(
    Eigen::VectorXd::Zero(g.size()));
  ADMMContactSolver admm_solver(int(delassus_matrix_plain.rows()));
  admm_solver.setAbsolutePrecision(1e-13);
  admm_solver.setRelativePrecision(1e-14);

  const bool has_converged =
    admm_solver.solve(G_expression, g, constraint_models, dt, preconditioner_vec, primal_solution);
  BOOST_CHECK(has_converged);

  dual_solution = (G * primal_solution.get()).cwiseProduct(time_scaling) + g;

  BOOST_CHECK(std::fabs(primal_solution.get().dot(dual_solution)) <= 1e-8);
  BOOST_CHECK(primal_solution.get().isZero());

  typedef TestBoxTpl<ConstraintModel> TestBox;

  // Test static motion with zero external force
  {
    TestBox test(model, constraint_models);
    test(q0, v0, tau0, Force::Zero(), dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.dual_solution.isZero(2e-10));
    BOOST_CHECK(test.v_next.isZero(2e-10));
    BOOST_CHECK(box_set.isInside(test.primal_solution));
  }

  for (int i = 0; i < 6; ++i)
  {
    TestBox test(model, constraint_models);
    test(q0, v0, tau0 + 2 * Force::Vector6::Unit(i) / dt, Force::Zero(), dt);

    //    std::cout << "test.dual_solution: " << test.dual_solution.transpose() << std::endl;
    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(!test.primal_solution.isZero(2e-10));
    BOOST_CHECK(!test.v_next.isZero(2e-10));
    BOOST_CHECK(box_set.isInside(test.primal_solution));
    BOOST_CHECK(std::fabs(test.primal_solution[i] - box_set.lb()[i]) < 1e-8);
  }

  // Sign reversed
  for (int i = 0; i < 6; ++i)
  {
    TestBox test(model, constraint_models);
    test(q0, v0, tau0 - 2 * Force::Vector6::Unit(i) / dt, Force::Zero(), dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(!test.dual_solution.isZero(2e-10));
    BOOST_CHECK(!test.v_next.isZero(2e-10));
    BOOST_CHECK(box_set.isInside(test.primal_solution));
    BOOST_CHECK(std::fabs(test.primal_solution[i] - box_set.ub()[i]) < 1e-8);
  }
}

BOOST_AUTO_TEST_CASE(joint_limit_slider)
{
  Model model;
  model.addJoint(0, JointModelPX(), SE3::Identity(), "slider");
  model.lowerPositionLimit[0] = 0.;

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  const double box_mass = 10;
  const Inertia box_inertia = Inertia::FromBox(box_mass, box_dims[0], box_dims[1], box_dims[2]);

  model.appendBodyToJoint(1, box_inertia);
  model.gravity.setZero();
  Data data(model);

  Eigen::VectorXd q0 = Eigen::VectorXd::Zero(model.nq);
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd tau_push_against_lower_bound = -Eigen::VectorXd::Ones(model.nv);

  const double dt = 1e-3;

  typedef JointLimitConstraintModel ConstraintModel;
  typedef ConstraintModel::ConstraintData ConstraintData;
  std::vector<ConstraintModel> constraint_models;
  std::vector<ConstraintData> constraint_datas;

  ConstraintModel joint_limit_constraint_model(model, ConstraintModel::JointIndexVector(1, 1));
  constraint_models.push_back(joint_limit_constraint_model);

  for (const auto & cm : constraint_models)
    constraint_datas.push_back(cm.createData());

  const Eigen::VectorXd v_free_against_lower_bound =
    v0 + dt * aba(model, data, q0, v0, tau_push_against_lower_bound, Convention::WORLD);
  const Eigen::VectorXd v_free_move_away =
    v0 + dt * aba(model, data, q0, v0, -tau_push_against_lower_bound, Convention::WORLD);

  // Cholesky of the Delassus matrix
  crba(model, data, q0, Convention::WORLD);

  data.q_in = q0;
  auto & cmodel = constraint_models[0];
  auto & cdata = constraint_datas[0];
  cmodel.resize(model, data, cdata);
  cmodel.calc(model, data, cdata);
  ContactCholeskyDecomposition chol(model, constraint_models);
  chol.resize(model, constraint_models);
  chol.compute(model, data, constraint_models, constraint_datas, 1e-10);

  auto G_expression = chol.getDelassusCholeskyExpression();
  const auto G_plain = G_expression.matrix();
  const Eigen::MatrixXd delassus_matrix_plain = G_expression.matrix();

  Eigen::MatrixXd constraint_jacobian(cmodel.activeSize(), model.nv);
  constraint_jacobian.setZero();
  getConstraintsJacobian(model, data, constraint_models, constraint_datas, constraint_jacobian);

  Eigen::VectorXd time_scaling(delassus_matrix_plain.rows());
  internal::getTimeScalingFromAccelerationToConstraints(constraint_models, dt, time_scaling);

  // External torques push the slider against the lower bound
  {
    const Eigen::VectorXd g_against_lower_bound = constraint_jacobian * v_free_against_lower_bound;
    const Eigen::VectorXd g_tilde_against_lower_bound =
      g_against_lower_bound + cdata.constraint_residual / dt;
    g_against_lower_bound.const_cast_derived().array() *= time_scaling.array() / dt;
    g_tilde_against_lower_bound.const_cast_derived().array() *= time_scaling.array() / dt;

    Eigen::VectorXd dual_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    Eigen::VectorXd primal_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    ADMMContactSolver admm_solver(int(delassus_matrix_plain.rows()));
    admm_solver.setAbsolutePrecision(1e-13);
    admm_solver.setRelativePrecision(1e-14);

    Eigen::VectorXd compliance = Eigen::VectorXd::Zero(g_tilde_against_lower_bound.size());
    G_expression.updateCompliance(compliance);
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> preconditioner_vec(
      Eigen::VectorXd::Ones(g_tilde_against_lower_bound.size()));
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution_constref(primal_solution);
    const bool has_converged = admm_solver.solve(
      G_expression, g_tilde_against_lower_bound, constraint_models, dt, preconditioner_vec,
      primal_solution_constref);
    primal_solution = admm_solver.getPrimalSolution();
    BOOST_CHECK(has_converged);

    dual_solution = (G_plain * primal_solution).cwiseProduct(time_scaling) + g_against_lower_bound;
    Eigen::VectorXd dual_solution2 = admm_solver.getDualSolution();

    BOOST_CHECK(std::fabs(primal_solution.dot(dual_solution)) <= 1e-8);
    BOOST_CHECK(dual_solution.isZero(1e-6));
    BOOST_CHECK(dual_solution2.isZero(1e-6));

    BOOST_CHECK((tau_push_against_lower_bound + constraint_jacobian.transpose() * primal_solution)
                  .isZero(1e-6));
  }

  // External torques push the slider away from the lower bound
  {
    const Eigen::VectorXd g_move_away = constraint_jacobian * v_free_move_away;
    const Eigen::VectorXd g_tilde_move_away = g_move_away + cdata.constraint_residual / dt;
    g_move_away.const_cast_derived().array() *= time_scaling.array() / dt;
    g_tilde_move_away.const_cast_derived().array() *= time_scaling.array() / dt;

    Eigen::VectorXd dual_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    Eigen::VectorXd primal_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    ADMMContactSolver admm_solver(int(delassus_matrix_plain.rows()));
    admm_solver.setAbsolutePrecision(1e-13);
    admm_solver.setRelativePrecision(1e-14);

    Eigen::VectorXd compliance = Eigen::VectorXd::Zero(g_tilde_move_away.size());
    G_expression.updateCompliance(compliance);
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> preconditioner_vec(
      Eigen::VectorXd::Ones(g_tilde_move_away.size()));
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution_constref(primal_solution);
    const bool has_converged = admm_solver.solve(
      G_expression, g_tilde_move_away, constraint_models, dt, preconditioner_vec,
      primal_solution_constref);
    primal_solution = admm_solver.getPrimalSolution();
    BOOST_CHECK(has_converged);

    dual_solution = (G_plain * primal_solution).cwiseProduct(time_scaling) + g_move_away;
    Eigen::VectorXd dual_solution2 = admm_solver.getDualSolution();

    BOOST_CHECK(std::fabs(primal_solution.dot(dual_solution)) <= 1e-8);
    BOOST_CHECK(primal_solution.isZero());
    BOOST_CHECK(dual_solution.isApprox(g_move_away));
  }
}

BOOST_AUTO_TEST_CASE(joint_limit_revolute_xyz)
{
  Model model;
  JointIndex joint_id_x = model.addJoint(0, JointModelRX(), SE3::Identity(), "revolute_x");
  JointIndex joint_id_y = model.addJoint(joint_id_x, JointModelRY(), SE3::Identity(), "revolute_y");
  JointIndex joint_id_z = model.addJoint(joint_id_y, JointModelRZ(), SE3::Identity(), "revolute_z");

  const SE3::Vector3 small_box_dims = SE3::Vector3::Ones() * 1e-3;
  const double small_box_mass = 1e-6;
  const Inertia small_box_inertia =
    Inertia::FromBox(small_box_mass, small_box_dims[0], small_box_dims[1], small_box_dims[2]);

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  const double box_mass = 10;
  const Inertia box_inertia = Inertia::FromBox(box_mass, box_dims[0], box_dims[1], box_dims[2]);

  model.appendBodyToJoint(joint_id_x, small_box_inertia);
  model.appendBodyToJoint(joint_id_y, small_box_inertia);
  model.appendBodyToJoint(joint_id_z, box_inertia);
  model.gravity.setZero();
  model.lowerPositionLimit[0] = 0.;
  model.lowerPositionLimit[1] = 0.;
  model.lowerPositionLimit[2] = 0.;
  Data data(model);

  Eigen::VectorXd q0 = Eigen::VectorXd::Zero(model.nq);
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd tau_push_against_lower_bound = -Eigen::VectorXd::Ones(model.nv);

  const double dt = 1e-3;

  typedef JointLimitConstraintModel ConstraintModel;
  typedef ConstraintModel::ConstraintData ConstraintData;
  std::vector<ConstraintModel> constraint_models;
  std::vector<ConstraintData> constraint_datas;
  ConstraintModel::JointIndexVector active_joints = {joint_id_x, joint_id_y, joint_id_z};

  ConstraintModel joint_limit_constraint_model(model, active_joints);
  constraint_models.push_back(joint_limit_constraint_model);

  for (const auto & cm : constraint_models)
    constraint_datas.push_back(cm.createData());

  const Eigen::VectorXd v_free_against_lower_bound =
    v0 + dt * aba(model, data, q0, v0, tau_push_against_lower_bound, Convention::WORLD);
  const Eigen::VectorXd v_free_move_away =
    v0 + dt * aba(model, data, q0, v0, -tau_push_against_lower_bound, Convention::WORLD);

  // Cholesky of the Delassus matrix
  crba(model, data, q0, Convention::WORLD);

  data.q_in = q0;
  auto & cmodel = constraint_models[0];
  auto & cdata = constraint_datas[0];
  cmodel.resize(model, data, cdata);
  cmodel.calc(model, data, cdata);
  ContactCholeskyDecomposition chol(model, constraint_models);
  chol.resize(model, constraint_models);
  chol.compute(model, data, constraint_models, constraint_datas, 1e-10);

  auto G_expression = chol.getDelassusCholeskyExpression();
  const auto G_plain = G_expression.matrix();
  const Eigen::MatrixXd delassus_matrix_plain = G_expression.matrix();

  Eigen::MatrixXd constraint_jacobian(cmodel.activeSize(), model.nv);
  constraint_jacobian.setZero();
  getConstraintsJacobian(model, data, constraint_models, constraint_datas, constraint_jacobian);

  Eigen::VectorXd time_scaling(delassus_matrix_plain.rows());
  internal::getTimeScalingFromAccelerationToConstraints(constraint_models, dt, time_scaling);

  // External torques push the slider against the lower bound
  {
    const Eigen::VectorXd g_against_lower_bound = constraint_jacobian * v_free_against_lower_bound;
    const Eigen::VectorXd g_tilde_against_lower_bound =
      g_against_lower_bound + cdata.constraint_residual / dt;
    g_against_lower_bound.const_cast_derived().array() *= time_scaling.array() / dt;
    g_tilde_against_lower_bound.const_cast_derived().array() *= time_scaling.array() / dt;

    Eigen::VectorXd dual_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    Eigen::VectorXd primal_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    ADMMContactSolver admm_solver(int(delassus_matrix_plain.rows()));
    admm_solver.setAbsolutePrecision(1e-13);
    admm_solver.setRelativePrecision(1e-14);

    Eigen::VectorXd compliance = Eigen::VectorXd::Zero(g_tilde_against_lower_bound.size());
    G_expression.updateCompliance(compliance);
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> preconditioner_vec(
      Eigen::VectorXd::Ones(g_tilde_against_lower_bound.size()));
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution_constref(primal_solution);
    const bool has_converged = admm_solver.solve(
      G_expression, g_tilde_against_lower_bound, constraint_models, dt, preconditioner_vec,
      primal_solution_constref);
    primal_solution = admm_solver.getPrimalSolution();
    BOOST_CHECK(has_converged);

    dual_solution = (G_plain * primal_solution).cwiseProduct(time_scaling) + g_against_lower_bound;
    Eigen::VectorXd dual_solution2 = admm_solver.getDualSolution();

    BOOST_CHECK(std::fabs(primal_solution.dot(dual_solution)) <= 1e-8);
    BOOST_CHECK(dual_solution.isZero(1e-6));
    BOOST_CHECK(dual_solution2.isZero(1e-6));

    // std::cout << "tau_push_against_lower_bound:   " << tau_push_against_lower_bound << std::endl;
    // std::cout << "constraint_jacobian.transpose() * primal_solution:   "
    //           << constraint_jacobian.transpose() * primal_solution << std::endl;
    // std::cout << "primal_solution:   " << primal_solution << std::endl;

    BOOST_CHECK((tau_push_against_lower_bound + constraint_jacobian.transpose() * primal_solution)
                  .isZero(1e-6));
  }

  // External torques push the slider away from the lower bound
  {
    const Eigen::VectorXd g_move_away = constraint_jacobian * v_free_move_away;
    const Eigen::VectorXd g_tilde_move_away = g_move_away + cdata.constraint_residual / dt;
    g_move_away.const_cast_derived().array() *= time_scaling.array() / dt;
    g_tilde_move_away.const_cast_derived().array() *= time_scaling.array() / dt;

    Eigen::VectorXd dual_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    Eigen::VectorXd primal_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    ADMMContactSolver admm_solver(int(delassus_matrix_plain.rows()));
    admm_solver.setAbsolutePrecision(1e-13);
    admm_solver.setRelativePrecision(1e-14);

    Eigen::VectorXd compliance = Eigen::VectorXd::Zero(g_tilde_move_away.size());
    G_expression.updateCompliance(compliance);
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> preconditioner_vec(
      Eigen::VectorXd::Ones(g_tilde_move_away.size()));
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution_constref(primal_solution);
    const bool has_converged = admm_solver.solve(
      G_expression, g_tilde_move_away, constraint_models, dt, preconditioner_vec,
      primal_solution_constref);
    primal_solution = admm_solver.getPrimalSolution();
    BOOST_CHECK(has_converged);

    dual_solution = (G_plain * primal_solution).cwiseProduct(time_scaling) + g_move_away;
    Eigen::VectorXd dual_solution2 = admm_solver.getDualSolution();

    BOOST_CHECK(std::fabs(primal_solution.dot(dual_solution)) <= 1e-8);
    BOOST_CHECK(primal_solution.isZero());
    BOOST_CHECK(dual_solution.isApprox(g_move_away));
  }
}

BOOST_AUTO_TEST_CASE(joint_limit_slider_xyz)
{
  Model model;
  JointIndex joint_id_x = model.addJoint(0, JointModelPX(), SE3::Identity(), "slider_x");
  JointIndex joint_id_y = model.addJoint(joint_id_x, JointModelPY(), SE3::Identity(), "slider_y");
  JointIndex joint_id_z = model.addJoint(joint_id_y, JointModelPZ(), SE3::Identity(), "slider_z");

  const SE3::Vector3 small_box_dims = SE3::Vector3::Ones() * 1e-3;
  const double small_box_mass = 1e-6;
  const Inertia small_box_inertia =
    Inertia::FromBox(small_box_mass, small_box_dims[0], small_box_dims[1], small_box_dims[2]);

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  const double box_mass = 10;
  const Inertia box_inertia = Inertia::FromBox(box_mass, box_dims[0], box_dims[1], box_dims[2]);

  model.appendBodyToJoint(joint_id_x, small_box_inertia);
  model.appendBodyToJoint(joint_id_y, small_box_inertia);
  model.appendBodyToJoint(joint_id_z, box_inertia);
  model.gravity.setZero();
  model.lowerPositionLimit[0] = 0.;
  model.lowerPositionLimit[1] = 0.;
  model.lowerPositionLimit[2] = 0.;
  Data data(model);

  Eigen::VectorXd q0 = Eigen::VectorXd::Zero(model.nq);
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd tau_push_against_lower_bound = -Eigen::VectorXd::Ones(model.nv);

  const double dt = 1e-3;

  typedef JointLimitConstraintModel ConstraintModel;
  typedef ConstraintModel::ConstraintData ConstraintData;
  std::vector<ConstraintModel> constraint_models;
  std::vector<ConstraintData> constraint_datas;
  ConstraintModel::JointIndexVector active_joints = {joint_id_x, joint_id_y, joint_id_z};

  ConstraintModel joint_limit_constraint_model(model, active_joints);
  constraint_models.push_back(joint_limit_constraint_model);

  for (const auto & cm : constraint_models)
    constraint_datas.push_back(cm.createData());

  const Eigen::VectorXd v_free_against_lower_bound =
    v0 + dt * aba(model, data, q0, v0, tau_push_against_lower_bound, Convention::WORLD);
  const Eigen::VectorXd v_free_move_away =
    v0 + dt * aba(model, data, q0, v0, -tau_push_against_lower_bound, Convention::WORLD);

  // Cholesky of the Delassus matrix
  crba(model, data, q0, Convention::WORLD);

  data.q_in = q0;
  auto & cmodel = constraint_models[0];
  auto & cdata = constraint_datas[0];
  cmodel.resize(model, data, cdata);
  cmodel.calc(model, data, cdata);
  ContactCholeskyDecomposition chol(model, constraint_models);
  chol.resize(model, constraint_models);
  chol.compute(model, data, constraint_models, constraint_datas, 1e-10);

  auto G_expression = chol.getDelassusCholeskyExpression();
  const auto G_plain = G_expression.matrix();
  const Eigen::MatrixXd delassus_matrix_plain = G_expression.matrix();

  Eigen::MatrixXd constraint_jacobian(cmodel.activeSize(), model.nv);
  constraint_jacobian.setZero();
  getConstraintsJacobian(model, data, constraint_models, constraint_datas, constraint_jacobian);

  Eigen::VectorXd time_scaling(delassus_matrix_plain.rows());
  internal::getTimeScalingFromAccelerationToConstraints(constraint_models, dt, time_scaling);

  // External torques push the slider against the lower bound
  {
    const Eigen::VectorXd g_against_lower_bound = constraint_jacobian * v_free_against_lower_bound;
    const Eigen::VectorXd g_tilde_against_lower_bound =
      g_against_lower_bound + cdata.constraint_residual / dt;
    g_against_lower_bound.const_cast_derived().array() *= time_scaling.array() / dt;
    g_tilde_against_lower_bound.const_cast_derived().array() *= time_scaling.array() / dt;

    Eigen::VectorXd dual_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    Eigen::VectorXd primal_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    ADMMContactSolver admm_solver(int(delassus_matrix_plain.rows()));
    admm_solver.setAbsolutePrecision(1e-13);
    admm_solver.setRelativePrecision(1e-14);

    Eigen::VectorXd compliance = Eigen::VectorXd::Zero(g_tilde_against_lower_bound.size());
    G_expression.updateCompliance(compliance);
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> preconditioner_vec(
      Eigen::VectorXd::Ones(g_tilde_against_lower_bound.size()));
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution_constref(primal_solution);
    const bool has_converged = admm_solver.solve(
      G_expression, g_tilde_against_lower_bound, constraint_models, dt, preconditioner_vec,
      primal_solution_constref);
    primal_solution = admm_solver.getPrimalSolution();
    BOOST_CHECK(has_converged);

    dual_solution = (G_plain * primal_solution).cwiseProduct(time_scaling) + g_against_lower_bound;
    Eigen::VectorXd dual_solution2 = admm_solver.getDualSolution();

    BOOST_CHECK(std::fabs(primal_solution.dot(dual_solution)) <= 1e-8);
    BOOST_CHECK(dual_solution.isZero(1e-6));
    BOOST_CHECK(dual_solution2.isZero(1e-6));

    // std::cout << "tau_push_against_lower_bound:   " << tau_push_against_lower_bound << std::endl;
    // std::cout << "constraint_jacobian.transpose() * primal_solution:   "
    //           << constraint_jacobian.transpose() * primal_solution << std::endl;
    // std::cout << "primal_solution:   " << primal_solution << std::endl;

    BOOST_CHECK((tau_push_against_lower_bound + constraint_jacobian.transpose() * primal_solution)
                  .isZero(1e-6));
  }

  // External torques push the slider away from the lower bound
  {
    const Eigen::VectorXd g_move_away = constraint_jacobian * v_free_move_away;
    const Eigen::VectorXd g_tilde_move_away = g_move_away + cdata.constraint_residual / dt;
    g_move_away.const_cast_derived().array() *= time_scaling.array() / dt;
    g_tilde_move_away.const_cast_derived().array() *= time_scaling.array() / dt;

    Eigen::VectorXd dual_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    Eigen::VectorXd primal_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    ADMMContactSolver admm_solver(int(delassus_matrix_plain.rows()));
    admm_solver.setAbsolutePrecision(1e-13);
    admm_solver.setRelativePrecision(1e-14);

    Eigen::VectorXd compliance = Eigen::VectorXd::Zero(g_tilde_move_away.size());
    G_expression.updateCompliance(compliance);
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> preconditioner_vec(
      Eigen::VectorXd::Ones(g_tilde_move_away.size()));
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution_constref(primal_solution);
    const bool has_converged = admm_solver.solve(
      G_expression, g_tilde_move_away, constraint_models, dt, preconditioner_vec,
      primal_solution_constref);
    primal_solution = admm_solver.getPrimalSolution();
    BOOST_CHECK(has_converged);

    dual_solution = (G_plain * primal_solution).cwiseProduct(time_scaling) + g_move_away;
    Eigen::VectorXd dual_solution2 = admm_solver.getDualSolution();

    BOOST_CHECK(std::fabs(primal_solution.dot(dual_solution)) <= 1e-8);
    BOOST_CHECK(primal_solution.isZero());
    BOOST_CHECK(dual_solution.isApprox(g_move_away));
  }
}

BOOST_AUTO_TEST_CASE(joint_limit_translation)
{
  // We test limits for a joint with nq>1
  Model model;
  model.addJoint(0, JointModelTranslation(), SE3::Identity(), "translation");
  model.lowerPositionLimit = Eigen::VectorXd::Ones(model.nq) * -10000;
  model.lowerPositionLimit[2] = 0;
  model.upperPositionLimit = Eigen::VectorXd::Ones(model.nq) * 10000;

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  const double box_mass = 10;
  const Inertia box_inertia = Inertia::FromBox(box_mass, box_dims[0], box_dims[1], box_dims[2]);

  model.appendBodyToJoint(1, box_inertia);
  Data data(model);

  Eigen::VectorXd q0 = neutral(model);
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  Eigen::VectorXd tau_gravity = Eigen::VectorXd::Zero(model.nv);
  tau_gravity(2) = 9.81 * box_mass;

  const double dt = 1e-3;

  typedef JointLimitConstraintModel ConstraintModel;
  typedef ConstraintModel::ConstraintData ConstraintData;
  std::vector<ConstraintModel> constraint_models;
  std::vector<ConstraintData> constraint_datas;

  ConstraintModel joint_limit_constraint_model(model, ConstraintModel::JointIndexVector(1, 1));
  constraint_models.push_back(joint_limit_constraint_model);

  for (const auto & cm : constraint_models)
    constraint_datas.push_back(cm.createData());

  const Eigen::VectorXd v_free_against_lower_bound =
    v0 + dt * aba(model, data, q0, v0, Eigen::VectorXd::Zero(model.nv), Convention::WORLD);
  const Eigen::VectorXd v_free_move_away =
    v0 + dt * aba(model, data, q0, v0, tau_gravity, Convention::WORLD);

  // Cholesky of the Delassus matrix
  crba(model, data, q0, Convention::WORLD);

  data.q_in = q0;
  auto & cmodel = constraint_models[0];
  auto & cdata = constraint_datas[0];
  cmodel.resize(model, data, cdata);
  cmodel.calc(model, data, cdata);
  ContactCholeskyDecomposition chol(model, constraint_models);
  chol.resize(model, constraint_models);
  chol.compute(model, data, constraint_models, constraint_datas, 1e-10);

  auto G_expression = chol.getDelassusCholeskyExpression();
  const auto G_plain = G_expression.matrix();
  const Eigen::MatrixXd delassus_matrix_plain = G_expression.matrix();

  Eigen::MatrixXd constraint_jacobian(cmodel.activeSize(), model.nv);
  constraint_jacobian.setZero();
  getConstraintsJacobian(model, data, constraint_models, constraint_datas, constraint_jacobian);

  Eigen::VectorXd time_scaling(delassus_matrix_plain.rows());
  internal::getTimeScalingFromAccelerationToConstraints(constraint_models, dt, time_scaling);

  // Gravity pushes the freeflyer against the lower bound
  {
    const Eigen::VectorXd g_against_lower_bound = constraint_jacobian * v_free_against_lower_bound;
    const Eigen::VectorXd g_tilde_against_lower_bound =
      g_against_lower_bound + cdata.constraint_residual / dt;
    g_against_lower_bound.const_cast_derived().array() *= time_scaling.array() / dt;
    g_tilde_against_lower_bound.const_cast_derived().array() *= time_scaling.array() / dt;

    Eigen::VectorXd constraint_velocity = Eigen::VectorXd::Zero(cmodel.activeSize());
    Eigen::VectorXd primal_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    ADMMContactSolver admm_solver(int(delassus_matrix_plain.rows()));
    admm_solver.setAbsolutePrecision(1e-13);
    admm_solver.setRelativePrecision(1e-14);

    Eigen::VectorXd compliance = Eigen::VectorXd::Zero(g_tilde_against_lower_bound.size());
    G_expression.updateCompliance(compliance);
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> preconditioner_vec(
      Eigen::VectorXd::Ones(g_tilde_against_lower_bound.size()));
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution_constref(primal_solution);
    const bool has_converged = admm_solver.solve(
      G_expression, g_tilde_against_lower_bound, constraint_models, dt, preconditioner_vec,
      primal_solution_constref);
    primal_solution = admm_solver.getPrimalSolution();
    // std::cout << " admm_solver.getAbsoluteConvergenceResidual():   " <<
    // admm_solver.getAbsoluteConvergenceResidual() << std::endl; std::cout << "
    // admm_solver.getRelativeConvergenceResidual():   " <<
    // admm_solver.getRelativeConvergenceResidual() << std::endl;
    BOOST_CHECK(has_converged);

    constraint_velocity =
      G_plain * primal_solution.cwiseProduct(time_scaling) + g_against_lower_bound;
    constraint_velocity /= dt;
    Eigen::VectorXd dual_solution = admm_solver.getDualSolution();

    BOOST_CHECK(std::fabs(primal_solution.dot(dual_solution)) <= 1e-8);
    BOOST_CHECK(constraint_velocity.isZero(1e-6));
    BOOST_CHECK(
      (dual_solution
       - ((G_plain * primal_solution).cwiseProduct(time_scaling) + g_tilde_against_lower_bound))
        .isZero(1e-6));

    BOOST_CHECK((-tau_gravity + constraint_jacobian.transpose() * primal_solution).isZero(1e-6));
  }

  // External torques compensate the gravity to push the freeflyer away from the lower bound
  {
    const Eigen::VectorXd g_move_away = constraint_jacobian * v_free_move_away;
    const Eigen::VectorXd g_tilde_move_away = g_move_away + cdata.constraint_residual / dt;
    g_move_away.const_cast_derived().array() *= time_scaling.array() / dt;
    g_tilde_move_away.const_cast_derived().array() *= time_scaling.array() / dt;

    Eigen::VectorXd dual_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    Eigen::VectorXd primal_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    ADMMContactSolver admm_solver(int(delassus_matrix_plain.rows()));
    admm_solver.setAbsolutePrecision(1e-13);
    admm_solver.setRelativePrecision(1e-14);

    Eigen::VectorXd compliance = Eigen::VectorXd::Zero(g_tilde_move_away.size());
    G_expression.updateCompliance(compliance);
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> preconditioner_vec(
      Eigen::VectorXd::Ones(g_tilde_move_away.size()));
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution_constref(primal_solution);
    const bool has_converged = admm_solver.solve(
      G_expression, g_tilde_move_away, constraint_models, dt, preconditioner_vec,
      primal_solution_constref);
    primal_solution = admm_solver.getPrimalSolution();
    BOOST_CHECK(has_converged);

    dual_solution = (G_plain * primal_solution).cwiseProduct(time_scaling) + g_move_away;
    Eigen::VectorXd dual_solution2 = admm_solver.getDualSolution();

    BOOST_CHECK(std::fabs(primal_solution.dot(dual_solution)) <= 1e-8);
    BOOST_CHECK(primal_solution.isZero());
    BOOST_CHECK(dual_solution.isApprox(g_move_away));
  }
}

BOOST_AUTO_TEST_CASE(joint_limit_freeflyer)
{
  // We test limits for a joint with nq>1
  Model model;
  model.addJoint(0, JointModelFreeFlyer(), SE3::Identity(), "freeflyer");
  model.lowerPositionLimit = Eigen::VectorXd::Ones(model.nq) * -10000;
  model.lowerPositionLimit[2] = 0;
  model.upperPositionLimit = Eigen::VectorXd::Ones(model.nq) * 10000;

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  const double box_mass = 10;
  const Inertia box_inertia = Inertia::FromBox(box_mass, box_dims[0], box_dims[1], box_dims[2]);

  model.appendBodyToJoint(1, box_inertia);
  Data data(model);

  Eigen::VectorXd q0 = neutral(model);
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  Eigen::VectorXd tau_gravity = Eigen::VectorXd::Zero(model.nv);
  tau_gravity(2) = 9.81 * box_mass;

  const double dt = 1e-3;

  typedef JointLimitConstraintModel ConstraintModel;
  typedef ConstraintModel::ConstraintData ConstraintData;
  std::vector<ConstraintModel> constraint_models;
  std::vector<ConstraintData> constraint_datas;

  ConstraintModel joint_limit_constraint_model(model, ConstraintModel::JointIndexVector(1, 1));
  constraint_models.push_back(joint_limit_constraint_model);

  for (const auto & cm : constraint_models)
    constraint_datas.push_back(cm.createData());

  const Eigen::VectorXd v_free_against_lower_bound =
    v0 + dt * aba(model, data, q0, v0, Eigen::VectorXd::Zero(model.nv), Convention::WORLD);
  const Eigen::VectorXd v_free_move_away =
    v0 + dt * aba(model, data, q0, v0, tau_gravity, Convention::WORLD);

  // Cholesky of the Delassus matrix
  crba(model, data, q0, Convention::WORLD);

  data.q_in = q0;
  auto & cmodel = constraint_models[0];
  auto & cdata = constraint_datas[0];
  cmodel.resize(model, data, cdata);
  cmodel.calc(model, data, cdata);
  ContactCholeskyDecomposition chol(model, constraint_models);
  chol.resize(model, constraint_models);
  chol.compute(model, data, constraint_models, constraint_datas, 1e-10);

  auto G_expression = chol.getDelassusCholeskyExpression();
  const auto G_plain = G_expression.matrix();
  const Eigen::MatrixXd delassus_matrix_plain = G_expression.matrix();

  Eigen::MatrixXd constraint_jacobian(cmodel.activeSize(), model.nv);
  constraint_jacobian.setZero();
  getConstraintsJacobian(model, data, constraint_models, constraint_datas, constraint_jacobian);

  Eigen::VectorXd time_scaling(delassus_matrix_plain.rows());
  internal::getTimeScalingFromAccelerationToConstraints(constraint_models, dt, time_scaling);

  // Gravity pushes the freeflyer against the lower bound
  {
    const Eigen::VectorXd g_against_lower_bound = constraint_jacobian * v_free_against_lower_bound;
    const Eigen::VectorXd g_tilde_against_lower_bound =
      g_against_lower_bound + cdata.constraint_residual / dt;
    g_against_lower_bound.const_cast_derived().array() *= time_scaling.array() / dt;
    g_tilde_against_lower_bound.const_cast_derived().array() *= time_scaling.array() / dt;

    Eigen::VectorXd constraint_velocity = Eigen::VectorXd::Zero(cmodel.activeSize());
    Eigen::VectorXd primal_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    ADMMContactSolver admm_solver(int(delassus_matrix_plain.rows()));
    admm_solver.setAbsolutePrecision(1e-13);
    admm_solver.setRelativePrecision(1e-14);

    Eigen::VectorXd compliance = Eigen::VectorXd::Zero(g_tilde_against_lower_bound.size());
    G_expression.updateCompliance(compliance);
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> preconditioner_vec(
      Eigen::VectorXd::Ones(g_tilde_against_lower_bound.size()));
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution_constref(primal_solution);
    const bool has_converged = admm_solver.solve(
      G_expression, g_tilde_against_lower_bound, constraint_models, dt, preconditioner_vec,
      primal_solution_constref);
    primal_solution = admm_solver.getPrimalSolution();
    BOOST_CHECK(has_converged);

    constraint_velocity =
      G_plain * primal_solution.cwiseProduct(time_scaling) + g_against_lower_bound;
    Eigen::VectorXd dual_solution = admm_solver.getDualSolution();

    BOOST_CHECK(std::fabs(primal_solution.dot(dual_solution)) <= 1e-8);
    BOOST_CHECK(constraint_velocity.isZero(1e-6));
    BOOST_CHECK(
      (dual_solution
       - (G_plain * primal_solution.cwiseProduct(time_scaling) + g_tilde_against_lower_bound))
        .isZero(1e-6));

    BOOST_CHECK((-tau_gravity + constraint_jacobian.transpose() * primal_solution).isZero(1e-6));
  }

  // External torques compensate the gravity to push the freeflyer away from the lower bound
  {
    const Eigen::VectorXd g_move_away = constraint_jacobian * v_free_move_away;
    const Eigen::VectorXd g_tilde_move_away = g_move_away + cdata.constraint_residual / dt;
    g_move_away.const_cast_derived().array() *= time_scaling.array() / dt;
    g_tilde_move_away.const_cast_derived().array() *= time_scaling.array() / dt;

    Eigen::VectorXd dual_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    Eigen::VectorXd primal_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    ADMMContactSolver admm_solver(int(delassus_matrix_plain.rows()));
    admm_solver.setAbsolutePrecision(1e-13);
    admm_solver.setRelativePrecision(1e-14);

    Eigen::VectorXd compliance = Eigen::VectorXd::Zero(g_tilde_move_away.size());
    G_expression.updateCompliance(compliance);
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> preconditioner_vec(
      Eigen::VectorXd::Ones(g_tilde_move_away.size()));
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution_constref(primal_solution);
    const bool has_converged = admm_solver.solve(
      G_expression, g_tilde_move_away, constraint_models, dt, preconditioner_vec,
      primal_solution_constref);
    primal_solution = admm_solver.getPrimalSolution();
    BOOST_CHECK(has_converged);

    dual_solution = (G_plain * primal_solution).cwiseProduct(time_scaling) + g_move_away;
    Eigen::VectorXd dual_solution2 = admm_solver.getDualSolution();

    BOOST_CHECK(std::fabs(primal_solution.dot(dual_solution)) <= 1e-8);
    BOOST_CHECK(primal_solution.isZero());
    BOOST_CHECK(dual_solution.isApprox(g_move_away));
  }
}

BOOST_AUTO_TEST_CASE(joint_limit_composite)
{
  // We test limits for a joint with nq>1
  JointModelComposite joint;
  joint.addJoint(JointModelRX());
  joint.addJoint(JointModelRY());
  Model model;
  model.addJoint(0, joint, SE3::Identity(), "composite");
  model.lowerPositionLimit = Eigen::VectorXd::Ones(model.nq) * -10000;
  model.lowerPositionLimit[1] = 0;
  model.upperPositionLimit = Eigen::VectorXd::Ones(model.nq) * 10000;

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  const double box_mass = 10;
  const Inertia box_inertia = Inertia::FromBox(box_mass, box_dims[0], box_dims[1], box_dims[2]);

  model.appendBodyToJoint(1, box_inertia);
  model.gravity.setZero();
  Data data(model);

  Eigen::VectorXd q0 = Eigen::VectorXd::Zero(model.nq);
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd tau_push_against_lower_bound = -Eigen::VectorXd::Ones(model.nv);

  const double dt = 1e-3;

  typedef JointLimitConstraintModel ConstraintModel;
  typedef ConstraintModel::ConstraintData ConstraintData;
  std::vector<ConstraintModel> constraint_models;
  std::vector<ConstraintData> constraint_datas;

  ConstraintModel joint_limit_constraint_model(model, ConstraintModel::JointIndexVector(1, 1));
  constraint_models.push_back(joint_limit_constraint_model);

  for (const auto & cm : constraint_models)
    constraint_datas.push_back(cm.createData());

  const Eigen::VectorXd v_free_against_lower_bound =
    v0 + dt * aba(model, data, q0, v0, tau_push_against_lower_bound, Convention::WORLD);
  const Eigen::VectorXd v_free_move_away =
    v0 + dt * aba(model, data, q0, v0, -tau_push_against_lower_bound, Convention::WORLD);

  // Cholesky of the Delassus matrix
  crba(model, data, q0, Convention::WORLD);

  data.q_in = q0;
  auto & cmodel = constraint_models[0];
  auto & cdata = constraint_datas[0];
  cmodel.resize(model, data, cdata);
  cmodel.calc(model, data, cdata);
  ContactCholeskyDecomposition chol(model, constraint_models);
  chol.resize(model, constraint_models);
  chol.compute(model, data, constraint_models, constraint_datas, 1e-10);

  auto G_expression = chol.getDelassusCholeskyExpression();
  const auto G_plain = G_expression.matrix();
  const Eigen::MatrixXd delassus_matrix_plain = G_expression.matrix();

  Eigen::MatrixXd constraint_jacobian(cmodel.activeSize(), model.nv);
  constraint_jacobian.setZero();
  getConstraintsJacobian(model, data, constraint_models, constraint_datas, constraint_jacobian);

  Eigen::VectorXd time_scaling(delassus_matrix_plain.rows());
  internal::getTimeScalingFromAccelerationToConstraints(constraint_models, dt, time_scaling);

  // External torques push the freeflyer against from the lower bound
  {
    const Eigen::VectorXd g_against_lower_bound = constraint_jacobian * v_free_against_lower_bound;
    const Eigen::VectorXd g_tilde_against_lower_bound =
      g_against_lower_bound + cdata.constraint_residual / dt;
    g_against_lower_bound.const_cast_derived().array() *= time_scaling.array() / dt;
    g_tilde_against_lower_bound.const_cast_derived().array() *= time_scaling.array() / dt;

    Eigen::VectorXd constraint_velocity = Eigen::VectorXd::Zero(cmodel.activeSize());
    Eigen::VectorXd primal_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    ADMMContactSolver admm_solver(int(delassus_matrix_plain.rows()));
    admm_solver.setAbsolutePrecision(1e-13);
    admm_solver.setRelativePrecision(1e-14);

    Eigen::VectorXd compliance = Eigen::VectorXd::Zero(g_tilde_against_lower_bound.size());
    G_expression.updateCompliance(compliance);
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> preconditioner_vec(
      Eigen::VectorXd::Ones(g_tilde_against_lower_bound.size()));
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution_constref(primal_solution);
    const bool has_converged = admm_solver.solve(
      G_expression, g_tilde_against_lower_bound, constraint_models, dt, preconditioner_vec,
      primal_solution_constref);
    primal_solution = admm_solver.getPrimalSolution();
    BOOST_CHECK(has_converged);

    constraint_velocity =
      G_plain * primal_solution.cwiseProduct(time_scaling) + g_against_lower_bound;

    Eigen::VectorXd dual_solution = admm_solver.getDualSolution();

    BOOST_CHECK(std::fabs(primal_solution.dot(dual_solution)) <= 1e-8);
    BOOST_CHECK(std::abs(constraint_velocity[0]) < 1e-6);
    BOOST_CHECK(
      (dual_solution
       - (G_plain * primal_solution.cwiseProduct(time_scaling) + g_tilde_against_lower_bound))
        .isZero(1e-6));

    BOOST_CHECK(
      std::abs(
        (tau_push_against_lower_bound + constraint_jacobian.transpose() * primal_solution)(1))
      < 1e-6);
  }

  // External torques push the freeflyer away from the lower bound
  {
    const Eigen::VectorXd g_move_away = constraint_jacobian * v_free_move_away;
    const Eigen::VectorXd g_tilde_move_away = g_move_away + cdata.constraint_residual / dt;
    g_move_away.const_cast_derived().array() *= time_scaling.array() / dt;
    g_tilde_move_away.const_cast_derived().array() *= time_scaling.array() / dt;

    Eigen::VectorXd dual_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    Eigen::VectorXd primal_solution = Eigen::VectorXd::Zero(cmodel.activeSize());
    ADMMContactSolver admm_solver(int(delassus_matrix_plain.rows()));
    admm_solver.setAbsolutePrecision(1e-13);
    admm_solver.setRelativePrecision(1e-14);

    Eigen::VectorXd compliance = Eigen::VectorXd::Zero(g_tilde_move_away.size());
    G_expression.updateCompliance(compliance);
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> preconditioner_vec(
      Eigen::VectorXd::Ones(g_tilde_move_away.size()));
    boost::optional<Eigen::Ref<const Eigen::VectorXd>> primal_solution_constref(primal_solution);
    const bool has_converged = admm_solver.solve(
      G_expression, g_tilde_move_away, constraint_models, dt, preconditioner_vec,
      primal_solution_constref);
    primal_solution = admm_solver.getPrimalSolution();
    BOOST_CHECK(has_converged);

    dual_solution = (G_plain * primal_solution).cwiseProduct(time_scaling) + g_move_away;
    Eigen::VectorXd dual_solution2 = admm_solver.getDualSolution();

    BOOST_CHECK(std::fabs(primal_solution.dot(dual_solution)) <= 1e-8);
    BOOST_CHECK(primal_solution.isZero());
    BOOST_CHECK(dual_solution.isApprox(g_move_away));
  }
}

BOOST_AUTO_TEST_SUITE_END()
