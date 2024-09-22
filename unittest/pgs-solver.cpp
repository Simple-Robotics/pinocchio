//
// Copyright (c) 2024 INRIA
//

#include "pinocchio/algorithm/contact-info.hpp"
#include "pinocchio/algorithm/contact-cholesky.hxx"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/contact-jacobian.hpp"
#include "pinocchio/algorithm/pgs-solver.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/crba.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;

double mu = 1e-4;

template<typename _ConstraintModel, typename _ConstraintSet>
struct TestBoxTpl
{
  typedef _ConstraintModel ConstraintModel;
  typedef _ConstraintSet ConstraintSet;

  typedef typename ConstraintModel::ConstraintData ConstraintData;

  TestBoxTpl(
    const Model & model,
    const std::vector<ConstraintModel> & constraint_models,
    const std::vector<ConstraintSet> & constraint_sets)
  : model(model)
  , data(model)
  , constraint_models(constraint_models)
  , constraint_sets(constraint_sets)
  , v_next(Eigen::VectorXd::Zero(model.nv))
  {
    for (const auto & cm : constraint_models)
    {
      constraint_datas.push_back(cm.createData());
    }

    const Eigen::DenseIndex constraint_size = getTotalConstraintSize(constraint_models);
    primal_solution = dual_solution = Eigen::VectorXd::Zero(constraint_size);
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
    //    std::cout << "G:\n" << delassus_matrix_plain << std::endl;

    Eigen::MatrixXd constraint_jacobian(12, model.nv);
    constraint_jacobian.setZero();
    getConstraintsJacobian(model, data, constraint_models, constraint_datas, constraint_jacobian);

    const Eigen::VectorXd g = constraint_jacobian * v_free;
    //    std::cout << "g: " << g.transpose() << std::endl;

    PGSContactSolver pgs_solver(int(delassus_matrix_plain.rows()));
    pgs_solver.setAbsolutePrecision(1e-10);
    pgs_solver.setRelativePrecision(1e-14);
    has_converged = pgs_solver.solve(G, g, constraint_sets, dual_solution);

    //    std::cout << "x_sol: " << x_sol.transpose() << std::endl;

    primal_solution = G * dual_solution + g;
    //    std::cout << "constraint_velocity: " << constraint_velocity.transpose() << std::endl;

    f_tot.setZero();
    for (int k = 0; k < 4; ++k)
    {
      f_tot += dual_solution.segment(3 * k, 3);
    }

    //    std::cout << "f_tot: " << f_tot.transpose() << std::endl;
    for (int k = 0; k < 4; ++k)
    {
      const auto & cm = constraint_models[size_t(k)];
      Force contact_force = Force::Zero();
      contact_force.linear() = dual_solution.segment(3 * k, 3) / dt;
      external_forces[cm.joint1_id] += cm.joint1_placement.act(contact_force);
      external_forces[cm.joint2_id] -= cm.joint2_placement.act(contact_force);
    }
    v_next = v0 + dt * aba(model, data, q0, v0, tau0, external_forces, Convention::WORLD);
  }

  Model model;
  Data data;
  std::vector<ConstraintModel> constraint_models;
  std::vector<ConstraintData> constraint_datas;
  std::vector<ConstraintSet> constraint_sets;
  Eigen::VectorXd v_next;

  Eigen::VectorXd primal_solution, dual_solution;
  Eigen::Vector3d f_tot;
  bool has_converged;
};

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

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

  typedef RigidConstraintModel ConstraintModel;
  typedef CoulombFrictionCone ConstraintSet;
  typedef TestBoxTpl<ConstraintModel, ConstraintSet> TestBox;
  std::vector<ConstraintModel> constraint_models;

  {
    const SE3 local_placement1(
      SE3::Matrix3::Identity(), 0.5 * SE3::Vector3(box_dims[0], box_dims[1], -box_dims[2]));
    SE3::Matrix3 rot = SE3::Matrix3::Identity();
    for (int i = 0; i < 4; ++i)
    {
      const SE3 local_placement(SE3::Matrix3::Identity(), rot * local_placement1.translation());
      ConstraintModel cm(CONTACT_3D, model, 1, local_placement);
      constraint_models.push_back(cm);
      rot = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitZ()).toRotationMatrix() * rot;
    }

    for (const auto & cm : constraint_models)
    {
      //      std::cout << "placement cm:\n" << cm.joint1_placement << std::endl;
    }
  }

  // Create stack of friction cones
  const double friction_value = 0.4;
  std::vector<ConstraintSet> constraint_sets;
  for (int k = 0; k < 4; ++k)
  {
    constraint_sets.push_back(CoulombFrictionCone(friction_value));
  }

  // Test static motion with zero external force
  {
    const Force fext = Force::Zero();

    TestBox test(model, constraint_models, constraint_sets);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.primal_solution.isZero(2e-10));
    BOOST_CHECK(test.f_tot.isApprox(-box_mass * Model::gravity981 * dt, 1e-8));
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

    TestBox test(model, constraint_models, constraint_sets);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.primal_solution.isZero(1e-8));
    const Force::Vector3 f_tot_ref = (-box_mass * Model::gravity981 - fext.linear()) * dt;
    BOOST_CHECK(test.f_tot.isApprox(f_tot_ref, 1e-6));
    BOOST_CHECK(test.v_next.isZero(1e-8));
  }

  // Test slidding motion
  for (int k = 0; k < num_tests; ++k)
  {
    const double scaling = 1.1;
    Force fext = Force::Zero();
    fext.linear().head<2>().setRandom().normalize();
    fext.linear() *= scaling * f_sliding;

    TestBox test(model, constraint_models, constraint_sets);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    const Force::Vector3 f_tot_ref =
      (-box_mass * Model::gravity981 - 1 / scaling * fext.linear()) * dt;
    BOOST_CHECK(test.f_tot.isApprox(f_tot_ref, 1e-6));
    BOOST_CHECK(
      math::fabs(Motion(test.v_next).linear().norm() - (f_sliding * 0.1 / box_mass * dt)) <= 1e-6);
    BOOST_CHECK(Motion(test.v_next).angular().isZero(1e-6));
  }
}

BOOST_AUTO_TEST_SUITE_END()
