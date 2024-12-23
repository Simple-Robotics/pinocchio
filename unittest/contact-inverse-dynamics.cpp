//
// Copyright (c) 2024 INRIA
//

#include <iostream>

#include "pinocchio/algorithm/contact-inverse-dynamics.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include "pinocchio/algorithm/constraints/constraints.hpp"
#include "pinocchio/multibody/sample-models.hpp"
#include "pinocchio/algorithm/contact-solver-utils.hpp"

#include <boost/test/unit_test.hpp>

using namespace Eigen;
using namespace pinocchio;
using namespace pinocchio::internal;

typedef PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(FrictionalPointConstraintModel)
  FrictionalPointConstraintModelVector;
typedef PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(FrictionalPointConstraintData)
  FrictionalPointConstraintDataVector;

void init(Model & model, FrictionalPointConstraintModelVector & constraint_models)
{
  pinocchio::buildModels::humanoidRandom(model, true);
  model.lowerPositionLimit.head<3>().fill(-1.);
  model.upperPositionLimit.head<3>().fill(1.);

  const std::string RF = "rleg6_joint";
  FrictionalPointConstraintModel ci_RF(model, model.getJointId(RF));
  ci_RF.set().mu = 0.4;
  constraint_models.push_back(ci_RF);

  const std::string LF = "lleg6_joint";
  FrictionalPointConstraintModel ci_LF(model, model.getJointId(LF));
  ci_LF.set().mu = 0.4;
  constraint_models.push_back(ci_LF);
}

FrictionalPointConstraintDataVector
createData(const FrictionalPointConstraintModelVector & constraint_models)
{
  FrictionalPointConstraintDataVector constraint_datas;

  for (const auto & cmodel : constraint_models)
  {
    constraint_datas.push_back(cmodel.createData());
  }

  return constraint_datas;
}

template<typename VectorLike>
typename PINOCCHIO_EIGEN_PLAIN_TYPE(VectorLike) abs(const Eigen::MatrixBase<VectorLike> & vec)
{
  return Eigen::abs(vec.array()).matrix();
}

template<typename VectorLike>
void makeIsotropic(
  FrictionalPointConstraintModelVector & constraint_models,
  const Eigen::MatrixBase<VectorLike> & vec_)
{
  auto & vec = vec_.const_cast_derived();

  Eigen::Index row_id = 0;
  for (const auto & cmodel : constraint_models)
  {
    const auto csize = cmodel.size();

    auto vec_seg = vec.segment(row_id, csize);
    vec_seg[1] = vec_seg[0];

    row_id += csize;
  }
}

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(test_contact_inverse_dynamics_3D)
{
#ifdef NDEBUG
  const int num_tests = int(1e4);
#else
  const int num_tests = int(1e2);
#endif

  Model model;
  FrictionalPointConstraintModelVector constraint_models;

  init(model, constraint_models);
  const double mu_prox = 1e-4;
  VectorXd q = randomConfiguration(model);
  VectorXd v = VectorXd::Random(model.nv);
  VectorXd tau = VectorXd::Random(model.nv);
  VectorXd a = Eigen::VectorXd::Zero(model.nv);

  Eigen::DenseIndex constraint_dim = 0;
  for (const auto & cmodel : constraint_models)
  {
    constraint_dim += cmodel.size();
  }

  BOOST_CHECK(constraint_dim > 0);

  for (int n = 0; n < num_tests; ++n)
  {
    Eigen::VectorXd R = abs(Eigen::VectorXd::Random(constraint_dim))
                        + Eigen::VectorXd::Constant(constraint_dim, 1e-10);
    makeIsotropic(constraint_models, R);

    const Eigen::VectorXd R_inv = Eigen::inverse(R.array());
    BOOST_CHECK(MatrixXd(R_inv.asDiagonal() * R.asDiagonal()).isIdentity());

    ProximalSettings prox_settings(1e-12, 1e-12, /*mu = */ 0, 200);

    const Eigen::VectorXd x_positive = abs(Eigen::VectorXd::Random(constraint_dim));
    const Eigen::VectorXd x_in_cone = Eigen::VectorXd::Zero(constraint_dim);

    computeConeProjection(constraint_models, x_positive, x_in_cone);

    const Eigen::VectorXd constraint_velocity_ref = -(R.asDiagonal() * x_in_cone).eval();
    const Eigen::VectorXd sigma_ref = (constraint_velocity_ref + R.asDiagonal() * x_in_cone);
    BOOST_CHECK(sigma_ref.isZero());

    Eigen::VectorXd x_sol = Eigen::VectorXd::Zero(constraint_dim);

    bool has_converged = computeInverseDynamicsConstraintForces(
      constraint_models, constraint_velocity_ref, R, x_sol, prox_settings);
    BOOST_CHECK(has_converged);

    Eigen::VectorXd sigma = constraint_velocity_ref + R.asDiagonal() * x_sol;

    Eigen::VectorXd sigma_correction(sigma);
    computeDeSaxeCorrection(constraint_models, sigma, sigma_correction);
    sigma += sigma_correction;

    BOOST_CHECK(sigma.isZero(1e-10));
    Eigen::VectorXd sigma_projected(sigma);
    computeDualConeProjection(constraint_models, sigma, sigma_projected);
    BOOST_CHECK((sigma_projected - sigma).lpNorm<Eigen::Infinity>() <= 1e-10);
  }
}

BOOST_AUTO_TEST_CASE(test_contact_whole_body_inverse_dynamics_3D)
{
}

BOOST_AUTO_TEST_SUITE_END()
