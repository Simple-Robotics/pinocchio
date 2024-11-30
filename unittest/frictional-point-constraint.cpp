//
// Copyright (c) 2024 INRIA
//

#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/contact-info.hpp"
#include "pinocchio/algorithm/contact-cholesky.hpp"
#include "pinocchio/algorithm/contact-jacobian.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/multibody/sample-models.hpp"
#include "pinocchio/utils/timer.hpp"
#include "pinocchio/spatial/classic-acceleration.hpp"
#include "pinocchio/algorithm/constraints/bilateral-point-constraint.hpp"

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;
using namespace Eigen;

template<typename T>
bool within(const T & elt, const std::vector<T> & vec)
{
  typename std::vector<T>::const_iterator it;

  it = std::find(vec.begin(), vec.end(), elt);
  if (it != vec.end())
    return true;
  else
    return false;
}

template<typename Matrix>
bool within(const typename Matrix::Scalar & elt, const Eigen::MatrixBase<Matrix> & mat)
{
  for (DenseIndex i = 0; i < mat.rows(); ++i)
    for (DenseIndex j = 0; j < mat.rows(); ++j)
    {
      if (elt == mat(i, j))
        return true;
    }

  return false;
}

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(basic_constructor)
{
  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);

  // Check complete constructor
  const SE3 M(SE3::Random());
  BilateralPointConstraintModel cmodel2(model, 0, M);
  BOOST_CHECK(cmodel2.joint1_id == 0);
  BOOST_CHECK(cmodel2.joint1_placement == M);
  BOOST_CHECK(cmodel2.size() == 3);

  // Check contructor with two arguments
  BilateralPointConstraintModel cmodel2prime(model, 0);
  BOOST_CHECK(cmodel2prime.joint1_id == 0);
  BOOST_CHECK(cmodel2prime.joint1_placement.isIdentity(0.));
  BOOST_CHECK(cmodel2prime.size() == 3);

  // Check default copy constructor
  BilateralPointConstraintModel cmodel3(cmodel2);
  BOOST_CHECK(cmodel3 == cmodel2);
}

void check_A1_and_A2(
  const Model & model,
  const Data & data,
  const BilateralPointConstraintModel & cmodel,
  BilateralPointConstraintData & cdata)
{
  const BilateralPointConstraintModel::Matrix36 A1_world = cmodel.getA1(cdata, WorldFrame());
  BilateralPointConstraintModel::Matrix36 A1_world_ref =
    -cdata.oMc1.toActionMatrixInverse().topRows<3>();
  A1_world_ref.rightCols<3>() +=
    skew(cdata.constraint_position_error) * cdata.oMc1.rotation().transpose();

  BOOST_CHECK(A1_world.isApprox(A1_world_ref));

  const BilateralPointConstraintModel::Matrix36 A2_world = cmodel.getA2(cdata, WorldFrame());
  const BilateralPointConstraintModel::Matrix36 A2_world_ref =
    cdata.c1Mc2.rotation() * cdata.oMc2.toActionMatrixInverse().topRows<3>();

  BOOST_CHECK(A2_world.isApprox(A2_world_ref));

  const BilateralPointConstraintModel::Matrix36 A1_local = cmodel.getA1(cdata, LocalFrame());
  BilateralPointConstraintModel::Matrix36 A1_local_ref =
    -cmodel.joint1_placement.toActionMatrixInverse().topRows<3>();
  A1_local_ref.rightCols<3>() +=
    skew(cdata.constraint_position_error) * cmodel.joint1_placement.rotation().transpose();

  BOOST_CHECK(A1_local.isApprox(A1_local_ref));

  const RigidConstraintModel::Matrix36 A2_local = cmodel.getA2(cdata, LocalFrame());
  const RigidConstraintModel::Matrix36 A2_local_ref =
    cdata.c1Mc2.rotation() * cmodel.joint2_placement.toActionMatrixInverse().topRows<3>();

  BOOST_CHECK(A2_local.isApprox(A2_local_ref));

  // Check Jacobians
  Data::MatrixXs J_ref(3, model.nv);
  J_ref.setZero();
  getConstraintJacobian(model, data, cmodel, cdata, J_ref);

  // World
  const Data::Matrix6x J1_world = getJointJacobian(model, data, cmodel.joint1_id, WORLD);
  const Data::Matrix6x J2_world = getJointJacobian(model, data, cmodel.joint2_id, WORLD);
  const Data::Matrix3x J_world = A1_world * J1_world + A2_world * J2_world;

  BOOST_CHECK(J_world.isApprox(J_ref));

  // Local
  const Data::Matrix6x J1_local = getJointJacobian(model, data, cmodel.joint1_id, LOCAL);
  const Data::Matrix6x J2_local = getJointJacobian(model, data, cmodel.joint2_id, LOCAL);
  const Data::Matrix3x J_local = A1_local * J1_local + A2_local * J2_local;

  BOOST_CHECK(J_local.isApprox(J_ref));

  // Check Jacobian matrix product
  const Eigen::DenseIndex m = 40;
  const Data::MatrixXs mat = Data::MatrixXs::Random(model.nv, m);

  Data::MatrixXs res(cmodel.size(), m);
  res.setZero();
  cmodel.jacobianMatrixProduct(model, data, cdata, mat, res);

  const Data::MatrixXs res_ref = J_ref * mat;

  BOOST_CHECK(res.isApprox(res_ref));
}

BOOST_AUTO_TEST_CASE(constraint3D_basic_operations)
{
  const pinocchio::Model model;
  const pinocchio::Data data(model);
  RigidConstraintModel cm(CONTACT_3D, model, 0, SE3::Random(), LOCAL);
  RigidConstraintData cd(cm);
  cm.calc(model, data, cd);

  const pinocchio::SE3 placement = cm.joint1_placement;

  {
    const Eigen::Vector3d diagonal_inertia(1, 2, 3);

    const pinocchio::SE3::Matrix6 spatial_inertia =
      cm.computeConstraintSpatialInertia(placement, diagonal_inertia);
    BOOST_CHECK(spatial_inertia.transpose().isApprox(spatial_inertia)); // check symmetric matrix

    const auto A1 = cm.getA1(cd, LocalFrame());
    const pinocchio::SE3::Matrix6 spatial_inertia_ref =
      A1.transpose() * diagonal_inertia.asDiagonal() * A1;

    BOOST_CHECK(spatial_inertia.isApprox(spatial_inertia_ref));
  }

  // Scalar
  {
    const double constant_value = 10;
    const Eigen::Vector3d diagonal_inertia = Eigen::Vector3d::Constant(constant_value);

    const pinocchio::SE3::Matrix6 spatial_inertia =
      cm.computeConstraintSpatialInertia(placement, diagonal_inertia);
    BOOST_CHECK(spatial_inertia.transpose().isApprox(spatial_inertia)); // check symmetric matrix

    const auto A1 = cm.getA1(cd, LocalFrame());
    const pinocchio::SE3::Matrix6 spatial_inertia_ref =
      A1.transpose() * diagonal_inertia.asDiagonal() * A1;

    BOOST_CHECK(spatial_inertia.isApprox(spatial_inertia_ref));

    const Inertia spatial_inertia_ref2(constant_value, placement.translation(), Symmetric3::Zero());
    BOOST_CHECK(spatial_inertia.isApprox(spatial_inertia_ref2.matrix()));
  }
}

template<typename VectorLike>
Eigen::MatrixXd compute_jacobian_fd(
  const Model & model,
  const BilateralPointConstraintModel & cmodel,
  const Eigen::MatrixBase<VectorLike> & q,
  const double eps)
{
  Data data_fd(model), data(model);
  BilateralPointConstraintData cdata(cmodel), cdata_fd(cmodel);

  Eigen::MatrixXd res(3, model.nv);
  res.setZero();

  forwardKinematics(model, data, q),

    cmodel.calc(model, data, cdata);
  Eigen::VectorXd v_plus(model.nv);
  v_plus.setZero();

  for (int i = 0; i < model.nv; ++i)
  {
    v_plus[i] = eps;
    const auto q_plus = integrate(model, q, v_plus);
    forwardKinematics(model, data_fd, q_plus),

      cmodel.calc(model, data_fd, cdata_fd);

    res.col(i) = (cdata_fd.constraint_position_error - cdata.constraint_position_error) / eps;

    v_plus[i] = 0;
  }

  return res;
}

Vector3d computeConstraintError(
  const Model & model, const Data & data, const BilateralPointConstraintModel & cm)
{
  PINOCCHIO_UNUSED_VARIABLE(model);

  const SE3 oMc1 = data.oMi[cm.joint1_id] * cm.joint1_placement;
  const SE3 oMc2 = data.oMi[cm.joint2_id] * cm.joint2_placement;

  const Vector3d error_world_frame = oMc2.translation() - oMc1.translation();
  const Vector3d error_local_frame1 = oMc1.rotation().transpose() * error_world_frame;

  return error_local_frame1;
}

BOOST_AUTO_TEST_CASE(contact_models_sparsity_and_jacobians)
{

  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);
  Data data(model);

  model.lowerPositionLimit.head<3>().fill(-1.);
  model.upperPositionLimit.head<3>().fill(1.);
  VectorXd q = randomConfiguration(model);
  VectorXd v = VectorXd::Random(model.nv);
  VectorXd a = VectorXd::Random(model.nv);

  forwardKinematics(model, data, q, v);
  computeJointJacobians(model, data, q);

  const std::string RF = "rleg6_joint";
  const std::string LF = "lleg6_joint";
  const double eps_fd = 1e-8;

  const BilateralPointConstraintModel cm_RF(model, model.getJointId(RF), SE3::Random());
  const BilateralPointConstraintModel cm_LF(model, model.getJointId(LF), SE3::Random());
  const BilateralPointConstraintModel clm_RF_LF(
    model, cm_RF.joint1_id, cm_RF.joint1_placement, cm_LF.joint1_id, cm_LF.joint1_placement);

  // Check errors values
  {
    Data data(model);

    BilateralPointConstraintData cd_RF(cm_RF);
    BilateralPointConstraintData cd_LF(cm_LF);
    BilateralPointConstraintData cld_RF_LF(clm_RF_LF);

    forwardKinematics(model, data, q);

    cm_RF.calc(model, data, cd_RF);
    const Vector3d position_error_RF_ref = computeConstraintError(model, data, cm_RF);
    BOOST_CHECK(cd_RF.constraint_position_error.isApprox(position_error_RF_ref));

    cm_LF.calc(model, data, cd_LF);
    const Vector3d position_error_LF_ref = computeConstraintError(model, data, cm_LF);
    BOOST_CHECK(cd_LF.constraint_position_error.isApprox(position_error_LF_ref));

    clm_RF_LF.calc(model, data, cld_RF_LF);
    const Vector3d position_error_RF_LF_ref = computeConstraintError(model, data, clm_RF_LF);
    BOOST_CHECK(cld_RF_LF.constraint_position_error.isApprox(position_error_RF_LF_ref));
  }
  {
    forwardKinematics(model, data, q, v, a);

    BilateralPointConstraintData cd_RF(cm_RF);
    cm_RF.calc(model, data, cd_RF);
    BilateralPointConstraintData cd_LF(cm_LF);
    cm_LF.calc(model, data, cd_LF);
    BilateralPointConstraintData cld_RF_LF(clm_RF_LF);
    clm_RF_LF.calc(model, data, cld_RF_LF);

    Data::Matrix6x J6_RF_LOCAL(6, model.nv);
    J6_RF_LOCAL.setZero();
    getFrameJacobian(model, data, cm_RF.joint1_id, cm_RF.joint1_placement, LOCAL, J6_RF_LOCAL);
    Data::Matrix3x J_RF_LOCAL(3, model.nv);
    J_RF_LOCAL = -J6_RF_LOCAL.middleRows<3>(SE3::LINEAR);
    J_RF_LOCAL += cross(cd_RF.constraint_position_error, J6_RF_LOCAL.middleRows<3>(SE3::ANGULAR));

    Data::Matrix6x J6_LF_LOCAL(6, model.nv);
    J6_LF_LOCAL.setZero();
    getFrameJacobian(model, data, cm_LF.joint1_id, cm_LF.joint1_placement, LOCAL, J6_LF_LOCAL);
    Data::Matrix3x J_LF_LOCAL(3, model.nv);
    J_LF_LOCAL = -J6_LF_LOCAL.middleRows<3>(SE3::LINEAR);
    J_LF_LOCAL += cross(cd_LF.constraint_position_error, J6_LF_LOCAL.middleRows<3>(SE3::ANGULAR));

    for (DenseIndex k = 0; k < model.nv; ++k)
    {
      BOOST_CHECK(
        J_RF_LOCAL.middleRows<3>(SE3::LINEAR).col(k).isZero() != cm_RF.colwise_joint1_sparsity[k]);
      BOOST_CHECK(
        J_LF_LOCAL.middleRows<3>(SE3::LINEAR).col(k).isZero() != cm_LF.colwise_joint1_sparsity[k]);
    }
    BOOST_CHECK(cm_RF.colwise_joint2_sparsity.isZero());
    BOOST_CHECK(cm_LF.colwise_joint2_sparsity.isZero());

    const SE3 oMc1 = data.oMi[clm_RF_LF.joint1_id] * clm_RF_LF.joint1_placement;
    const SE3 oMc2 = data.oMi[clm_RF_LF.joint2_id] * clm_RF_LF.joint2_placement;
    const SE3 c1Mc2 = oMc1.actInv(oMc2);
    Data::Matrix3x J_clm_LOCAL = c1Mc2.rotation() * J6_LF_LOCAL.middleRows<3>(SE3::LINEAR)
                                 - J6_RF_LOCAL.middleRows<3>(SE3::LINEAR);
    J_clm_LOCAL +=
      cross(cld_RF_LF.constraint_position_error, J6_RF_LOCAL.middleRows<3>(SE3::ANGULAR));

    for (DenseIndex k = 0; k < model.nv; ++k)
    {
      BOOST_CHECK(J_clm_LOCAL.col(k).isZero(0) != within(k, clm_RF_LF.colwise_span_indexes));
    }

    // Check Jacobian vs sparse Jacobian computation
    Data::MatrixXs J_RF_sparse(3, model.nv);
    J_RF_sparse.setZero(); // TODO: change input type when all the API would be refactorized
                           // with CRTP on contact constraints
    getConstraintJacobian(model, data, cm_RF, cd_RF, J_RF_sparse);
    BOOST_CHECK(J_RF_LOCAL.isApprox(J_RF_sparse));

    const auto J_RF_fd = compute_jacobian_fd(model, cm_RF, q, eps_fd);
    BOOST_CHECK(J_RF_sparse.isApprox(J_RF_fd, sqrt(eps_fd)));

    Data::MatrixXs J_LF_sparse(3, model.nv);
    J_LF_sparse.setZero(); // TODO: change input type when all the API would be refactorized
                           // with CRTP on contact constraints
    getConstraintJacobian(model, data, cm_LF, cd_LF, J_LF_sparse);
    BOOST_CHECK(J_LF_LOCAL.middleRows<3>(SE3::LINEAR).isApprox(J_LF_sparse));

    const auto J_LF_fd = compute_jacobian_fd(model, cm_LF, q, eps_fd);
    BOOST_CHECK(J_LF_sparse.isApprox(J_LF_fd, sqrt(eps_fd)));

    Data::MatrixXs J_clm_sparse(3, model.nv);
    J_clm_sparse.setZero(); // TODO: change input type when all the API would be refactorized
                            // with CRTP on contact constraints
    getConstraintJacobian(model, data, clm_RF_LF, cld_RF_LF, J_clm_sparse);
    BOOST_CHECK(J_clm_LOCAL.isApprox(J_clm_sparse));

    const auto J_clm_fd = compute_jacobian_fd(model, clm_RF_LF, q, eps_fd);
    BOOST_CHECK(J_clm_sparse.isApprox(J_clm_fd, sqrt(eps_fd)));

    // Check velocity and acceleration
    {
      const double dt = eps_fd;
      BilateralPointConstraintData cd_RF(cm_RF), cd_RF_plus(cm_RF);
      cm_RF.calc(model, data, cd_RF);

      Data data_plus(model);
      const VectorXd v_plus = v + a * dt;
      const VectorXd q_plus = integrate(model, q, v_plus * dt);
      forwardKinematics(model, data_plus, q_plus, v_plus);

      {
        BilateralPointConstraintData cd_RF(cm_RF), cd_RF_plus(cm_RF);
        cm_RF.calc(model, data, cd_RF);
        BOOST_CHECK(cd_RF.constraint_velocity_error.isApprox(J_RF_sparse * v));

        cm_RF.calc(model, data_plus, cd_RF_plus);
        const Vector3d constraint_velocity_error_fd =
          (cd_RF_plus.constraint_position_error - cd_RF.constraint_position_error) / dt;
        BOOST_CHECK(
          cd_RF.constraint_velocity_error.isApprox(constraint_velocity_error_fd, sqrt(dt)));

        const Vector3d constraint_acceleration_error_fd =
          (cd_RF_plus.constraint_velocity_error - cd_RF.constraint_velocity_error) / dt;
        BOOST_CHECK(
          cd_RF.constraint_acceleration_error.isApprox(constraint_acceleration_error_fd, sqrt(dt)));
      }

      {
        BilateralPointConstraintData cd_LF(cm_LF), cd_LF_plus(cm_LF);
        cm_LF.calc(model, data, cd_LF);
        BOOST_CHECK(cd_LF.constraint_velocity_error.isApprox(J_LF_sparse * v));

        cm_LF.calc(model, data_plus, cd_LF_plus);
        const Vector3d constraint_velocity_error_fd =
          (cd_LF_plus.constraint_position_error - cd_LF.constraint_position_error) / dt;
        BOOST_CHECK(
          cd_LF.constraint_velocity_error.isApprox(constraint_velocity_error_fd, sqrt(dt)));

        const Vector3d constraint_acceleration_error_fd =
          (cd_LF_plus.constraint_velocity_error - cd_LF.constraint_velocity_error) / dt;
        BOOST_CHECK(
          cd_LF.constraint_acceleration_error.isApprox(constraint_acceleration_error_fd, sqrt(dt)));
      }

      {
        BilateralPointConstraintData cld_RF_LF(clm_RF_LF), cld_RF_LF_plus(clm_RF_LF);
        clm_RF_LF.calc(model, data, cld_RF_LF);
        BOOST_CHECK(cld_RF_LF.constraint_velocity_error.isApprox(J_clm_sparse * v));

        clm_RF_LF.calc(model, data_plus, cld_RF_LF_plus);
        const Vector3d constraint_velocity_error_fd =
          (cld_RF_LF_plus.constraint_position_error - cld_RF_LF.constraint_position_error) / dt;
        BOOST_CHECK(
          cld_RF_LF.constraint_velocity_error.isApprox(constraint_velocity_error_fd, sqrt(dt)));

        const Vector3d constraint_acceleration_error_fd =
          (cld_RF_LF_plus.constraint_velocity_error - cld_RF_LF.constraint_velocity_error) / dt;
        BOOST_CHECK(cld_RF_LF.constraint_acceleration_error.isApprox(
          constraint_acceleration_error_fd, sqrt(dt)));
      }
    }

    check_A1_and_A2(model, data, cm_RF, cd_RF);
    check_A1_and_A2(model, data, cm_LF, cd_LF);
    check_A1_and_A2(model, data, clm_RF_LF, cld_RF_LF);

    // Check acceleration contributions
    {
      Data data(model), data_zero_acc(model);
      forwardKinematics(model, data, q, v, a);
      computeJointJacobians(model, data, q);
      forwardKinematics(model, data_zero_acc, q, v, VectorXd::Zero(model.nv));

      // RF
      BilateralPointConstraintData cd_RF(cm_RF), cd_RF_zero_acc(cm_RF);
      cm_RF.calc(model, data, cd_RF);
      cm_RF.calc(model, data_zero_acc, cd_RF_zero_acc);

      Data::MatrixXs J_RF_sparse(3, model.nv);
      J_RF_sparse.setZero();
      cm_RF.jacobian(model, data, cd_RF, J_RF_sparse);

      BOOST_CHECK((J_RF_sparse * a + cd_RF_zero_acc.constraint_acceleration_error)
                    .isApprox(cd_RF.constraint_acceleration_error));

      // LF
      BilateralPointConstraintData cd_LF(cm_LF), cd_LF_zero_acc(cm_LF);
      cm_LF.calc(model, data, cd_LF);
      cm_LF.calc(model, data_zero_acc, cd_LF_zero_acc);

      Data::MatrixXs J_LF_sparse(3, model.nv);
      J_LF_sparse.setZero();
      cm_LF.jacobian(model, data, cd_LF, J_LF_sparse);

      BOOST_CHECK((J_LF_sparse * a + cd_LF_zero_acc.constraint_acceleration_error)
                    .isApprox(cd_LF.constraint_acceleration_error));

      // Close loop
      BilateralPointConstraintData cld_RF_LF(clm_RF_LF), cld_RF_LF_zero_acc(clm_RF_LF);
      clm_RF_LF.calc(model, data, cld_RF_LF);
      clm_RF_LF.calc(model, data_zero_acc, cld_RF_LF_zero_acc);

      Data::MatrixXs J_clm_sparse(3, model.nv);
      J_clm_sparse.setZero();
      clm_RF_LF.jacobian(model, data, cld_RF_LF, J_clm_sparse);

      BOOST_CHECK((J_clm_sparse * a + cld_RF_LF_zero_acc.constraint_acceleration_error)
                    .isApprox(cld_RF_LF.constraint_acceleration_error));
    }
  }
}

BOOST_AUTO_TEST_CASE(cast)
{

  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);

  const std::string RF = "rleg6_joint";

  const BilateralPointConstraintModel cm_RF(model, model.getJointId(RF), SE3::Random());
  const auto cm_RF_cast_double = cm_RF.cast<double>();
  BOOST_CHECK(cm_RF_cast_double == cm_RF);

  const auto cm_RF_cast_long_double = cm_RF.cast<long double>();
  BOOST_CHECK(cm_RF_cast_long_double.cast<double>() == cm_RF);
}

BOOST_AUTO_TEST_CASE(cholesky)
{

  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);
  Data data(model), data_ref(model);

  model.lowerPositionLimit.head<3>().fill(-1.);
  model.upperPositionLimit.head<3>().fill(1.);
  VectorXd q = randomConfiguration(model);
  VectorXd v = VectorXd::Random(model.nv);
  VectorXd a = VectorXd::Random(model.nv);

  crba(model, data, q, Convention::WORLD);

  const std::string RF = "rleg6_joint";
  const std::string LF = "lleg6_joint";

  const BilateralPointConstraintModel cm_RF(model, model.getJointId(RF), SE3::Random());
  const BilateralPointConstraintModel cm_LF(model, model.getJointId(LF), SE3::Random());
  const BilateralPointConstraintModel clm_RF_LF(
    model, cm_RF.joint1_id, cm_RF.joint1_placement, cm_LF.joint1_id, cm_LF.joint1_placement);

  std::vector<BilateralPointConstraintModel> constraint_models;
  constraint_models.push_back(cm_RF);
  constraint_models.push_back(cm_LF);
  constraint_models.push_back(clm_RF_LF);

  std::vector<BilateralPointConstraintData> constraint_datas, constraint_datas_ref;
  for (const auto & cm : constraint_models)
  {
    constraint_datas.push_back(cm.createData());
    constraint_datas_ref.push_back(cm.createData());
  }

  const double mu = 1e-10;
  ContactCholeskyDecomposition cholesky(model, constraint_models);
  cholesky.compute(model, data, constraint_models, constraint_datas, mu);

  crba(model, data_ref, q, Convention::WORLD);
  make_symmetric(data_ref.M);
  const auto total_size = getTotalConstraintSize(constraint_models);
  Eigen::MatrixXd J_constraints(total_size, model.nv);
  J_constraints.setZero();
  getConstraintsJacobian(model, data_ref, constraint_models, constraint_datas, J_constraints);

  Eigen::MatrixXd H_ref = Eigen::MatrixXd::Zero(total_size + model.nv, total_size + model.nv);
  H_ref.topLeftCorner(total_size, total_size).diagonal().fill(-mu);
  H_ref.bottomRightCorner(model.nv, model.nv) = data_ref.M;
  H_ref.topRightCorner(total_size, model.nv) = J_constraints;
  H_ref.bottomLeftCorner(model.nv, total_size) = J_constraints.transpose();

  BOOST_CHECK(cholesky.matrix().isApprox(H_ref));
}

BOOST_AUTO_TEST_SUITE_END()
