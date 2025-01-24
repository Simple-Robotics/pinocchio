//
// Copyright (c) 2024-2025 INRIA
//

#include "pinocchio/math/lanczos-decomposition.hpp"
#include "pinocchio/algorithm/delassus-operator-dense.hpp"
#include "pinocchio/algorithm/contact-cholesky.hpp"
#include "pinocchio/algorithm/crba.hpp"

#include <boost/variant.hpp> // to avoid C99 warnings

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

using namespace pinocchio;

int line;
#define SET_LINE line = __LINE__ + 1
#define PINOCCHIO_CHECK(cond) BOOST_CHECK_MESSAGE(cond, "from line " << line << ": " #cond)
#define PINOCCHIO_CHECK_EQUAL(a, b)                                                                \
  BOOST_CHECK_MESSAGE(                                                                             \
    (a) == (b), "from line " << line << ": " #a "[" << (a) << "] != " #b "[" << (b) << "].")

BOOST_AUTO_TEST_CASE(test_basic_constructor)
{
  const Eigen::DenseIndex mat_size = 20;
  const Eigen::DenseIndex decomposition_size = 10;

  typedef LanczosDecompositionTpl<Eigen::MatrixXd> LanczosDecomposition;
  LanczosDecomposition lanczos_decomposition(mat_size, decomposition_size);

  BOOST_CHECK(lanczos_decomposition.size() == mat_size);
  BOOST_CHECK(lanczos_decomposition.decompositionSize() == decomposition_size);
}

BOOST_AUTO_TEST_CASE(test_identity)
{
  const Eigen::DenseIndex mat_size = 20;
  const Eigen::MatrixXd identity_matrix = Eigen::MatrixXd::Identity(mat_size, mat_size);

  typedef LanczosDecompositionTpl<Eigen::MatrixXd> LanczosDecomposition;
  LanczosDecomposition lanczos_decomposition(identity_matrix, 10);

  const auto residual = lanczos_decomposition.computeDecompositionResidual(identity_matrix);
  BOOST_CHECK(residual.isZero());

  BOOST_CHECK(lanczos_decomposition.rank() == 1);
  BOOST_CHECK((lanczos_decomposition.Qs().transpose() * lanczos_decomposition.Qs())
                .topLeftCorner(lanczos_decomposition.rank(), lanczos_decomposition.rank())
                .isIdentity());
}

BOOST_AUTO_TEST_CASE(test_diagonal_matrix)
{
  const Eigen::DenseIndex mat_size = 20;
  const Eigen::VectorXd diagonal_terms = Eigen::VectorXd::LinSpaced(mat_size, 0.0, mat_size - 1);
  const Eigen::MatrixXd diagonal_matrix = diagonal_terms.asDiagonal();

  typedef LanczosDecompositionTpl<Eigen::MatrixXd> LanczosDecomposition;
  LanczosDecomposition lanczos_decomposition(diagonal_matrix, 10);

  const auto residual = lanczos_decomposition.computeDecompositionResidual(diagonal_matrix);
  BOOST_CHECK(residual.isZero());

  for (Eigen::DenseIndex col_id = 0; col_id < lanczos_decomposition.Ts().cols(); ++col_id)
  {
    BOOST_CHECK(math::fabs(lanczos_decomposition.Qs().col(col_id).norm() - 1.) <= 1e-12);
  }

  BOOST_CHECK(lanczos_decomposition.rank() == lanczos_decomposition.Ts().cols());
  BOOST_CHECK((lanczos_decomposition.Qs().transpose() * lanczos_decomposition.Qs()).isIdentity());
}

void checkDecomposition(
  const LanczosDecompositionTpl<Eigen::MatrixXd> & lanczos_decomposition,
  const Eigen::MatrixXd & matrix,
  const bool full_rank = true)
{
  const auto residual = lanczos_decomposition.computeDecompositionResidual(matrix);

  const auto & Ts = lanczos_decomposition.Ts();
  const auto & Qs = lanczos_decomposition.Qs();
  const Eigen::MatrixXd residual_bis = matrix * Qs - Qs * Ts.matrix();
  const Eigen::MatrixXd residual_tierce = Qs.transpose() * matrix * Qs - Ts.matrix();
  const Eigen::MatrixXd residual_fourth = matrix - Qs * Ts.matrix() * Qs.transpose();
  PINOCCHIO_CHECK(residual.isZero(1e-10));

  const auto rank = lanczos_decomposition.rank();
  if (full_rank)
    PINOCCHIO_CHECK_EQUAL(rank, lanczos_decomposition.Ts().cols());

  for (Eigen::DenseIndex col_id = 0; col_id < rank; ++col_id)
  {
    PINOCCHIO_CHECK(math::fabs(lanczos_decomposition.Qs().col(col_id).norm() - 1.) <= 1e-12);
  }

  for (Eigen::DenseIndex col_id = rank; lanczos_decomposition.decompositionSize() < rank; ++col_id)
  {
    PINOCCHIO_CHECK(lanczos_decomposition.Qs().col(col_id).isZero(0));
  }

  const auto Qs_rank = lanczos_decomposition.Qs().leftCols(rank);
  PINOCCHIO_CHECK((Qs_rank.transpose() * Qs_rank).isIdentity());
}

BOOST_AUTO_TEST_CASE(test_random)
{
  typedef LanczosDecompositionTpl<Eigen::MatrixXd> LanczosDecomposition;
  const Eigen::DenseIndex mat_size = 20;

  for (int it = 0; it < 1000; ++it)
  {
    const Eigen::MatrixXd A = Eigen::MatrixXd::Random(mat_size, mat_size);
    const Eigen::MatrixXd matrix = A.transpose() * A;

    LanczosDecomposition lanczos_decomposition(matrix, 10);

    SET_LINE;
    checkDecomposition(lanczos_decomposition, matrix);
  }
}

BOOST_AUTO_TEST_CASE(test_low_rank)
{
  typedef LanczosDecompositionTpl<Eigen::MatrixXd> LanczosDecomposition;
  const Eigen::DenseIndex mat_size = 20;
  Eigen::MatrixXd A = Eigen::MatrixXd::Identity(mat_size, mat_size);
  A.row(mat_size - 1).setZero();
  A.col(mat_size - 1).setZero();

  LanczosDecomposition lanczos_decomposition_10(A, 10);
  SET_LINE;
  checkDecomposition(lanczos_decomposition_10, A, false);
}

BOOST_AUTO_TEST_CASE(test_delassus)
{
  typedef LanczosDecompositionTpl<Eigen::MatrixXd> LanczosDecomposition;
  const Eigen::DenseIndex mat_size = 20;

  for (int it = 0; it < 1000; ++it)
  {
    const Eigen::MatrixXd A = Eigen::MatrixXd::Random(mat_size, mat_size);
    const Eigen::MatrixXd matrix = A.transpose() * A;

    DelassusOperatorDense delassus(matrix);

    LanczosDecomposition lanczos_decomposition(delassus, 10);

    SET_LINE;
    checkDecomposition(lanczos_decomposition, matrix);
  }
}

void buildStackOfCubeModel(
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

BOOST_AUTO_TEST_CASE(test_delassus_cube)
{
  typedef LanczosDecompositionTpl<Eigen::MatrixXd> LanczosDecomposition;
  ::pinocchio::Model model;
  typedef FrictionalPointConstraintModel ConstraintModel;
  typedef typename ConstraintModel::ConstraintData ConstraintData;
  std::vector<ConstraintModel> constraint_models;

  const double box_mass = 10;
  std::vector<double> masses = {box_mass};

  buildStackOfCubeModel(masses, model, constraint_models);

  BOOST_CHECK(model.check(model.createData()));

  Data data(model);

  Eigen::VectorXd q0 = neutral(model);

  crba(model, data, q0, Convention::WORLD);
  ContactCholeskyDecomposition chol(model, constraint_models);
  std::vector<ConstraintData> constraint_datas;
  for (const auto & cm : constraint_models)
  {
    constraint_datas.push_back(cm.createData());
  }
  chol.compute(model, data, constraint_models, constraint_datas, 1e-10);

  const Eigen::MatrixXd delassus_matrix_plain = chol.getDelassusCholeskyExpression().matrix(true);
  auto G_expression = chol.getDelassusCholeskyExpression();

  BOOST_CHECK(delassus_matrix_plain.isApprox(delassus_matrix_plain.transpose(), 0));

  {
    LanczosDecomposition lanczos_decomposition(G_expression, 3);
    SET_LINE;
    checkDecomposition(lanczos_decomposition, delassus_matrix_plain);
  }

  {
    LanczosDecomposition lanczos_decomposition(G_expression, 4);
    SET_LINE;
    checkDecomposition(lanczos_decomposition, delassus_matrix_plain);
  }
}

BOOST_AUTO_TEST_CASE(test_delassus_light_cube)
{
  typedef LanczosDecompositionTpl<Eigen::MatrixXd> LanczosDecomposition;
  ::pinocchio::Model model;
  typedef FrictionalPointConstraintModel ConstraintModel;
  typedef typename ConstraintModel::ConstraintData ConstraintData;
  std::vector<ConstraintModel> constraint_models;

  const double box_mass = 1e-3;
  std::vector<double> masses = {box_mass};

  buildStackOfCubeModel(masses, model, constraint_models);

  BOOST_CHECK(model.check(model.createData()));

  Data data(model);

  Eigen::VectorXd q0 = neutral(model);

  crba(model, data, q0, Convention::WORLD);
  ContactCholeskyDecomposition chol(model, constraint_models);
  std::vector<ConstraintData> constraint_datas;
  for (const auto & cm : constraint_models)
  {
    constraint_datas.push_back(cm.createData());
  }
  chol.compute(model, data, constraint_models, constraint_datas, 1e-10);

  const Eigen::MatrixXd delassus_matrix_plain = chol.getDelassusCholeskyExpression().matrix(true);
  auto G_expression = chol.getDelassusCholeskyExpression();

  BOOST_CHECK(delassus_matrix_plain.isApprox(delassus_matrix_plain.transpose(), 1e-12));

  {
    LanczosDecomposition lanczos_decomposition(G_expression, 3);
    SET_LINE;
    checkDecomposition(lanczos_decomposition, delassus_matrix_plain);
  }
  BOOST_CHECK(delassus_matrix_plain.isApprox(delassus_matrix_plain.transpose(), 0));

  {
    LanczosDecomposition lanczos_decomposition(G_expression, 4);
    SET_LINE;
    checkDecomposition(lanczos_decomposition, delassus_matrix_plain);
  }
}

BOOST_AUTO_TEST_SUITE_END()
