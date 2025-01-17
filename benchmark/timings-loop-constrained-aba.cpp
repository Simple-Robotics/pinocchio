//
// Copyright (c) 2024-2025 INRIA
//

#include "pinocchio/algorithm/contact-info.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
#include "pinocchio/algorithm/model.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/contact-dynamics.hpp"
#include "pinocchio/algorithm/constrained-dynamics.hpp"
#include "pinocchio/algorithm/cholesky.hpp"
#include "pinocchio/multibody/fwd.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/algorithm/pv.hpp"
#include "pinocchio/algorithm/loop-constrained-aba.hpp"

#include <iostream>
#include <numeric>

#include "pinocchio/utils/timer.hpp"

// usage ./benchmark/timings-CLCaba case_num prox-iters mu

template<typename T>
void print_first_n_vector_elements(std::vector<T> v, int n)
{
  for (int i = 0; i < n; i++)
  {
    std::cout << v[i] << ", ";
  }
  std::cout << "\n";
}

int main(int argc, const char ** argv)
{
  using namespace Eigen;
  using namespace pinocchio;

  PinocchioTicToc timer(PinocchioTicToc::US);
#ifdef NDEBUG
  const int NBT = 10000;
#else
  const int NBT = 10;
  std::cout << "(the time score in debug mode is not relevant) " << std::endl;
#endif

  // Build model
  Model model1;

  bool with_ff = false;
  int case_num = 0;

  if (argc > 1)
  {
    case_num = std::stoi(std::string{argv[1]});
    if (case_num >= 1)
      with_ff = true;
  }

  // std::string filename = PINOCCHIO_MODEL_DIR + std::string("/shadow_hand.urdf");
  std::string filename;
  if (case_num == 0)
    filename = PINOCCHIO_MODEL_DIR
               + std::string(
                 "/example-robot-data/robots/allegro_hand_description/urdf/allegro_left_hand.urdf");
  else if (case_num == 1)
    filename = PINOCCHIO_MODEL_DIR
               + std::string("/example-robot-data/robots/allegro_hand_description/urdf/"
                             "allegro_left_and_right_hands.urdf");
  else if (case_num >= 2)
    filename = PINOCCHIO_MODEL_DIR
               + std::string("/example-robot-data/robots/allegro_hand_description/urdf/"
                             "allegro_left_right_simple_humanoid.urdf");
  // std::string filename2 = PINOCCHIO_MODEL_DIR +
  // std::string("/example-robot-data/robots/allegro_hand_description/urdf/allegro_right_hand.urdf");

  // std::string filename3 = PINOCCHIO_MODEL_DIR + std::string("/simple_humanoid.urdf");
  // if(argc>2) filename = argv[2];

  if (with_ff)
    pinocchio::urdf::buildModel(filename, JointModelFreeFlyer(), model1);
  else
    pinocchio::urdf::buildModel(filename, model1);

  Model & model = model1;

  std::vector<std::string> allegro_hand_links{
    "link_3.0_tip", "link_7.0_tip", "link_11.0_tip", "link_15.0_tip"};
  std::vector<std::string> allegro_right_hand_links{
    "link_3.0_tip_r", "link_7.0_tip_r", "link_11.0_tip_r", "link_15.0_tip_r"};

  const std::string RF = allegro_hand_links[3];
  const JointIndex RF_id = model.frames[model.getFrameId(RF)].parent;
  const std::string MF = allegro_hand_links[2];
  const JointIndex MF_id = model.frames[model.getFrameId(MF)].parent;
  const std::string TF = allegro_hand_links[0];
  const JointIndex TF_id = model.frames[model.getFrameId(TF)].parent;
  const std::string IF = allegro_hand_links[1];
  const JointIndex IF_id = model.frames[model.getFrameId(IF)].parent;

  const std::string RF_r = allegro_right_hand_links[3];
  const JointIndex RF_id_r = model.frames[model.getFrameId(RF)].parent;
  const std::string MF_r = allegro_right_hand_links[2];
  const JointIndex MF_id_r = model.frames[model.getFrameId(MF)].parent;
  const std::string TF_r = allegro_right_hand_links[0];
  const JointIndex TF_id_r = model.frames[model.getFrameId(TF)].parent;
  const std::string IF_r = allegro_right_hand_links[1];
  const JointIndex IF_id_r = model.frames[model.getFrameId(IF)].parent;

  JointIndex cube_joint = model.addJoint(0, JointModelFreeFlyer(), SE3::Random(), "joint_cube");
  model.appendBodyToJoint(cube_joint, Inertia::Random(), SE3::Random());

  ProximalSettings prox_settings;
  if (argc > 2)
    prox_settings.max_iter = std::stoi(std::string{argv[2]});
  else
    prox_settings.max_iter = 20;

  prox_settings.mu = 1e-6;
  if (argc > 3)
    prox_settings.mu = std::pow(10, -1 * std::stoi(std::string{argv[3]}));

  prox_settings.relative_accuracy = 1e-18;
  prox_settings.absolute_accuracy = 1e-10;

  std::cout << "nq = " << model.nq << std::endl;
  std::cout << "nv = " << model.nv << std::endl;
  std::cout << "--" << std::endl;

  Data data(model), data_caba(model), data_caba_ref(model);
  VectorXd qmax = Eigen::VectorXd::Ones(model.nq);

  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(VectorXd) qs(NBT);
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(VectorXd) qdots(NBT);
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(VectorXd) qddots(NBT);
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(VectorXd) taus(NBT);

  for (size_t i = 0; i < NBT; ++i)
  {
    qs[i] = randomConfiguration(model, -qmax, qmax);
    qdots[i] = Eigen::VectorXd::Random(model.nv) * 0;
    qddots[i] = Eigen::VectorXd::Random(model.nv);
    taus[i] = Eigen::VectorXd::Random(model.nv);
  }

  timer.tic();
  SMOOTH(NBT)
  {
    aba(model, data, qs[_smooth], qdots[_smooth], taus[_smooth]);
  }

  double total_time = 0;

  RigidConstraintModel ci_CL(CONTACT_3D, model, TF_id, cube_joint, LOCAL_WORLD_ALIGNED);
  // ci_CL.joint1_placement.setRandom();
  ci_CL.joint2_placement.setRandom();

  RigidConstraintModel ci_CL2(CONTACT_3D, model, IF_id, cube_joint, LOCAL_WORLD_ALIGNED);
  // ci_CL2.joint1_placement.setRandom();
  ci_CL2.joint2_placement.setRandom();

  RigidConstraintModel ci_CL3(CONTACT_3D, model, MF_id, cube_joint, LOCAL_WORLD_ALIGNED);
  // ci_CL3.joint1_placement.setRandom();
  ci_CL3.joint2_placement.setRandom();

  RigidConstraintModel ci_CL4(CONTACT_3D, model, RF_id, cube_joint, LOCAL_WORLD_ALIGNED);
  // ci_CL4.joint1_placement.setRandom();
  ci_CL4.joint2_placement.setRandom();

  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel) contact_model_CL;
  PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData) contact_data_CL;

  if (case_num >= 0)
  {

    contact_model_CL.push_back(ci_CL);
    contact_data_CL.push_back(RigidConstraintData(ci_CL));
    contact_model_CL.push_back(ci_CL2);
    contact_data_CL.push_back(RigidConstraintData(ci_CL2));
    contact_model_CL.push_back(ci_CL3);
    contact_data_CL.push_back(RigidConstraintData(ci_CL3));
    contact_model_CL.push_back(ci_CL4);
    contact_data_CL.push_back(RigidConstraintData(ci_CL4));
  }

  if (case_num >= 1)
  {
    RigidConstraintModel ci_CL_r(CONTACT_3D, model, TF_id_r, cube_joint, LOCAL_WORLD_ALIGNED);
    // ci_CL_r.joint1_placement.setRandom();
    ci_CL_r.joint2_placement.setRandom();
    RigidConstraintModel ci_CL2_r(CONTACT_3D, model, IF_id_r, cube_joint, LOCAL_WORLD_ALIGNED);
    // ci_CL2_r.joint1_placement.setRandom();
    ci_CL2_r.joint2_placement.setRandom();
    RigidConstraintModel ci_CL3_r(CONTACT_3D, model, MF_id_r, cube_joint, LOCAL_WORLD_ALIGNED);
    // ci_CL3_r.joint1_placement.setRandom();
    ci_CL3_r.joint2_placement.setRandom();
    RigidConstraintModel ci_CL4_r(CONTACT_3D, model, RF_id_r, cube_joint, LOCAL_WORLD_ALIGNED);
    // ci_CL4_r.joint1_placement.setRandom();
    ci_CL4_r.joint2_placement.setRandom();
    contact_model_CL.push_back(ci_CL_r);
    contact_data_CL.push_back(RigidConstraintData(ci_CL_r));
    contact_model_CL.push_back(ci_CL2_r);
    contact_data_CL.push_back(RigidConstraintData(ci_CL2_r));
    contact_model_CL.push_back(ci_CL3_r);
    contact_data_CL.push_back(RigidConstraintData(ci_CL3_r));
    contact_model_CL.push_back(ci_CL4_r);
    contact_data_CL.push_back(RigidConstraintData(ci_CL4_r));
  }

  JointIndex Foot_id;
  // FEET constraints
  if (case_num >= 3)
  {
    const std::string RFoot = "RLEG_LINK6";
    const JointIndex RFoot_id = model.frames[model.getFrameId(RFoot)].parent;
    const std::string LFoot = "LLEG_LINK6";
    const JointIndex LFoot_id = model.frames[model.getFrameId(LFoot)].parent;
    Foot_id = RFoot_id;
    RigidConstraintModel ci_LFoot(CONTACT_6D, model, LFoot_id, LOCAL_WORLD_ALIGNED);
    RigidConstraintModel ci_RFoot(CONTACT_6D, model, RFoot_id, LOCAL_WORLD_ALIGNED);
    contact_model_CL.push_back(ci_LFoot);
    contact_data_CL.push_back(RigidConstraintData(ci_LFoot));
    contact_model_CL.push_back(ci_RFoot);
    contact_data_CL.push_back(RigidConstraintData(ci_LFoot));
  }

  std::vector<double> lcaba_residual(NBT);
  std::vector<double> proxLTL_residual(NBT);

  std::vector<double> lcaba_ddq_residual(NBT);
  std::vector<double> proxLTL_ddq_residual(NBT);

  std::vector<int> lcaba_iter_count(NBT);
  std::vector<int> ltl_iter_count(NBT);

  bool residual_benchmarking = false; // Set to false to do timings benchmarking
  if (residual_benchmarking)
  {
    // Investigate the convergence of the algorithms over different proximal iterations
    initLcaba(model, data_caba, contact_model_CL);
    initConstraintDynamics(model, data_caba_ref, contact_model_CL);
    initConstraintDynamics(model, data, contact_model_CL);

    double old_mu = prox_settings.mu;
    int old_max_iter = prox_settings.max_iter;

    for (size_t _smooth = 0; _smooth < NBT; _smooth++)
    {
      prox_settings.max_iter = old_max_iter;
      prox_settings.mu = old_mu;
      lcaba(
        model, data_caba, qs[_smooth], qdots[_smooth], taus[_smooth], contact_model_CL,
        contact_data_CL, prox_settings);
      lcaba_iter_count[_smooth] = prox_settings.iter;
      // std::cout << "lcaba iter count = " << prox_settings.iter << std::endl;

      // ProximalSettings prox_settings_ltl{prox_settings};
      constraintDynamics(
        model, data, qs[_smooth], qdots[_smooth], taus[_smooth], contact_model_CL, contact_data_CL,
        prox_settings);
      ltl_iter_count[_smooth] = prox_settings.iter;

      // std::cout << "proxLTL iter count = " << prox_settings.iter << std::endl;
      // std::cout << "Constraint forces =  " << data.lambda_c.transpose() << std::endl;

      prox_settings.max_iter = 50;
      prox_settings.mu = 1e-3;
      constraintDynamics(
        model, data_caba_ref, qs[_smooth], qdots[_smooth], taus[_smooth], contact_model_CL,
        contact_data_CL, prox_settings);

      long constraint_dim = data.lambda_c_prox.rows();
      const Eigen::MatrixXd & J_ref =
        data.contact_chol.matrix().topRightCorner(constraint_dim, model.nv);
      const Eigen::MatrixXd & rhs = data.primal_rhs_contact.topRows(constraint_dim);
      Eigen::MatrixXd constraint_error = J_ref * data.ddq - rhs;

      lcaba_residual[_smooth] = (J_ref * data_caba.ddq - rhs).template lpNorm<Eigen::Infinity>();
      proxLTL_residual[_smooth] = (J_ref * data.ddq - rhs).template lpNorm<Eigen::Infinity>();

      // std::cout << "ltl ddq = " << data.ddq.transpose() << std::endl;
      // std::cout << "caba - ltl ddq = " << (data.ddq - data_caba.ddq).template
      // lpNorm<Eigen::Infinity>() << std::endl; lcaba_residual[_smooth] =
      // (J_ref.transpose()*(J_ref*data_caba.ddq - rhs)).template lpNorm<Eigen::Infinity>();
      // proxLTL_residual[_smooth] = (J_ref.transpose()*(J_ref*data.ddq - rhs)).template
      // lpNorm<Eigen::Infinity>();

      lcaba_ddq_residual[_smooth] =
        (data_caba_ref.ddq - data_caba.ddq).template lpNorm<Eigen::Infinity>()
        / data_caba_ref.ddq.template lpNorm<Eigen::Infinity>();
      proxLTL_ddq_residual[_smooth] =
        (data_caba_ref.ddq - data.ddq).template lpNorm<Eigen::Infinity>()
        / data_caba_ref.ddq.template lpNorm<Eigen::Infinity>();
    }

    double lcaba_residual_sum = std::accumulate(lcaba_residual.begin(), lcaba_residual.end(), 0.0);
    double lcaba_residual_mean = lcaba_residual_sum / NBT;
    double lcaba_sq_sum =
      std::inner_product(lcaba_residual.begin(), lcaba_residual.end(), lcaba_residual.begin(), 0.0);
    double lcaba_residual_stddev =
      std::sqrt(lcaba_sq_sum / NBT - lcaba_residual_mean * lcaba_residual_mean);

    // std::cout << "LCABA residual mean = " << lcaba_residual_mean << ", and std = " <<
    // lcaba_residual_stddev << ", max = " << *std::max_element(lcaba_residual.begin(),
    // lcaba_residual.end()) << std::endl;

    // std::cout << "Constraint residual stats:" << lcaba_residual_mean << ", " <<
    // lcaba_residual_stddev << ", " << *std::max_element(lcaba_residual.begin(),
    // lcaba_residual.end()) << ", " << *std::min_element(lcaba_residual.begin(),
    // lcaba_residual.end());

    double proxLTL_residual_sum =
      std::accumulate(proxLTL_residual.begin(), proxLTL_residual.end(), 0.0);
    double proxLTL_residual_mean = proxLTL_residual_sum / NBT;
    double proxLTL_sq_sum = std::inner_product(
      proxLTL_residual.begin(), proxLTL_residual.end(), proxLTL_residual.begin(), 0.0);
    double proxLTL_residual_stddev =
      std::sqrt(proxLTL_sq_sum / NBT - proxLTL_residual_mean * proxLTL_residual_mean);

    // std::cout << ", " << proxLTL_residual_mean << ", " << proxLTL_residual_stddev << ", " <<
    // *std::max_element(proxLTL_residual.begin(), proxLTL_residual.end()) << ", " <<
    // *std::min_element(proxLTL_residual.begin(), proxLTL_residual.end()) << std::endl;

    double lcaba_ddq_residual_sum =
      std::accumulate(lcaba_ddq_residual.begin(), lcaba_ddq_residual.end(), 0.0);
    double lcaba_ddq_residual_mean = lcaba_ddq_residual_sum / NBT;
    lcaba_sq_sum = std::inner_product(
      lcaba_ddq_residual.begin(), lcaba_ddq_residual.end(), lcaba_ddq_residual.begin(), 0.0);
    double lcaba_ddq_residual_stddev =
      std::sqrt(lcaba_sq_sum / NBT - lcaba_ddq_residual_mean * lcaba_ddq_residual_mean);

    // std::cout << "ddq residual mean = " << lcaba_ddq_residual_mean << ", and std = " <<
    // lcaba_ddq_residual_stddev << ", max = " << *std::max_element(lcaba_ddq_residual.begin(),
    // lcaba_ddq_residual.end()) << std::endl;

    // std::cout << "LCABA constraint residuals (first 5)";
    // print_first_n_vector_elements(lcaba_residual, 5); std::cout << "LCABA ddq residuals (first
    // 5)"; print_first_n_vector_elements(lcaba_ddq_residual, 5); std::cout << "LCABA iter count
    // (first 5)"; print_first_n_vector_elements(lcaba_iter_count, 5); std::cout << "LCABA ddq
    // residuals (first 5)"; print_first_n_vector_elements(proxLTL_ddq_residual, 5);

    double proxLTL_ddq_residual_sum =
      std::accumulate(proxLTL_ddq_residual.begin(), proxLTL_ddq_residual.end(), 0.0);
    double proxLTL_ddq_residual_mean = proxLTL_ddq_residual_sum / NBT;
    proxLTL_sq_sum = std::inner_product(
      proxLTL_ddq_residual.begin(), proxLTL_ddq_residual.end(), proxLTL_ddq_residual.begin(), 0.0);
    double proxLTL_ddq_residual_stddev =
      std::sqrt(proxLTL_sq_sum / NBT - proxLTL_ddq_residual_mean * proxLTL_ddq_residual_mean);

    // std::cout << "proxLTL ddq residual mean = " << proxLTL_ddq_residual_mean << ", and std = " <<
    // proxLTL_ddq_residual_stddev << ", max = " << *std::max_element(proxLTL_ddq_residual.begin(),
    // proxLTL_ddq_residual.end()) << std::endl;

    std::cout << "ddq residual stats:" << lcaba_ddq_residual_mean << ", "
              << lcaba_ddq_residual_stddev << ", "
              << *std::max_element(lcaba_ddq_residual.begin(), lcaba_ddq_residual.end()) << ", "
              << *std::min_element(lcaba_ddq_residual.begin(), lcaba_ddq_residual.end()) << ", "
              << proxLTL_ddq_residual_mean << ", " << proxLTL_ddq_residual_stddev << ", "
              << *std::max_element(proxLTL_ddq_residual.begin(), proxLTL_ddq_residual.end()) << ", "
              << *std::min_element(proxLTL_ddq_residual.begin(), proxLTL_ddq_residual.end())
              << std::endl;
  }
  else
  {
    initLcaba(model, data_caba, contact_model_CL);
    timer.tic();
    SMOOTH(NBT)
    {
      lcaba(
        model, data_caba, qs[_smooth], qdots[_smooth], taus[_smooth], contact_model_CL,
        contact_data_CL, prox_settings);
    }
    std::cout << "CL-constrainedABA closed loops {6D} = \t";
    timer.toc(std::cout, NBT);

    initConstraintDynamics(model, data, contact_model_CL);
    timer.tic();
    SMOOTH(NBT)
    {
      constraintDynamics(
        model, data, qs[_smooth], qdots[_smooth], taus[_smooth], contact_model_CL, contact_data_CL,
        prox_settings);
    }
    std::cout << "constraintDynamics closed loops {6D} = \t";
    timer.toc(std::cout, NBT);

    std::cout << "--" << std::endl;

    long constraint_dim = data.lambda_c_prox.rows();
    const Eigen::MatrixXd & J_ref =
      data.contact_chol.matrix().topRightCorner(constraint_dim, model.nv);
    const Eigen::MatrixXd & rhs = data.primal_rhs_contact.topRows(constraint_dim);
    Eigen::VectorXd constraint_error = J_ref * data.ddq - rhs;
    std::cout << "Constraint residual LTL = " << constraint_error.template lpNorm<Eigen::Infinity>()
              << std::endl;
    std::cout << "Constraint residual LCABA = "
              << (J_ref * data_caba.ddq - rhs).template lpNorm<Eigen::Infinity>() << std::endl;

    std::cout << "constraint accelerations = " << (J_ref * data.ddq).transpose() << std::endl;

    std::cout << "Gravity = " << model.gravity.linear().transpose() << std::endl;

    std::cout << "--" << std::endl;
  }

  return 0;
}
