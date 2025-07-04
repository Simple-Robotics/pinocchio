//
// Copyright (c) 2015-2018 CNRS
// Copyright (c) 2018-2025 INRIA
// Copyright (c) 2015 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#ifndef __pinocchio_multibody_data_hxx__
#define __pinocchio_multibody_data_hxx__

#include "pinocchio/spatial/fwd.hpp"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/utils/string-generator.hpp"
#include "pinocchio/multibody/liegroup/liegroup-algo.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

/// @cond DEV

namespace pinocchio
{

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  DataTpl<Scalar, Options, JointCollectionTpl>::DataTpl(const Model & model)
  : q_in(neutral(model))
  , v_in(VectorXs::Zero(model.nv))
  , a_in(VectorXs::Zero(model.nv))
  , tau_in(VectorXs::Zero(model.nv))
  , a((std::size_t)model.njoints, Motion::Zero())
  , oa((std::size_t)model.njoints, Motion::Zero())
  , oa_drift((std::size_t)model.njoints, Motion::Zero())
  , oa_augmented((std::size_t)model.njoints, Motion::Zero())
  , a_gf((std::size_t)model.njoints, Motion::Zero())
  , oa_gf((std::size_t)model.njoints, Motion::Zero())
  , v((std::size_t)model.njoints, Motion::Zero())
  , ov((std::size_t)model.njoints, Motion::Zero())
  , f((std::size_t)model.njoints, Force::Zero())
  , of((std::size_t)model.njoints, Force::Zero())
  , of_augmented((std::size_t)model.njoints, Force::Zero())
  , h((std::size_t)model.njoints, Force::Zero())
  , oh((std::size_t)model.njoints, Force::Zero())
  , oMi((std::size_t)model.njoints, SE3::Identity())
  , liMi((std::size_t)model.njoints, SE3::Identity())
  , tau(VectorXs::Zero(model.nv))
  , nle(VectorXs::Zero(model.nv))
  , g(VectorXs::Zero(model.nv))
  , oMf((std::size_t)model.nframes, SE3::Identity())
  , Ycrb((std::size_t)model.njoints, Inertia::Zero())
  , dYcrb((std::size_t)model.njoints, Inertia::Zero())
  , M(MatrixXs::Zero(model.nv, model.nv))
  , Minv(MatrixXs::Zero(model.nv, model.nv))
  , C(MatrixXs::Zero(model.nv, model.nv))
  , dHdq(Matrix6x::Zero(6, model.nv))
  , dFdq(Matrix6x::Zero(6, model.nv))
  , dFdv(Matrix6x::Zero(6, model.nv))
  , dFda(Matrix6x::Zero(6, model.nv))
  , SDinv(Matrix6x::Zero(6, model.nv))
  , UDinv(Matrix6x::Zero(6, model.nv))
  , IS(MatrixXs::Zero(6, model.nv))
  , vxI((std::size_t)model.njoints, Inertia::Matrix6::Zero())
  , Ivx((std::size_t)model.njoints, Inertia::Matrix6::Zero())
  , B((std::size_t)model.njoints, Inertia::Matrix6::Zero())
  , oinertias((std::size_t)model.njoints, Inertia::Zero())
  , oYcrb((std::size_t)model.njoints, Inertia::Zero())
  , doYcrb((std::size_t)model.njoints, Inertia::Matrix6::Zero())
  , ddq(VectorXs::Zero(model.nv))
  , Yaba((std::size_t)model.njoints, Inertia::Matrix6::Zero())
  , oYaba((std::size_t)model.njoints, Inertia::Matrix6::Zero())
  , oYaba_augmented((std::size_t)model.njoints, Inertia::Matrix6::Zero())
  , oL((std::size_t)model.njoints, Inertia::Matrix6::Zero())
  , oK((std::size_t)model.njoints, Inertia::Matrix6::Zero())
  , u(VectorXs::Zero(model.nv))
  , Ag(Matrix6x::Zero(6, model.nv))
  , dAg(Matrix6x::Zero(6, model.nv))
  , hg(Force::Zero())
  , dhg(Force::Zero())
  , Ig(Inertia::Zero())
  , Fcrb((std::size_t)model.njoints, Matrix6x::Zero(6, model.nv))
  , nvSubtree((std::size_t)model.njoints, -1)
  , start_idx_v_fromRow((std::size_t)model.nvExtended, -1)
  , end_idx_v_fromRow((std::size_t)model.nvExtended, -1)
  , idx_vExtended_to_idx_v_fromRow((std::size_t)model.nvExtended, -1)
  , U(MatrixXs::Identity(model.nv, model.nv))
  , D(VectorXs::Zero(model.nv))
  , Dinv(VectorXs::Zero(model.nv))
  , tmp(VectorXs::Zero(model.nv))
  , parents_fromRow((std::size_t)model.nvExtended, -1)
  , mimic_parents_fromRow((std::size_t)model.nvExtended, -1)
  , non_mimic_parents_fromRow((std::size_t)model.nvExtended, -1)
  , supports_fromRow((std::size_t)model.nv)
  , nvSubtree_fromRow((std::size_t)model.nvExtended, -1)
  , J(Matrix6x::Zero(6, model.nvExtended))
  , dJ(Matrix6x::Zero(6, model.nvExtended))
  , ddJ(Matrix6x::Zero(6, model.nvExtended))
  , psid(Matrix6x::Zero(6, model.nv))
  , psidd(Matrix6x::Zero(6, model.nv))
  , dVdq(Matrix6x::Zero(6, model.nv))
  , dAdq(Matrix6x::Zero(6, model.nv))
  , dAdv(Matrix6x::Zero(6, model.nv))
  , dtau_dq(RowMatrixXs::Zero(model.nv, model.nv))
  , dtau_dv(RowMatrixXs::Zero(model.nv, model.nv))
  , ddq_dq(RowMatrixXs::Zero(model.nv, model.nv))
  , ddq_dv(RowMatrixXs::Zero(model.nv, model.nv))
  , ddq_dtau(RowMatrixXs::Zero(model.nv, model.nv))
  , iMf((std::size_t)model.njoints, SE3::Identity())
  , com((std::size_t)model.njoints, Vector3::Zero())
  , vcom((std::size_t)model.njoints, Vector3::Zero())
  , acom((std::size_t)model.njoints, Vector3::Zero())
  , mass((std::size_t)model.njoints, (Scalar)(-1))
  , Jcom(Matrix3x::Zero(3, model.nv))
  , kinetic_energy(Scalar(0))
  , potential_energy(Scalar(0))
  , mechanical_energy(Scalar(0))
  , JMinvJt()
  , llt_JMinvJt()
  , lambda_c()
  , lambda_c_prox()
  , diff_lambda_c()
  , sDUiJt(MatrixXs::Zero(model.nv, model.nv))
  , torque_residual(VectorXs::Zero(model.nv))
  , dq_after(VectorXs::Zero(model.nv))
  , impulse_c()
  , staticRegressor(Matrix3x::Zero(3, 4 * (model.njoints - 1)))
  , bodyRegressor(BodyRegressorType::Zero())
  , jointTorqueRegressor(MatrixXs::Zero(model.nv, 10 * (model.njoints - 1)))
  , kineticEnergyRegressor(RowVectorXs::Zero(10 * (model.njoints - 1)))
  , potentialEnergyRegressor(RowVectorXs::Zero(10 * (model.njoints - 1)))
  , KA((std::size_t)model.njoints, Matrix6x::Zero(6, 0))
  , LA((std::size_t)model.njoints, MatrixXs::Zero(0, 0))
  , lA((std::size_t)model.njoints, VectorXs::Zero(0))
  , lambdaA((std::size_t)model.njoints, VectorXs::Zero(0))
  , par_cons_ind((std::size_t)model.njoints, 0)
  , a_bias((std::size_t)model.njoints, Motion::Zero())
  , KAS((std::size_t)model.njoints, MatrixXs::Zero(0, 0))

#if EIGEN_VERSION_AT_LEAST(3, 2, 90) && !EIGEN_VERSION_AT_LEAST(3, 2, 93)
  , kinematic_hessians(
      6,
      std::max(1, model.nv),
      std::max(1, model.nv)) // the minimum size should be 1 for compatibility reasons
  , d2tau_dqdq(
      std::max(1, model.nv),
      std::max(1, model.nv),
      std::max(1, model.nv)) // the minimum size should be 1 for compatibility reasons
  , d2tau_dvdv(
      std::max(1, model.nv),
      std::max(1, model.nv),
      std::max(1, model.nv)) // the minimum size should be 1 for compatibility reasons
  , d2tau_dqdv(
      std::max(1, model.nv),
      std::max(1, model.nv),
      std::max(1, model.nv)) // the minimum size should be 1 for compatibility reasons
  , d2tau_dadq(
      std::max(1, model.nv),
      std::max(1, model.nv),
      std::max(1, model.nv)) // the minimum size should be 1 for compatibility reasons
#else
  , kinematic_hessians(6, model.nv, model.nv)
  , d2tau_dqdq(model.nv, model.nv, model.nv)
  , d2tau_dvdv(model.nv, model.nv, model.nv)
  , d2tau_dqdv(model.nv, model.nv, model.nv)
  , d2tau_dadq(model.nv, model.nv, model.nv)
#endif
  , extended_motion_propagator((std::size_t)model.njoints, Matrix6::Zero())
  , extended_motion_propagator2((std::size_t)model.njoints, Matrix6::Zero())
  , spatial_inv_inertia((std::size_t)model.njoints, Matrix6::Zero())
  , accumulation_descendant((std::size_t)model.njoints, 0)
  , accumulation_ancestor((std::size_t)model.njoints, 0)
  , constraints_supported_dim((std::size_t)model.njoints, 0)
  , constraints_supported((std::size_t)model.njoints)
  , constraints_on_joint((std::size_t)model.njoints)
  , neighbour_links((std::size_t)model.njoints)
  , joint_cross_coupling(model.njoints, model.njoints)
  , joint_apparent_inertia(VectorXs::Zero(model.nv))
  {
    typedef typename Model::JointIndex JointIndex;

    /* Create data structure associated to the joints */
    joints.reserve(JointIndex(model.njoints));
    for (JointIndex i = 0; i < JointIndex(model.njoints); ++i)
    {
      typedef CreateJointData<Scalar, Options, JointCollectionTpl> Constructor;
      joints.push_back(Constructor::run(model.joints[i]));
    }
    joints_augmented = joints;

    /* Init for CRBA */
    M.setZero();
    Minv.setZero();

    computeNvSubtree(model);

    /* Init for Cholesky */
    computeParents_fromRow(model);
    computeSupports_fromRow(model);

    /* Init universe states relatively to itself */
    a_gf[0] = -model.gravity;

    kinematic_hessians.setZero();
    d2tau_dqdq.setZero();
    d2tau_dvdv.setZero();
    d2tau_dqdv.setZero();
    d2tau_dadq.setZero();
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline void DataTpl<Scalar, Options, JointCollectionTpl>::computeNvSubtree(const Model & model)
  {
    for (JointIndex joint_id = 0; joint_id < JointIndex(model.njoints); ++joint_id)
    {
      // Build a "correct" representation of mimic nvSubtree by using nvExtended, which will cover
      // its children nv, and allow for a simple check
      if (boost::get<JointModelMimicTpl<Scalar, Options, JointCollectionTpl>>(
            &model.joints[joint_id]))
        nvSubtree[joint_id] = 0;
      else
      {
        int nv_;
        const JointIndex last_child =
          model.subtrees[joint_id].size() > 0 ? model.subtrees[joint_id].back() : JointIndex(0);
        if (boost::get<JointModelMimicTpl<Scalar, Options, JointCollectionTpl>>(
              &model.joints[last_child]))
          nv_ = model.joints[last_child].nvExtended();
        else
          nv_ = model.joints[last_child].nv();
        nvSubtree[joint_id] =
          model.joints[last_child].idx_v() + nv_ - model.joints[joint_id].idx_v();
      }
    }
    // fill mimic data
    for (const JointIndex mimicking_id : model.mimicking_joints)
    {
      const auto & mimicking_sub = model.subtrees[mimicking_id];
      size_t j = 1;
      for (; j < mimicking_sub.size(); j++)
      {
        if (model.nvs[mimicking_sub[j]] != 0)
          break;
      }
      if (mimicking_sub.size() == 1)
        mimic_subtree_joint.push_back(0);
      else
        mimic_subtree_joint.push_back(mimicking_sub[j]);
    }
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline void
  DataTpl<Scalar, Options, JointCollectionTpl>::computeParents_fromRow(const Model & model)
  {
    typedef typename Model::Index Index;

    for (Index joint = 1; joint < (Index)(model.njoints); joint++)
    {
      const Index & parent = model.parents[joint];
      const int idx_vj = model.joints[joint].idx_v();
      const int nvExtended_j = model.joints[joint].nvExtended();
      const int idx_vExtended_j = model.joints[joint].idx_vExtended();

      assert(idx_vExtended_j >= 0 && idx_vExtended_j < model.nvExtended);
      assert(idx_vj >= 0 && idx_vj < model.nv);

      if (parent > 0)
        parents_fromRow[(Index)idx_vExtended_j] =
          model.joints[parent].idx_vExtended() + model.joints[parent].nvExtended() - 1;
      else
        parents_fromRow[(Index)idx_vExtended_j] = -1;

      JointIndex first_non_mimic_parent_id = parent;
      while (first_non_mimic_parent_id > 0
             && boost::get<JointModelMimicTpl<Scalar, Options, JointCollectionTpl>>(
               &model.joints[first_non_mimic_parent_id]))
      {
        first_non_mimic_parent_id = model.parents[first_non_mimic_parent_id];
      }

      if (first_non_mimic_parent_id > 0)
        non_mimic_parents_fromRow[(Index)idx_vExtended_j] =
          model.joints[first_non_mimic_parent_id].idx_vExtended()
          + model.joints[first_non_mimic_parent_id].nvExtended() - 1;
      else
        non_mimic_parents_fromRow[(Index)idx_vExtended_j] = -1;

      JointIndex first_mimic_parent_id = parent;
      while (first_mimic_parent_id > 0
             && !boost::get<JointModelMimicTpl<Scalar, Options, JointCollectionTpl>>(
               &model.joints[first_mimic_parent_id]))
      {
        first_mimic_parent_id = model.parents[first_mimic_parent_id];
      }

      if (first_mimic_parent_id > 0)
        mimic_parents_fromRow[(Index)idx_vExtended_j] =
          model.joints[first_mimic_parent_id].idx_vExtended()
          + model.joints[first_mimic_parent_id].nvExtended() - 1;
      else
        mimic_parents_fromRow[(Index)idx_vExtended_j] = -1;

      nvSubtree_fromRow[(Index)idx_vExtended_j] = nvSubtree[joint];
      start_idx_v_fromRow[(size_t)idx_vj] = idx_vj;
      end_idx_v_fromRow[(size_t)idx_vj] = idx_vj + nvExtended_j - 1;
      idx_vExtended_to_idx_v_fromRow[(size_t)idx_vExtended_j] = idx_vj;

      for (int row = 1; row < nvExtended_j; ++row)
      {
        parents_fromRow[(size_t)(idx_vExtended_j + row)] = idx_vExtended_j + row - 1;
        mimic_parents_fromRow[(size_t)(idx_vExtended_j + row)] = idx_vExtended_j + row - 1;
        non_mimic_parents_fromRow[(size_t)(idx_vExtended_j + row)] = idx_vExtended_j + row - 1;
        nvSubtree_fromRow[(size_t)(idx_vExtended_j + row)] = nvSubtree[joint] - row;
        start_idx_v_fromRow[(size_t)(idx_vExtended_j + row)] =
          start_idx_v_fromRow[(size_t)idx_vExtended_j];
        end_idx_v_fromRow[(size_t)(idx_vExtended_j + row)] =
          end_idx_v_fromRow[(size_t)idx_vExtended_j];
        idx_vExtended_to_idx_v_fromRow[(size_t)(idx_vExtended_j + row)] =
          idx_vExtended_to_idx_v_fromRow[(size_t)idx_vExtended_j] + row;
      }
    }
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline void
  DataTpl<Scalar, Options, JointCollectionTpl>::computeSupports_fromRow(const Model & model)
  {
    typedef typename Model::JointIndex JointIndex;

    for (JointIndex joint_id = 1; joint_id < (JointIndex)(model.njoints); joint_id++)
    {
      const int nvj = nv(model.joints[joint_id]);
      const int idx_vj = idx_v(model.joints[joint_id]);

      assert(idx_vj >= 0 && idx_vj < model.nv);

      const int parent_fromRow = parents_fromRow[(size_t)idx_vj];

      if (parent_fromRow >= 0)
        supports_fromRow[(size_t)idx_vj] = supports_fromRow[(size_t)parent_fromRow];

      supports_fromRow[(size_t)idx_vj].push_back(idx_vj);

      for (int row = 1; row < nvj; ++row)
      {
        supports_fromRow[(size_t)(idx_vj + row)] = supports_fromRow[(size_t)(idx_vj + row - 1)];
        supports_fromRow[(size_t)(idx_vj + row)].push_back(idx_vj + row);
      }
    }
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  bool operator==(
    const DataTpl<Scalar, Options, JointCollectionTpl> & data1,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data2)
  {
    bool value =
      data1.joints == data2.joints && data1.joints_augmented == data2.joints_augmented
      && data1.q_in == data2.q_in && data1.v_in == data2.v_in && data1.a_in == data2.a_in
      && data1.tau_in == data2.tau_in && data1.a == data2.a && data1.oa == data2.oa
      && data1.oa_drift == data2.oa_drift && data1.oa_augmented == data2.oa_augmented
      && data1.a_gf == data2.a_gf && data1.oa_gf == data2.oa_gf && data1.v == data2.v
      && data1.ov == data2.ov && data1.f == data2.f && data1.of == data2.of
      && data1.of_augmented == data2.of_augmented && data1.h == data2.h && data1.oh == data2.oh
      && data1.oMi == data2.oMi && data1.liMi == data2.liMi && data1.tau == data2.tau
      && data1.nle == data2.nle && data1.g == data2.g && data1.oMf == data2.oMf
      && data1.Ycrb == data2.Ycrb && data1.dYcrb == data2.dYcrb && data1.M == data2.M
      && data1.Minv == data2.Minv && data1.C == data2.C && data1.dHdq == data2.dHdq
      && data1.dFdq == data2.dFdq && data1.dFdv == data2.dFdv && data1.dFda == data2.dFda
      && data1.SDinv == data2.SDinv && data1.UDinv == data2.UDinv && data1.IS == data2.IS
      && data1.vxI == data2.vxI && data1.Ivx == data2.Ivx && data1.oinertias == data2.oinertias
      && data1.oYcrb == data2.oYcrb && data1.doYcrb == data2.doYcrb && data1.ddq == data2.ddq
      && data1.Yaba == data2.Yaba && data1.oYaba == data2.oYaba
      && data1.oYaba_augmented == data2.oYaba_augmented && data1.oL == data2.oL
      && data1.oK == data2.oK && data1.u == data2.u && data1.Ag == data2.Ag
      && data1.dAg == data2.dAg && data1.hg == data2.hg && data1.dhg == data2.dhg
      && data1.Ig == data2.Ig && data1.Fcrb == data2.Fcrb && data1.nvSubtree == data2.nvSubtree
      && data1.start_idx_v_fromRow == data2.start_idx_v_fromRow
      && data1.end_idx_v_fromRow == data2.end_idx_v_fromRow && data1.U == data2.U
      && data1.D == data2.D && data1.Dinv == data2.Dinv
      && data1.parents_fromRow == data2.parents_fromRow
      && data1.mimic_parents_fromRow == data2.mimic_parents_fromRow
      && data1.non_mimic_parents_fromRow == data2.non_mimic_parents_fromRow
      && data1.idx_vExtended_to_idx_v_fromRow == data2.idx_vExtended_to_idx_v_fromRow
      && data1.mimic_subtree_joint == data2.mimic_subtree_joint
      && data1.supports_fromRow == data2.supports_fromRow
      && data1.nvSubtree_fromRow == data2.nvSubtree_fromRow && data1.J == data2.J
      && data1.dJ == data2.dJ && data1.ddJ == data2.ddJ && data1.psid == data2.psid
      && data1.psidd == data2.psidd && data1.dVdq == data2.dVdq && data1.dAdq == data2.dAdq
      && data1.dAdv == data2.dAdv && data1.dtau_dq == data2.dtau_dq
      && data1.dtau_dv == data2.dtau_dv && data1.ddq_dq == data2.ddq_dq
      && data1.ddq_dv == data2.ddq_dv && data1.dvc_dq == data2.dvc_dq
      && data1.dac_dq == data2.dac_dq && data1.dac_dv == data2.dac_dv
      && data1.dac_da == data2.dac_da && data1.osim == data2.osim
      && data1.dlambda_dq == data2.dlambda_dq && data1.dlambda_dv == data2.dlambda_dv
      && data1.dlambda_dtau == data2.dlambda_dtau && data1.dlambda_dx_prox == data2.dlambda_dx_prox
      && data1.drhs_prox == data2.drhs_prox && data1.iMf == data2.iMf && data1.com == data2.com
      && data1.vcom == data2.vcom && data1.acom == data2.acom && data1.mass == data2.mass
      && data1.Jcom == data2.Jcom && data1.kinetic_energy == data2.kinetic_energy
      && data1.potential_energy == data2.potential_energy
      && data1.mechanical_energy == data2.mechanical_energy && data1.JMinvJt == data2.JMinvJt
      && data1.lambda_c == data2.lambda_c && data1.lambda_c_prox == data2.lambda_c_prox
      && data1.diff_lambda_c == data2.diff_lambda_c && data1.sDUiJt == data2.sDUiJt
      && data1.torque_residual == data2.torque_residual && data1.dq_after == data2.dq_after
      && data1.impulse_c == data2.impulse_c && data1.staticRegressor == data2.staticRegressor
      && data1.bodyRegressor == data2.bodyRegressor
      && data1.jointTorqueRegressor == data2.jointTorqueRegressor
      //    && data1.contact_chol == data2.contact_chol
      && data1.primal_dual_contact_solution == data2.primal_dual_contact_solution
      && data1.extended_motion_propagator == data2.extended_motion_propagator
      && data1.extended_motion_propagator2 == data2.extended_motion_propagator2
      && data1.spatial_inv_inertia == data2.spatial_inv_inertia
      && data1.accumulation_descendant == data2.accumulation_descendant
      && data1.accumulation_ancestor == data2.accumulation_ancestor
      && data1.constraints_supported_dim == data2.constraints_supported_dim
      && data1.constraints_supported == data2.constraints_supported
      && data1.constraints_on_joint == data2.constraints_on_joint
      && data1.neighbour_links == data2.neighbour_links
      && data1.joint_cross_coupling == data2.joint_cross_coupling
      && data1.joint_apparent_inertia == data2.joint_apparent_inertia;

    // operator== for Eigen::Tensor provides an Expression which might be not evaluated as a boolean
    value &= Tensor<bool, 0>((data1.kinematic_hessians == data2.kinematic_hessians).all())(0)
             && Tensor<bool, 0>((data1.d2tau_dqdq == data2.d2tau_dqdq).all())(0)
             && Tensor<bool, 0>((data1.d2tau_dvdv == data2.d2tau_dvdv).all())(0)
             && Tensor<bool, 0>((data1.d2tau_dqdv == data2.d2tau_dqdv).all())(0)
             && Tensor<bool, 0>((data1.d2tau_dadq == data2.d2tau_dadq).all())(0);

    return value;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  bool operator!=(
    const DataTpl<Scalar, Options, JointCollectionTpl> & data1,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data2)
  {
    return !(data1 == data2);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::Data
  ModelTpl<Scalar, Options, JointCollectionTpl>::createData() const
  {
    return Data(*this);
  }

} // namespace pinocchio

/// @endcond

#endif // ifndef __pinocchio_multibody_data_hxx__
