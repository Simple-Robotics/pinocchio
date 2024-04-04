//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_contact_inverse_dynamics__hpp__
#define __pinocchio_algorithm_contact_inverse_dynamics__hpp__

#include "pinocchio/context.hpp"
#include "pinocchio/context/generic.hpp"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/constraints/coulomb-friction-cone.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include <Eigen/src/Core/util/Meta.h>
#include <boost/optional/optional.hpp>
#include <pinocchio/algorithm/contact-cholesky.hpp>
#include <pinocchio/algorithm/contact-jacobian.hpp>
#include "pinocchio/algorithm/proximal.hpp"

#include <boost/optional.hpp>
  
namespace pinocchio
{   

  ///
  /// \brief Compute the contact impulses given a target velocity of contact points.
  ///
  /// \tparam JointCollection Collection of Joint types.
  /// \tparam ConfigVectorType Type of the joint configuration vector.
  /// \tparam TangentVectorType1 Type of the joint velocity vector.
  /// \tparam TangentVectorType2 Type of the joint acceleration vector.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] c_ref The contact point velocity
  /// \param[in] contact_models The list of contact models.
  /// \param[in] contact_datas The list of contact_datas.
  /// \param[in] cones list of friction cones.
  /// \param[in] R vector representing the diagonal of the compliance matrix.
  /// \param[in] constraint_correction vector representing the constraint correction.
  /// \param[in] settings The settings for the proximal algorithm.
  /// \param[in] impulse_guess initial guess for the contact impulses.
  ///
  /// \return The desired joint torques stored in data.tau.
  ///
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, typename VectorLikeC, typename VectorLikeR,typename VectorLikeImp>
  const typename DataTpl<Scalar,Options,JointCollectionTpl>::TangentVectorType &
  computeContactImpulses(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
       DataTpl<Scalar,Options,JointCollectionTpl> & data,
       const Eigen::MatrixBase<VectorLikeC> & c_ref,
       const std::vector<RigidConstraintModelTpl<Scalar,Options>,ConstraintModelAllocator> & contact_models,
       std::vector<RigidConstraintDataTpl<Scalar,Options>,ConstraintDataAllocator> & contact_datas,
       const std::vector<CoulombFrictionConeTpl<Scalar>,CoulombFrictionConelAllocator> & cones,
       const Eigen::MatrixBase<VectorLikeR> & R,
      //  const Eigen::MatrixBase<VectorLikeGamma> & constraint_correction,
       ProximalSettingsTpl<Scalar> & settings,
       const boost::optional<VectorLikeImp > &impulse_guess= boost::none){

    const std::size_t problem_size  = R.size();
    const std::size_t n_contacts = cones.size();
    // PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_correction.size(), problem_size);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(contact_models.size(), n_contacts);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(contact_datas.size(), n_contacts);
    PINOCCHIO_CHECK_INPUT_ARGUMENT(check_expression_if_real<Scalar>(settings.mu > Scalar(0)),
                                   "mu has to be strictly positive");
    context::VectorXs R_prox; // TODO: malloc
    R_prox = R + context::VectorXs::Constant(problem_size,settings.mu);
    if(impulse_guess)
    {
      data.impulse_c = impulse_guess.get();
      PINOCCHIO_CHECK_ARGUMENT_SIZE(data.impulse_c.size(), problem_size);
    }
    else
    {
      data.impulse_c.setZero();
    }
    Scalar impulse_c_prev_norm_inf = data.impulse_c.template lpNorm<Eigen::Infinity>();
    Scalar complementarity, dual_feasibility;
    bool abs_prec_reached = false, rel_prec_reached = false;
    settings.iter = 1;
    for(; settings.iter <= settings.max_iter; ++settings.iter)
    {
      settings.absolute_residual = Scalar(0);
      settings.relative_residual = Scalar(0);
      for(std::size_t cone_id = 0; cone_id < n_contacts; ++cone_id)
      {
        const Eigen::DenseIndex row_id = 3*cone_id;
        const CoulombFrictionCone & cone = cones[cone_id];
        context::Vector3 impulse_c_prev = data.impulse_c.template segment<3>(row_id);
        auto impulse_segment = data.impulse_c.template segment<3>(row_id);
        auto R_prox_segment = R_prox.template segment<3>(row_id);
        auto c_ref_segment = c_ref.template segment<3>(row_id);
        context::Vector3 desaxce_segment = cone.computeNormalCorrection(c_ref_segment + (R.template segment<3>(row_id).array()*impulse_segment.array()).matrix());
        context::Vector3 c_cor_segment = c_ref_segment + desaxce_segment;
        impulse_segment = -((c_cor_segment -settings.mu * impulse_c_prev).array()/R_prox_segment.array()).matrix();
        impulse_segment = cone.weightedProject(impulse_segment, R_prox_segment);
        // evaluate convergence criteria
        Scalar contact_complementarity = cone.computeContactComplementarity(c_cor_segment, impulse_segment);
        settings.absolute_residual = math::max(settings.absolute_residual,contact_complementarity);
        context::Vector3 dimpulse_c = impulse_segment - impulse_c_prev;
        Scalar proximal_metric = dimpulse_c.template lpNorm<Eigen::Infinity>();
        settings.relative_residual = math::max(settings.relative_residual,proximal_metric);
      }
      
      const Scalar impulse_c_norm_inf = data.impulse_c.template lpNorm<Eigen::Infinity>();

      if(check_expression_if_real<Scalar,false>(settings.absolute_residual <= settings.absolute_accuracy))
        abs_prec_reached = true;
      else
        abs_prec_reached = false;

      if(check_expression_if_real<Scalar,false>(settings.relative_residual <= settings.relative_accuracy * math::max(impulse_c_norm_inf,impulse_c_prev_norm_inf)))
        rel_prec_reached = true;
      else
        rel_prec_reached = false;

      if(abs_prec_reached || rel_prec_reached)
        break;
      
      impulse_c_prev_norm_inf = impulse_c_norm_inf;
    }
    return data.impulse_c;
  }

  
  ///
  /// \brief The Contact Inverse Dynamics algorithm. It computes the inverse dynamics in the presence of contacts, aka the joint torques according to the current state of the system and the desired joint accelerations.
  ///
  /// \tparam JointCollection Collection of Joint types.
  /// \tparam ConfigVectorType Type of the joint configuration vector.
  /// \tparam TangentVectorType1 Type of the joint velocity vector.
  /// \tparam TangentVectorType2 Type of the joint acceleration vector.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] q The joint configuration vector (dim model.nq).
  /// \param[in] v The joint velocity vector (dim model.nv).
  /// \param[in] a The joint acceleration vector (dim model.nv).
  /// \param[in] dt The time step.
  /// \param[in] contact_models The list of contact models.
  /// \param[in] contact_datas The list of contact_datas.
  /// \param[in] cones list of friction cones.
  /// \param[in] R vector representing the diagonal of the compliance matrix.
  /// \param[in] constraint_correction vector representing the constraint correction.
  /// \param[in] settings The settings for the proximal algorithm.
  /// \param[in] lambda_guess initial guess for the contact forces.
  ///
  /// \return The desired joint torques stored in data.tau.
  ///
template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, typename ConfigVectorType, typename TangentVectorType1, typename TangentVectorType2, class ConstraintModelAllocator, class ConstraintDataAllocator, class CoulombFrictionConelAllocator, typename VectorLikeR, typename VectorLikeGamma,typename VectorLikeLam>
  const typename DataTpl<Scalar,Options,JointCollectionTpl>::TangentVectorType &
  contactInverseDynamics(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
       DataTpl<Scalar,Options,JointCollectionTpl> & data, 
       const Eigen::MatrixBase<ConfigVectorType> & q,
       const Eigen::MatrixBase<TangentVectorType1> & v,
       const Eigen::MatrixBase<TangentVectorType2> & a,
       Scalar dt,
       const std::vector<RigidConstraintModelTpl<Scalar,Options>,ConstraintModelAllocator> & contact_models,
       std::vector<RigidConstraintDataTpl<Scalar,Options>,ConstraintDataAllocator> & contact_datas,
       const std::vector<CoulombFrictionConeTpl<Scalar>,CoulombFrictionConelAllocator> & cones,
       const Eigen::MatrixBase<VectorLikeR> & R,
       const Eigen::MatrixBase<VectorLikeGamma> & constraint_correction,
       ProximalSettingsTpl<Scalar> & settings,
       const boost::optional<VectorLikeLam> &lambda_guess= boost::none){
    const std::size_t problem_size  = R.size();
    const std::size_t n_contacts = cones.size();
    context::MatrixXs J = context::MatrixXs::Zero(problem_size,model.nv); // TODO: malloc
    getConstraintsJacobian(model, data, contact_models, contact_datas, J);
    context::VectorXs v_ref, c_ref, tau_c;
    v_ref = v + dt*a;
    c_ref.noalias() = J* v_ref; //TODO should rather use the displacement
    c_ref += constraint_correction;
    boost::optional<context::VectorXs> impulse_guess = boost::none;
    if (lambda_guess){
      data.impulse_c = lambda_guess.get();
      data.impulse_c *= dt;
      impulse_guess = boost::make_optional(data.impulse_c);
    }
    computeContactImpulses(model, data, c_ref, contact_models, contact_datas, cones, R, settings, impulse_guess);
    // computeContactImpulses(model, data, c_ref, contact_models, contact_datas, cones, R, constraint_correction, settings, impulse_guess);
    data.lambda_c = data.impulse_c/dt;
    container::aligned_vector<context::Force> fext(model.njoints);
    for(std::size_t i = 0; i<model.njoints; i++){
      fext[i] = context::Force::Zero();
    }
    for(std::size_t i = 0; i<n_contacts; i++){
      const RigidConstraintModel & cmodel = contact_models[i];
      const Eigen::DenseIndex row_id = 3*i;
      auto lambda_segment = data.lambda_c.template segment<3>(row_id);
      typename RigidConstraintData::Matrix6 actInv_transpose1 = cmodel.joint1_placement.toActionMatrixInverse();
      actInv_transpose1.transposeInPlace();
      fext[cmodel.joint1_id] += Force(actInv_transpose1.template leftCols<3>() * lambda_segment);
      typename RigidConstraintData::Matrix6 actInv_transpose2 = cmodel.joint2_placement.toActionMatrixInverse();
      actInv_transpose2.transposeInPlace();
      fext[cmodel.joint2_id] += Force(actInv_transpose2.template leftCols<3>() * lambda_segment);
    }
    rnea(model, data, q, v, a, fext);
    return data.tau;
  }
  

} // namespace pinocchio 


// #if PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION
// #include "pinocchio/algorithm/contact-inverse-dynamics.txx"
// #endif // PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION

#endif // ifndef __pinocchio_algorithm_contact_inverse_dynamics_hpp__
