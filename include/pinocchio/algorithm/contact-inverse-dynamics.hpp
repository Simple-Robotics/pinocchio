//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_contact_inverse_dynamics_hpp__
#define __pinocchio_algorithm_contact_inverse_dynamics_hpp__

#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/contact-jacobian.hpp"
#include "pinocchio/algorithm/proximal.hpp"

#include <boost/optional.hpp>

namespace pinocchio
{

  ///
  /// \brief Compute the contact forces given a target velocity of contact points.
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
  /// \param[in] R vector representing the diagonal of the compliance matrix.
  /// \param[in] constraint_correction vector representing the constraint correction.
  /// \param[in] settings The settings for the proximal algorithm.
  /// \param[in] lambda_guess initial guess for the contact forces.
  ///
  /// \return The desired joint torques stored in data.tau.
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename VectorLikeC,
    template<typename T> class Holder,
    class ConstraintModelAllocator,
    class ConstraintDataAllocator,
    typename VectorLikeR,
    typename VectorLikeImp>
  const typename DataTpl<Scalar, Options, JointCollectionTpl>::TangentVectorType &
  computeContactForces(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const Eigen::MatrixBase<VectorLikeC> & c_ref,
    const std::vector<
      Holder<const FrictionalPointConstraintModelTpl<Scalar, Options>>,
      ConstraintModelAllocator> & contact_models,
    std::vector<
      Holder<FrictionalPointConstraintDataTpl<Scalar, Options>>,
      ConstraintDataAllocator> & contact_datas,
    const Eigen::MatrixBase<VectorLikeR> & R,
    // const Eigen::MatrixBase<VectorLikeGamma> & constraint_correction,
    ProximalSettingsTpl<Scalar> & settings,
    const boost::optional<VectorLikeImp> & lambda_guess = boost::none)
  {
    PINOCCHIO_UNUSED_VARIABLE(model);
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    using VectorXs = typename Model::VectorXs;
    using Vector3 = typename Model::Vector3;
    typedef FrictionalPointConstraintModelTpl<Scalar, Options> ConstraintModel;

    const Eigen::Index problem_size = R.size();
    const std::size_t n_contacts = contact_models.size();
    // PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_correction.size(), problem_size);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(contact_models.size(), n_contacts);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(contact_datas.size(), n_contacts);
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      check_expression_if_real<Scalar>(settings.mu > Scalar(0)), "mu has to be strictly positive");

    const VectorXs R_prox = R + VectorXs::Constant(problem_size, settings.mu);
    if (lambda_guess)
    {
      data.lambda_c = lambda_guess.get();
      PINOCCHIO_CHECK_ARGUMENT_SIZE(data.lambda_c.size(), problem_size);
    }
    else
    {
      data.lambda_c.setZero();
    }

    Scalar lambda_c_prev_norm_inf = data.lambda_c.template lpNorm<Eigen::Infinity>();
    bool abs_prec_reached = false, rel_prec_reached = false;
    settings.iter = 1;
    for (; settings.iter <= settings.max_iter; ++settings.iter)
    {
      settings.absolute_residual = Scalar(0);
      settings.relative_residual = Scalar(0);
      Eigen::DenseIndex row_id = 0;
      for (std::size_t constraint_id = 0; constraint_id < n_contacts; ++constraint_id)
      {

        const ConstraintModel & cmodel = contact_models[constraint_id];
        const auto constraint_size = cmodel.size();

        const CoulombFrictionCone & cone = cmodel.set();
        const Vector3 lambda_c_previous = data.lambda_c.template segment<3>(row_id);
        auto lambda_segment = data.lambda_c.template segment<3>(row_id);
        const auto R_prox_segment = R_prox.template segment<3>(row_id);
        const auto c_ref_segment = c_ref.template segment<3>(row_id);
        const Vector3 sigma_segment =
          c_ref_segment + (R.template segment<3>(row_id).array() * lambda_segment.array()).matrix();
        const Vector3 desaxce_correction = cone.computeNormalCorrection(sigma_segment);
        const Vector3 c_cor_segment = c_ref_segment + desaxce_correction;
        lambda_segment =
          -((c_cor_segment - settings.mu * lambda_c_previous).array() / R_prox_segment.array())
             .matrix();
        lambda_segment = cone.weightedProject(Vector3(lambda_segment), R_prox_segment);

        // Compute convergence criteria
        Scalar contact_complementarity =
          cone.computeConicComplementarity(sigma_segment + desaxce_correction, lambda_segment);
        Scalar dual_feasibility =
          std::abs(math::min(0., sigma_segment(2))); // proxy of dual feasibility
        settings.absolute_residual = math::max(
          settings.absolute_residual, math::max(contact_complementarity, dual_feasibility));
        const Vector3 dlambda_c = lambda_segment - lambda_c_previous;
        Scalar proximal_metric = dlambda_c.template lpNorm<Eigen::Infinity>();
        settings.relative_residual = math::max(settings.relative_residual, proximal_metric);

        row_id += constraint_size;
      }

      const Scalar lambda_c_norm_inf = data.lambda_c.template lpNorm<Eigen::Infinity>();

      if (check_expression_if_real<Scalar, false>(
            settings.absolute_residual <= settings.absolute_accuracy))
        abs_prec_reached = true;
      else
        abs_prec_reached = false;

      if (check_expression_if_real<Scalar, false>(
            settings.relative_residual
            <= settings.relative_accuracy * math::max(lambda_c_norm_inf, lambda_c_prev_norm_inf)))
        rel_prec_reached = true;
      else
        rel_prec_reached = false;

      if (abs_prec_reached || rel_prec_reached)
        break;

      lambda_c_prev_norm_inf = lambda_c_norm_inf;
    }

    return data.lambda_c;
  }

  ///
  /// \brief Compute the contact forces given a target velocity of contact points.
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
  /// \param[in] R vector representing the diagonal of the compliance matrix.
  /// \param[in] constraint_correction vector representing the constraint correction.
  /// \param[in] settings The settings for the proximal algorithm.
  /// \param[in] lambda_guess initial guess for the contact forces.
  ///
  /// \return The desired joint torques stored in data.tau.
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename VectorLikeC,
    class ConstraintModelAllocator,
    class ConstraintDataAllocator,
    typename VectorLikeR,
    typename VectorLikeImp>
  const typename DataTpl<Scalar, Options, JointCollectionTpl>::TangentVectorType &
  computeContactForces(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const Eigen::MatrixBase<VectorLikeC> & c_ref,
    const std::vector<
      FrictionalPointConstraintModelTpl<Scalar, Options>,
      ConstraintModelAllocator> & contact_models,
    std::vector<FrictionalPointConstraintDataTpl<Scalar, Options>, ConstraintDataAllocator> &
      contact_datas,
    const Eigen::MatrixBase<VectorLikeR> & R,
    // const Eigen::MatrixBase<VectorLikeGamma> & constraint_correction,
    ProximalSettingsTpl<Scalar> & settings,
    const boost::optional<VectorLikeImp> & lambda_guess = boost::none)
  {
    typedef std::reference_wrapper<const FrictionalPointConstraintModelTpl<Scalar, Options>>
      WrappedConstraintModelType;
    typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

    WrappedConstraintModelVector wrapped_constraint_models(
      contact_models.cbegin(), contact_models.cend());

    typedef std::reference_wrapper<FrictionalPointConstraintDataTpl<Scalar, Options>>
      WrappedConstraintDataType;
    typedef std::vector<WrappedConstraintDataType> WrappedConstraintDataVector;

    WrappedConstraintDataVector wrapped_constraint_datas(
      contact_datas.begin(), contact_datas.end());

    return computeContactForces(
      model, data, c_ref, wrapped_constraint_models, wrapped_constraint_datas, R, settings,
      lambda_guess);
  }

  ///
  /// \brief The Contact Inverse Dynamics algorithm. It computes the inverse dynamics in the
  /// presence of contacts, aka the joint torques according to the current state of the system and
  /// the desired joint accelerations.
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
  /// \param[in] R vector representing the diagonal of the compliance matrix.
  /// \param[in] constraint_correction vector representing the constraint correction.
  /// \param[in] settings The settings for the proximal algorithm.
  /// \param[in] lambda_guess initial guess for the contact forces.
  ///
  /// \return The desired joint torques stored in data.tau.
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType,
    typename TangentVectorType1,
    typename TangentVectorType2,
    template<typename T> class Holder,
    class ConstraintModelAllocator,
    class ConstraintDataAllocator,
    typename VectorLikeR,
    typename VectorLikeGamma,
    typename VectorLikeLam>
  const typename DataTpl<Scalar, Options, JointCollectionTpl>::TangentVectorType &
  contactInverseDynamics(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<TangentVectorType1> & v,
    const Eigen::MatrixBase<TangentVectorType2> & a,
    Scalar dt,
    const std::vector<
      Holder<const FrictionalPointConstraintModelTpl<Scalar, Options>>,
      ConstraintModelAllocator> & contact_models,
    std::vector<
      Holder<FrictionalPointConstraintDataTpl<Scalar, Options>>,
      ConstraintDataAllocator> & contact_datas,
    const Eigen::MatrixBase<VectorLikeR> & R,
    const Eigen::MatrixBase<VectorLikeGamma> & constraint_correction,
    ProximalSettingsTpl<Scalar> & settings,
    const boost::optional<VectorLikeLam> & lambda_guess = boost::none)
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    //    typedef FrictionalPointConstraintDataTpl<Scalar, Options> ConstraintData;

    typedef typename Model::MatrixXs MatrixXs;
    typedef typename Model::VectorXs VectorXs;
    typedef typename Model::Force Force;

    const Eigen::Index problem_size = R.size();
    const std::size_t n_contacts = contact_models.size();
    MatrixXs J = MatrixXs::Zero(problem_size, model.nv); // TODO: malloc
    getConstraintsJacobian(model, data, contact_models, contact_datas, J);
    VectorXs v_ref, c_ref, tau_c;
    v_ref.noalias() = v + dt * a;
    c_ref.noalias() = J * v_ref; // TODO should rather use the displacement
    c_ref += constraint_correction;
    c_ref /= dt; // we work with a formulation on forces
    computeContactForces(
      model, data, c_ref, contact_models, contact_datas, R, settings, lambda_guess);
    PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(Force)
    fext(std::size_t(model.njoints), Force::Zero());

    Eigen::DenseIndex row_id = 0;
    for (std::size_t i = 0; i < n_contacts; i++)
    {
      const FrictionalPointConstraintModel & cmodel = contact_models[i];
      const auto constraint_size = cmodel.size();

      auto lambda_segment = data.lambda_c.segment(row_id, constraint_size);
      typename RigidConstraintData::Matrix6 actInv_transpose1 =
        cmodel.joint1_placement.toActionMatrixInverse();
      actInv_transpose1.transposeInPlace();
      fext[cmodel.joint1_id] += Force(actInv_transpose1.template leftCols<3>() * lambda_segment);
      typename RigidConstraintData::Matrix6 actInv_transpose2 =
        cmodel.joint2_placement.toActionMatrixInverse();
      actInv_transpose2.transposeInPlace();
      fext[cmodel.joint2_id] += Force(actInv_transpose2.template leftCols<3>() * lambda_segment);

      row_id += constraint_size;
    }
    rnea(model, data, q, v, a, fext);
    return data.tau;
  }

  ///
  /// \brief The Contact Inverse Dynamics algorithm. It computes the inverse dynamics in the
  /// presence of contacts, aka the joint torques according to the current state of the system and
  /// the desired joint accelerations.
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
  /// \param[in] R vector representing the diagonal of the compliance matrix.
  /// \param[in] constraint_correction vector representing the constraint correction.
  /// \param[in] settings The settings for the proximal algorithm.
  /// \param[in] lambda_guess initial guess for the contact forces.
  ///
  /// \return The desired joint torques stored in data.tau.
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType,
    typename TangentVectorType1,
    typename TangentVectorType2,
    class ConstraintModelAllocator,
    class ConstraintDataAllocator,
    typename VectorLikeR,
    typename VectorLikeGamma,
    typename VectorLikeLam>
  const typename DataTpl<Scalar, Options, JointCollectionTpl>::TangentVectorType &
  contactInverseDynamics(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<TangentVectorType1> & v,
    const Eigen::MatrixBase<TangentVectorType2> & a,
    Scalar dt,
    const std::vector<
      FrictionalPointConstraintModelTpl<Scalar, Options>,
      ConstraintModelAllocator> & contact_models,
    std::vector<FrictionalPointConstraintDataTpl<Scalar, Options>, ConstraintDataAllocator> &
      contact_datas,
    const Eigen::MatrixBase<VectorLikeR> & R,
    const Eigen::MatrixBase<VectorLikeGamma> & constraint_correction,
    ProximalSettingsTpl<Scalar> & settings,
    const boost::optional<VectorLikeLam> & lambda_guess = boost::none)
  {

    typedef std::reference_wrapper<const FrictionalPointConstraintModelTpl<Scalar, Options>>
      WrappedConstraintModelType;
    typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

    WrappedConstraintModelVector wrapped_constraint_models(
      contact_models.cbegin(), contact_models.cend());

    typedef std::reference_wrapper<FrictionalPointConstraintDataTpl<Scalar, Options>>
      WrappedConstraintDataType;
    typedef std::vector<WrappedConstraintDataType> WrappedConstraintDataVector;

    WrappedConstraintDataVector wrapped_constraint_datas(
      contact_datas.begin(), contact_datas.end());

    return contactInverseDynamics(
      model, data, q, v, a, dt, wrapped_constraint_models, wrapped_constraint_datas, R,
      constraint_correction, settings, lambda_guess);
  }

} // namespace pinocchio

// #if PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION
// #include "pinocchio/algorithm/contact-inverse-dynamics.txx"
// #endif // PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION

#endif // ifndef __pinocchio_algorithm_contact_inverse_dynamics_hpp__
