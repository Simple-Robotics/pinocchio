//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_contact_inverse_dynamics_hpp__
#define __pinocchio_algorithm_contact_inverse_dynamics_hpp__

#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/contact-jacobian.hpp"
#include "pinocchio/algorithm/proximal.hpp"

namespace pinocchio
{

  ///
  /// \brief Compute the contact forces given a target velocity of contact points.
  ///
  /// \param[in] contact_models The vector of constraint models.
  /// \param[in] c_ref The desired constraint velocity.
  /// \param[in] R vector representing the diagonal of the compliance matrix.
  /// \param[in,out] lambda Vector of solution. Should be initialized with zeros or from an initial
  /// estimate.
  /// \param[in,out] settings The settings for the proximal algorithm
  /// \param[in] solve_ncp whether to solve the NCP (true) or CCP (false).
  ///
  template<
    typename Scalar,
    int Options,
    template<typename T> class Holder,
    class ConstraintModelAllocator,
    typename VectorLikeC,
    typename VectorLikeR,
    typename VectorLikeResult>
  bool computeInverseDynamicsConstraintForces(
    const std::vector<
      Holder<const FrictionalPointConstraintModelTpl<Scalar, Options>>,
      ConstraintModelAllocator> & contact_models,
    const Eigen::MatrixBase<VectorLikeC> & c_ref,
    const Eigen::MatrixBase<VectorLikeR> & R,
    const Eigen::MatrixBase<VectorLikeResult> & _lambda,
    ProximalSettingsTpl<Scalar> & settings,
    bool solve_ncp = true)
  {
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;
    typedef FrictionalPointConstraintModelTpl<Scalar, Options> ConstraintModel;

    auto & lambda = _lambda.const_cast_derived();

    // PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_correction.size(), problem_size);
    const std::size_t n_constraints = contact_models.size();
    PINOCCHIO_CHECK_ARGUMENT_SIZE(contact_models.size(), n_constraints);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(lambda.size(), R.size());
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      check_expression_if_real<Scalar>(settings.mu >= Scalar(0)), "mu has to be strictly positive");

    const Eigen::Index problem_size = R.size();
    const VectorXs R_prox = R + VectorXs::Constant(problem_size, settings.mu);

    assert(
      R.size() > 0 && check_expression_if_real<Scalar>(R_prox.minCoeff() >= Scalar(0))
      && "The minimal value of R_prox should strictly positive");

    Scalar lambda_c_prev_norm_inf = lambda.template lpNorm<Eigen::Infinity>();

    bool has_converged = false;
    settings.iter = 1;
    for (; settings.iter <= settings.max_iter; ++settings.iter)
    {
      bool abs_prec_reached = false, rel_prec_reached = false;
      settings.relative_residual = settings.absolute_residual = Scalar(0);

      Eigen::DenseIndex row_id = 0;
      for (std::size_t constraint_id = 0; constraint_id < n_constraints; ++constraint_id)
      {
        const ConstraintModel & cmodel = contact_models[constraint_id];
        const auto constraint_size = cmodel.size();

        const auto & cone = cmodel.set();
        auto lambda_segment = lambda.segment(row_id, constraint_size);
        const Vector3 lambda_c_previous = lambda_segment;

        const auto R_segment = R.segment(row_id, constraint_size);
        const auto R_prox_segment = R_prox.segment(row_id, constraint_size);
        const auto c_ref_segment = c_ref.segment(row_id, constraint_size);

        const Vector3 sigma_segment =
          c_ref_segment + (R_segment.array() * lambda_segment.array()).matrix();
        Vector3 desaxce_correction = Vector3::Zero();
        if (solve_ncp)
          desaxce_correction = cone.computeNormalCorrection(sigma_segment);
        const Vector3 c_cor_segment = c_ref_segment + desaxce_correction;

        // Update segment value
        const Vector3 lambda_ref =
          -(Vector3(c_cor_segment - settings.mu * lambda_c_previous).array()
            / R_prox_segment.array());
        lambda_segment = cone.weightedProject(lambda_ref, R_prox_segment);

        // Compute convergence criteria
        const Scalar contact_complementarity = cone.computeConicComplementarity(
          Vector3(sigma_segment + desaxce_correction), lambda_segment);
        const Scalar dual_feasibility =
          std::abs(math::min(0., sigma_segment(2))); // proxy of dual feasibility
        settings.absolute_residual = math::max(
          settings.absolute_residual, math::max(contact_complementarity, dual_feasibility));

        const Vector3 dlambda_c = lambda_segment - lambda_c_previous;
        const Scalar proximal_metric = dlambda_c.template lpNorm<Eigen::Infinity>();
        settings.relative_residual = math::max(settings.relative_residual, proximal_metric);

        row_id += constraint_size;
      }

      const Scalar lambda_c_norm_inf = lambda.template lpNorm<Eigen::Infinity>();

      if (check_expression_if_real<Scalar, false>(
            settings.absolute_residual <= settings.absolute_accuracy))
        abs_prec_reached = true;

      if (check_expression_if_real<Scalar, false>(
            settings.relative_residual
            <= settings.relative_accuracy * math::max(lambda_c_norm_inf, lambda_c_prev_norm_inf)))
        rel_prec_reached = true;

      if (abs_prec_reached || rel_prec_reached)
      {
        has_converged = true;
        break;
      }

      lambda_c_prev_norm_inf = lambda_c_norm_inf;
    }

    return has_converged;
  }

  ///
  /// \brief Compute the contact forces given a target velocity of contact points.
  ///
  /// \param[in] contact_models The vector of constraint models.
  /// \param[in] c_ref The desired constraint velocity.
  /// \param[in] R vector representing the diagonal of the compliance matrix.
  /// \param[in,out] lambda_sol Vector of solution. Should be initialized with zeros or from an
  /// initial estimate
  /// \param[in,out] settings The settings for the proximal algorithm
  /// \param[in] solve_ncp whether to solve the NCP (true) or CCP (false).
  ///
  template<
    typename Scalar,
    int Options,
    class ConstraintModelAllocator,
    typename VectorLikeC,
    typename VectorLikeR,
    typename VectorLikeResult>
  bool computeInverseDynamicsConstraintForces(
    const std::vector<
      FrictionalPointConstraintModelTpl<Scalar, Options>,
      ConstraintModelAllocator> & contact_models,
    const Eigen::MatrixBase<VectorLikeC> & c_ref,
    const Eigen::MatrixBase<VectorLikeR> & R,
    const Eigen::MatrixBase<VectorLikeResult> & lambda_sol,
    ProximalSettingsTpl<Scalar> & settings,
    bool solve_ncp = true)
  {
    typedef FrictionalPointConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
    typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

    WrappedConstraintModelVector wrapped_constraint_models(
      contact_models.cbegin(), contact_models.cend());

    return computeInverseDynamicsConstraintForces(
      wrapped_constraint_models, c_ref.derived(), R.derived(), lambda_sol.const_cast_derived(),
      settings, solve_ncp);
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
  /// \param[in] lambda_sol initial guess for the contact forces
  /// \param[in] settings The settings for the proximal algorithm
  /// \param[in] solve_ncp whether to solve the NCP (true) or CCP (false).
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
    typename VectorLikeGamma,
    typename VectorLikeR,
    typename VectorLikeLam>
  const typename DataTpl<Scalar, Options, JointCollectionTpl>::TangentVectorType &
  contactInverseDynamics(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<TangentVectorType1> & v,
    const Eigen::MatrixBase<TangentVectorType2> & a,
    const Scalar dt,
    const std::vector<
      Holder<const FrictionalPointConstraintModelTpl<Scalar, Options>>,
      ConstraintModelAllocator> & contact_models,
    std::vector<
      Holder<FrictionalPointConstraintDataTpl<Scalar, Options>>,
      ConstraintDataAllocator> & contact_datas,
    const Eigen::MatrixBase<VectorLikeGamma> & constraint_correction,
    const Eigen::MatrixBase<VectorLikeR> & R,
    const Eigen::MatrixBase<VectorLikeLam> & _lambda_sol,
    ProximalSettingsTpl<Scalar> & settings,
    bool solve_ncp = true)
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef FrictionalPointConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef typename ConstraintModel::ConstraintData ConstraintData;

    typedef typename Model::MatrixXs MatrixXs;
    typedef typename Model::VectorXs VectorXs;

    auto & lambda_sol = _lambda_sol.const_cast_derived();

    const Eigen::Index problem_size = R.size();
    const std::size_t n_constraints = contact_models.size();

    MatrixXs J = MatrixXs::Zero(problem_size, model.nv); // TODO: malloc
    getConstraintsJacobian(model, data, contact_models, contact_datas, J);
    VectorXs v_ref, c_ref, tau_c;
    v_ref.noalias() = v + dt * a;
    c_ref.noalias() = J * v_ref; // TODO should rather use the displacement
    c_ref += constraint_correction;
    c_ref /= dt; // we work with a formulation on forces
    computeInverseDynamicsConstraintForces(
      contact_models, c_ref, R, lambda_sol, settings, solve_ncp);

    {
      rnea(model, data, q, v, a);
      auto & tau = data.tau;
      Eigen::DenseIndex row_id = 0;
      for (std::size_t i = 0; i < n_constraints; i++)
      {
        const ConstraintModel & cmodel = contact_models[i];
        ConstraintData & cdata = contact_datas[i];
        const auto constraint_size = cmodel.size();

        const auto lambda_segment = lambda_sol.segment(row_id, constraint_size);
        cmodel.jacobianTransposeMatrixProduct(model, data, cdata, lambda_segment, tau, RmTo());

        row_id += constraint_size;
      }
    }
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
  /// \param[in] lambda_sol initial guess for the contact forces
  /// \param[in] settings The settings for the proximal algorithm
  /// \param[in] solve_ncp whether to solve the NCP (true) or CCP (false).
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
    typename VectorLikeGamma,
    typename VectorLikeR,
    typename VectorLikeLam>
  const typename DataTpl<Scalar, Options, JointCollectionTpl>::TangentVectorType &
  contactInverseDynamics(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<TangentVectorType1> & v,
    const Eigen::MatrixBase<TangentVectorType2> & a,
    const Scalar dt,
    const std::vector<
      FrictionalPointConstraintModelTpl<Scalar, Options>,
      ConstraintModelAllocator> & contact_models,
    std::vector<FrictionalPointConstraintDataTpl<Scalar, Options>, ConstraintDataAllocator> &
      contact_datas,
    const Eigen::MatrixBase<VectorLikeGamma> & constraint_correction,
    const Eigen::MatrixBase<VectorLikeR> & R,
    const Eigen::MatrixBase<VectorLikeLam> & lambda_sol,
    ProximalSettingsTpl<Scalar> & settings,
    bool solve_ncp = true)
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
      model, data, q, v, a, dt, wrapped_constraint_models, wrapped_constraint_datas,
      constraint_correction, R, lambda_sol.const_cast_derived(), settings, solve_ncp);
  }

} // namespace pinocchio

// #if PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION
// #include "pinocchio/algorithm/contact-inverse-dynamics.txx"
// #endif // PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION

#endif // ifndef __pinocchio_algorithm_contact_inverse_dynamics_hpp__
