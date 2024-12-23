//
// Copyright (c) 2024 INRIA
//

#define BOOST_PYTHON_MAX_ARITY 24

#include "pinocchio/bindings/python/algorithm/algorithms.hpp"
#include "pinocchio/bindings/python/utils/std-vector.hpp"
#include "pinocchio/algorithm/contact-inverse-dynamics.hpp"

#include <boost/optional.hpp>

namespace pinocchio
{
  namespace python
  {

#ifndef PINOCCHIO_PYTHON_SKIP_ALGORITHM_CONSTRAINED_DYNAMICS
    typedef context::Scalar Scalar;
    typedef context::VectorXs VectorXs;
    typedef const Eigen::Ref<const VectorXs> ConstRefVectorXs;

    static bp::tuple computeInverseDynamicsConstraintForces_wrapper(
      const VectorXs & c_ref,
      const context::FrictionalPointConstraintModelVector & contact_models,
      const VectorXs & R,
      const boost::optional<VectorXs> & lambda_guess,
      ProximalSettingsTpl<Scalar> & settings,
      bool solve_ncp)
    {
      VectorXs lambda_sol;
      if (lambda_guess)
        lambda_sol = lambda_guess.get();
      else
        lambda_sol = VectorXs::Zero(R.size());

      const bool has_converged = computeInverseDynamicsConstraintForces(
        contact_models, c_ref, R, lambda_sol, settings, solve_ncp);
      return bp::make_tuple(has_converged, bp::object(lambda_sol));
    }

//    static ConstRefVectorXs contactInverseDynamics_wrapper(
//      const context::Model & model,
//                                                           context::Data & data,
//      ConstRefVectorXs & q,
//      ConstRefVectorXs & v,
//      ConstRefVectorXs & a,
//      Scalar dt,
//      const context::FrictionalPointConstraintModelVector & contact_models,
//      context::FrictionalPointConstraintDataVector & contact_datas,
//      ConstRefVectorXs & R,
//      ConstRefVectorXs & constraint_correction,
//      ProximalSettingsTpl<Scalar> & settings,
//      const boost::optional<ConstRefVectorXs> & lambda_guess = boost::none)
//    {
//      return contactInverseDynamics(
//        model, data, q, v, a, dt, contact_models, contact_datas, R, constraint_correction,
//        settings, lambda_guess);
//    }
#endif // PINOCCHIO_PYTHON_SKIP_ALGORITHM_CONSTRAINED_DYNAMICS

    void exposeContactInverseDynamics()
    {
#ifndef PINOCCHIO_PYTHON_SKIP_ALGORITHM_CONSTRAINED_DYNAMICS
      bp::def(
        "computeInverseDynamicsConstraintForces", computeInverseDynamicsConstraintForces_wrapper,
        (bp::args("contact_models", "c_ref", "R"), bp::arg("lambda_guess") = boost::none,
         bp::arg("settings"), bp::arg("solve_ncp") = true),
        "Computes the inverse dynamics with frictional contacts. Returns a tuple containing "
        "(has_converged, lambda_sol).\n\n"
        "Parameters:\n"
        "\tcontact_models: list of contact models\n"
        "\tc_ref: the reference velocity of contact points\n"
        "\tR: vector representing the diagonal of the compliance matrix\n"
        "\tlambda_guess: optional initial guess for contact forces\n"
        "\tsettings: the settings of the proximal algorithm\n"
        "\tsolve_ncp: whether to solve the NCP (true) or CCP (false)\n");

//      bp::def(
//        "contactInverseDynamics", contactInverseDynamics_wrapper,
//        (bp::arg("model"), "data", "q", "v", "a", "dt", "contact_models", "contact_datas", "R",
//         "constraint_correction", bp::arg("settings"), bp::arg("lambda_guess") = boost::none),
//        "Compute the inverse dynamics with frictional contacts, store the result in Data and "
//        "return it.\n\n"
//        "Parameters:\n"
//        "\tmodel: model of the kinematic tree\n"
//        "\tdata: data related to the model\n"
//        "\tq: the joint configuration vector (size model.nq)\n"
//        "\tv: the joint velocity vector (size model.nv)\n"
//        "\ta: the joint acceleration vector (size model.nv)\n"
//        "\tdt: the time step\n"
//        "\tcontact_models: list of contact models\n"
//        "\tcontact_datas: list of contact datas\n"
//        "\tR: vector representing the diagonal of the compliance matrix\n"
//        "\tconstraint_correction: vector representing the constraint correction\n"
//        "\tsettings: the settings of the proximal algorithm\n"
//        "\tlambda_guess: initial guess for contact forces\n");
#endif // PINOCCHIO_PYTHON_SKIP_ALGORITHM_CONSTRAINED_DYNAMICS
    }
  } // namespace python
} // namespace pinocchio
