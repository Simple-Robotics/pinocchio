//
// Copyright (c) 2021-2025 INRIA
//

#include "pinocchio/bindings/python/algorithm/algorithms.hpp"
#include "pinocchio/algorithm/contact-jacobian.hpp"

#include "pinocchio/bindings/python/utils/std-vector.hpp"
#include "pinocchio/bindings/python/utils/model-checker.hpp"

namespace bp = boost::python;

namespace pinocchio
{
  namespace python
  {

    template<typename ConstraintModel, typename ConstraintData>
    static context::MatrixXs getConstraintJacobian_proxy(
      const context::Model & model,
      const context::Data & data,
      const ConstraintModel & constraint_model,
      ConstraintData & constraint_data)
    {
      context::MatrixXs J(constraint_model.size(), model.nv);
      J.setZero();
      getConstraintJacobian(model, data, constraint_model, constraint_data, J);
      return J;
    }

    template<typename ConstraintModelVector, typename ConstraintDataVector>
    static context::MatrixXs getConstraintsJacobian_proxy(
      const context::Model & model,
      const context::Data & data,
      const ConstraintModelVector & constraint_models,
      ConstraintDataVector & constraint_datas)
    {
      const Eigen::DenseIndex constraint_size = getTotalConstraintSize(constraint_models);
      context::MatrixXs J(constraint_size, model.nv);
      J.setZero();
      getConstraintsJacobian(model, data, constraint_models, constraint_datas, J);
      return J;
    }

    void exposeContactJacobian()
    {
      bp::def(
        "getConstraintJacobian",
        getConstraintJacobian_proxy<context::RigidConstraintModel, context::RigidConstraintData>,
        bp::args("model", "data", "constraint_model", "constraint_data"),
        "Computes the kinematic Jacobian associatied with a given constraint model.",
        mimic_not_supported_function<>(0));
      bp::def(
        "getConstraintsJacobian",
        getConstraintsJacobian_proxy<
          context::RigidConstraintModelVector, context::RigidConstraintDataVector>,
        bp::args("model", "data", "constraint_models", "constraint_datas"),
        "Computes the kinematic Jacobian associatied with a given set of constraint models.",
        mimic_not_supported_function<>(0));

      bp::def(
        "getConstraintJacobian",
        getConstraintJacobian_proxy<context::ConstraintModel, context::ConstraintData>,
        bp::args("model", "data", "constraint_model", "constraint_data"),
        "Computes the kinematic Jacobian associatied with a given constraint model.",
        mimic_not_supported_function<>(0));
      bp::def(
        "getConstraintsJacobian",
        getConstraintsJacobian_proxy<context::ConstraintModelVector, context::ConstraintDataVector>,
        bp::args("model", "data", "constraint_models", "constraint_datas"),
        "Computes the kinematic Jacobian associatied with a given set of constraint models.",
        mimic_not_supported_function<>(0));
    }
  } // namespace python
} // namespace pinocchio
