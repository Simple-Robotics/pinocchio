//
// Copyright (c) 2020-2025 INRIA
//

#include "pinocchio/bindings/python/algorithm/algorithms.hpp"
#include "pinocchio/bindings/python/algorithm/contact-info.hpp"
#include "pinocchio/bindings/python/algorithm/proximal.hpp"
#include "pinocchio/bindings/python/algorithm/contact-cholesky.hpp"

#include "pinocchio/bindings/python/utils/std-vector.hpp"
#include "pinocchio/bindings/python/utils/model-checker.hpp"

#include "pinocchio/algorithm/constrained-dynamics.hpp"

namespace bp = boost::python;

namespace pinocchio
{
  namespace python
  {

    typedef PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(context::RigidConstraintModel)
      RigidConstraintModelVector;
    typedef PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(context::RigidConstraintData)
      RigidConstraintDataVector;

    template<typename ConstraintModel, typename ConstraintData>
    static const context::VectorXs constraintDynamics_proxy(
      const context::Model & model,
      context::Data & data,
      const context::VectorXs & q,
      const context::VectorXs & v,
      const context::VectorXs & tau,
      const std::vector<ConstraintModel, Eigen::aligned_allocator<ConstraintModel>> &
        contact_models,
      std::vector<ConstraintData, Eigen::aligned_allocator<ConstraintData>> & contact_datas,
      context::ProximalSettings & prox_settings)
    {
      return constraintDynamics(
        model, data, q, v, tau, contact_models, contact_datas, prox_settings);
    }

    template<typename ConstraintModel, typename ConstraintData>
    static const context::VectorXs constraintDynamics_proxy_default(
      const context::Model & model,
      context::Data & data,
      const context::VectorXs & q,
      const context::VectorXs & v,
      const context::VectorXs & tau,
      const std::vector<ConstraintModel, Eigen::aligned_allocator<ConstraintModel>> &
        contact_models,
      std::vector<ConstraintData, Eigen::aligned_allocator<ConstraintData>> & contact_datas)
    {
      return constraintDynamics(model, data, q, v, tau, contact_models, contact_datas);
    }

    template<typename ConstraintModel>
    void exposeConstraintDynamicsFor()
    {

      typedef typename ConstraintModel::ConstraintData ConstraintData;
      typedef Eigen::aligned_allocator<ConstraintModel> ConstraintModelAllocator;

      bp::def(
        "initConstraintDynamics",
        &initConstraintDynamics<
          context::Scalar, context::Options, ConstraintModel, JointCollectionDefaultTpl,
          ConstraintModelAllocator>,
        bp::args("model", "data", "contact_models"),
        "This function allows to allocate the memory before hand for contact dynamics algorithms.\n"
        "This allows to avoid online memory allocation when running these algorithms.",
        mimic_not_supported_function<>(0));

      bp::def(
        "constraintDynamics", constraintDynamics_proxy<ConstraintModel, ConstraintData>,
        bp::args(
          "model", "data", "q", "v", "tau", "contact_models", "contact_datas", "prox_settings"),
        "Computes the forward dynamics with contact constraints according to a given list of "
        "Contact information.\n"
        "When using constraintDynamics for the first time, you should call first "
        "initConstraintDynamics to initialize the internal memory used in the algorithm.\n"
        "This function returns joint acceleration of the system. The contact forces are "
        "stored in the list data.contact_forces.",
        mimic_not_supported_function<>(0));

      bp::def(
        "constraintDynamics", constraintDynamics_proxy_default<ConstraintModel, ConstraintData>,
        bp::args("model", "data", "q", "v", "tau", "contact_models", "contact_datas"),
        "Computes the forward dynamics with contact constraints according to a given list of "
        "Contact information.\n"
        "When using constraintDynamics for the first time, you should call first "
        "initConstraintDynamics to initialize the internal memory used in the algorithm.\n"
        "This function returns joint acceleration of the system. The contact forces are stored in "
        "the list data.contact_forces.",
        mimic_not_supported_function<>(0));
    }

    void exposeConstraintDynamics()
    {
      using namespace Eigen;

      // Expose type of contacts
      if (!register_symbolic_link_to_registered_type<ContactType>())
      {
        bp::enum_<ContactType>("ContactType")
          .value("CONTACT_3D", CONTACT_3D)
          .value("CONTACT_6D", CONTACT_6D)
          .value("CONTACT_UNDEFINED", CONTACT_UNDEFINED);
      }

      ProximalSettingsPythonVisitor<context::ProximalSettings>::expose();

      RigidConstraintModelPythonVisitor<context::RigidConstraintModel>::expose();
      RigidConstraintDataPythonVisitor<context::RigidConstraintData>::expose();

      StdVectorPythonVisitor<RigidConstraintModelVector>::expose("StdVec_RigidConstraintModel");

      StdVectorPythonVisitor<RigidConstraintDataVector>::expose("StdVec_RigidConstraintData");

      ContactCholeskyDecompositionPythonVisitor<context::ContactCholeskyDecomposition>::expose();

      exposeConstraintDynamicsFor<context::RigidConstraintModel>();
      // exposeConstraintDynamicsFor<WeldConstraintModel>();
    }
  } // namespace python
} // namespace pinocchio
