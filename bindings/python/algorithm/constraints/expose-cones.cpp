//
// Copyright (c) 2022-2024 INRIA
//

#include "pinocchio/serialization/aligned-vector.hpp"

#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/set-coulomb-friction-cone.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/set-box-set.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/set-joint-limit-cone.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/set-trivial-cones.hpp"

// #include "pinocchio/bindings/python/serialization/serialization.hpp"
#include "pinocchio/bindings/python/utils/std-aligned-vector.hpp"

#ifdef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS
namespace eigenpy
{
  // has_operator_equal return true for CoulombFrictionCone and DualCoulombFrictionCone with casadi
  // but it should not.
  // We provide specialization to enforce good result.
  template<>
  struct has_operator_equal<
    pinocchio::python::context::CoulombFrictionCone,
    pinocchio::python::context::CoulombFrictionCone>
  {
    typedef std::false_type type;
  };
  template<>
  struct has_operator_equal<
    pinocchio::python::context::DualCoulombFrictionCone,
    pinocchio::python::context::DualCoulombFrictionCone>
  {
    typedef std::false_type type;
  };
} // namespace eigenpy
#endif // ifdef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS

namespace pinocchio
{
  namespace python
  {
    void exposeCones()
    {
      CoulombFrictionConePythonVisitor<context::CoulombFrictionCone>::expose();
      StdVectorPythonVisitor<context::CoulombFrictionConeVector>::expose(
        "StdVec_CoulombFrictionCone");
      // #ifndef PINOCCHIO_PYTHON_NO_SERIALIZATION
      //       serialize<StdAlignedVectorPythonVisitor<context::CoulombFrictionCone>::vector_type>();
      // #endif

      DualCoulombFrictionConePythonVisitor<context::DualCoulombFrictionCone>::expose();
      StdVectorPythonVisitor<context::DualCoulombFrictionConeVector>::expose(
        "StdVec_DualCoulombFrictionCone");
      // #ifndef PINOCCHIO_PYTHON_NO_SERIALIZATION
      //       serialize<StdAlignedVectorPythonVisitor<context::DualCoulombFrictionCone>::vector_type>();
      // #endif

      BoxSetPythonVisitor<context::BoxSet>::expose();
      JointLimitConstraintConePythonVisitor<context::JointLimitConstraintCone>::expose();
      TrivialConePythonVisitor<context::NullSet>::expose(
        "NullSet", "Set reduce to 0 singleton in R^d.");
      TrivialConePythonVisitor<context::UnboundedSet>::expose("UnboundedSet", "Set R^d.");
      TrivialConePythonVisitor<context::PositiveOrthantCone>::expose(
        "PositiveOrthantCone", "Set R_+^d.");
      TrivialConePythonVisitor<context::NegativeOrthantCone>::expose(
        "NegativeOrthantCone", "Set R_-^d.");
    }
  } // namespace python
} // namespace pinocchio
