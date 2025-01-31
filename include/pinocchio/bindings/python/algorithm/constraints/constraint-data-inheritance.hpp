//
// Copyright (c) 2015-2025 CNRS INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_data_inheritance_hpp__
#define __pinocchio_python_algorithm_constraints_data_inheritance_hpp__

#include <boost/python.hpp>
#include <eigenpy/exception.hpp>
#include <eigenpy/eigen-to-python.hpp>

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/frame-constraint-data-base.hpp"
#include "pinocchio/algorithm/constraints/point-constraint-data-base.hpp"
#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/bindings/python/utils/macros.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    // Default inheritance Visitor Template
    template<class ConstraintDataDerived, class ConstraintDataBase>
    struct ConstraintDataInheritancePythonVisitor
    : public bp::def_visitor<ConstraintDataInheritancePythonVisitor<ConstraintDataDerived, ConstraintDataBase>>
    {
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl
          ;
      }
    };

    // Specialize
    template<class ConstraintDataDerived>
    struct ConstraintDataInheritancePythonVisitor<ConstraintDataDerived, FrameConstraintDataBase<ConstraintDataDerived>>
    : public bp::def_visitor<
      ConstraintDataInheritancePythonVisitor<ConstraintDataDerived, FrameConstraintDataBase<ConstraintDataDerived>>
    >
    {
      typedef ConstraintDataDerived Self;
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl
          .def(bp::init<>(bp::arg("self"), "Default constructor"))
          .def(
            bp::init<const typename Self::ConstraintModel &>(bp::args("self", "constraint_model")))
          .PINOCCHIO_ADD_PROPERTY(Self, constraint_force, "Resulting force.")
          .PINOCCHIO_ADD_PROPERTY(Self, oMc1, "Placement of the constraint frame 1 wrt WORLD.")
          .PINOCCHIO_ADD_PROPERTY(Self, oMc2, "Placement of the constraint frame 2 wrt WORLD.")
          .PINOCCHIO_ADD_PROPERTY(Self, c1Mc2, "Placement of the constraint frame 2 wrt frame 1.")
          .PINOCCHIO_ADD_PROPERTY(Self, constraint_position_error, "Constraint placement (6D) error.")
          .PINOCCHIO_ADD_PROPERTY(Self, constraint_velocity_error, "Constraint velocity (6D) error.")
          .PINOCCHIO_ADD_PROPERTY(Self, constraint_acceleration_error, "Constraint acceleration (6D) error.")
          .PINOCCHIO_ADD_PROPERTY(Self, constraint_acceleration_biais_term, "Constraint acceleration (6D) term.")
          ;
      }
    };

    template<class ConstraintDataDerived>
    struct ConstraintDataInheritancePythonVisitor<ConstraintDataDerived, PointConstraintDataBase<ConstraintDataDerived>>
    : public bp::def_visitor<
      ConstraintDataInheritancePythonVisitor<ConstraintDataDerived, PointConstraintDataBase<ConstraintDataDerived>>
    >
    {
      typedef ConstraintDataDerived Self;
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl
          .def(bp::init<>(bp::arg("self"), "Default constructor"))
          .def(
            bp::init<const typename Self::ConstraintModel &>(bp::args("self", "constraint_model")))
          .PINOCCHIO_ADD_PROPERTY(Self, constraint_force, "Resulting force.")
          .PINOCCHIO_ADD_PROPERTY(Self, oMc1, "Placement of the constraint frame 1 wrt WORLD.")
          .PINOCCHIO_ADD_PROPERTY(Self, oMc2, "Placement of the constraint frame 2 wrt WORLD.")
          .PINOCCHIO_ADD_PROPERTY(Self, c1Mc2, "Placement of the constraint frame 2 wrt frame 1.")
          .PINOCCHIO_ADD_PROPERTY(Self, constraint_position_error, "Constraint position (3D) error.")
          .PINOCCHIO_ADD_PROPERTY(Self, constraint_velocity_error, "Constraint velocity (3D) error.")
          .PINOCCHIO_ADD_PROPERTY(Self, constraint_acceleration_error, "Constraint acceleration (3D) error.")
          .PINOCCHIO_ADD_PROPERTY(Self, constraint_acceleration_biais_term, "Constraint acceleration (3D) term.")
          ;
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_data_inheritance_hpp__
