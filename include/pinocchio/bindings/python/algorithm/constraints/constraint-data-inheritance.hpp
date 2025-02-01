//
// Copyright (c) 2025 INRIA
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
    template<class T, class TBase>
    struct ConstraintDataInheritancePythonVisitor
    : public bp::def_visitor<ConstraintDataInheritancePythonVisitor<T, TBase>>
    {
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl;
      }
    };

    // Specialize
    template<class T>
    struct ConstraintDataInheritancePythonVisitor<T, FrameConstraintDataBase<T>>
    : public bp::def_visitor<ConstraintDataInheritancePythonVisitor<T, FrameConstraintDataBase<T>>>
    {
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def(bp::init<>(bp::arg("self"), "Default constructor"))
          .def(bp::init<const typename T::ConstraintModel &>(bp::args("self", "constraint_model")))
          .PINOCCHIO_ADD_PROPERTY(T, constraint_force, "Resulting force.")
          .PINOCCHIO_ADD_PROPERTY(T, oMc1, "Placement of the constraint frame 1 wrt WORLD.")
          .PINOCCHIO_ADD_PROPERTY(T, oMc2, "Placement of the constraint frame 2 wrt WORLD.")
          .PINOCCHIO_ADD_PROPERTY(T, c1Mc2, "Placement of the constraint frame 2 wrt frame 1.")
          .PINOCCHIO_ADD_PROPERTY(T, constraint_position_error, "Constraint placement (6D) error.")
          .PINOCCHIO_ADD_PROPERTY(T, constraint_velocity_error, "Constraint velocity (6D) error.")
          .PINOCCHIO_ADD_PROPERTY(
            T, constraint_acceleration_error, "Constraint acceleration (6D) error.")
          .PINOCCHIO_ADD_PROPERTY(
            T, constraint_acceleration_biais_term, "Constraint acceleration (6D) term.");
      }
    };

    template<class T>
    struct ConstraintDataInheritancePythonVisitor<T, PointConstraintDataBase<T>>
    : public bp::def_visitor<ConstraintDataInheritancePythonVisitor<T, PointConstraintDataBase<T>>>
    {
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def(bp::init<>(bp::arg("self"), "Default constructor"))
          .def(bp::init<const typename T::ConstraintModel &>(bp::args("self", "constraint_model")))
          .PINOCCHIO_ADD_PROPERTY(T, constraint_force, "Resulting force.")
          .PINOCCHIO_ADD_PROPERTY(T, oMc1, "Placement of the constraint frame 1 wrt WORLD.")
          .PINOCCHIO_ADD_PROPERTY(T, oMc2, "Placement of the constraint frame 2 wrt WORLD.")
          .PINOCCHIO_ADD_PROPERTY(T, c1Mc2, "Placement of the constraint frame 2 wrt frame 1.")
          .PINOCCHIO_ADD_PROPERTY(T, constraint_position_error, "Constraint position (3D) error.")
          .PINOCCHIO_ADD_PROPERTY(T, constraint_velocity_error, "Constraint velocity (3D) error.")
          .PINOCCHIO_ADD_PROPERTY(
            T, constraint_acceleration_error, "Constraint acceleration (3D) error.")
          .PINOCCHIO_ADD_PROPERTY(
            T, constraint_acceleration_biais_term, "Constraint acceleration (3D) term.");
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_data_inheritance_hpp__
