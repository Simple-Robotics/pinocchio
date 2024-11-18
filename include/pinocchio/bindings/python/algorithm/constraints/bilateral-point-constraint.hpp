//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_python_algorithm_bilateral_box_constraint_hpp__
#define __pinocchio_python_algorithm_bilateral_box_constraint_hpp__

#include <eigenpy/eigenpy.hpp>

#include "pinocchio/algorithm/constraints/bilateral-point-constraint.hpp"

#include "pinocchio/bindings/python/utils/cast.hpp"
#include "pinocchio/bindings/python/utils/copyable.hpp"
#include "pinocchio/bindings/python/utils/macros.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<typename BilateralPointConstraintModel>
    struct BilateralPointConstraintModelPythonVisitor
    : public boost::python::def_visitor<
        BilateralPointConstraintModelPythonVisitor<BilateralPointConstraintModel>>
    {
      typedef BilateralPointConstraintModel Self;
      typedef typename BilateralPointConstraintModel::Scalar Scalar;

      typedef ModelTpl<Scalar, BilateralPointConstraintModel::Options, JointCollectionDefaultTpl>
        Model;

      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl

          .def(bp::init<const Model &, JointIndex, const SE3 &, JointIndex, const SE3 &>(
            (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id"), bp::arg("joint1_placement"),
             bp::arg("joint2_id"), bp::arg("joint2_placement")),
            "Contructor from given joint index and placement for the two joints "
            "implied in the constraint."))

          .def(bp::init<const Model &, JointIndex, const SE3 &>(
            (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id"), bp::arg("joint1_placement")),
            "Contructor from given joint index and placement for the two joints "
            "implied in the constraint."))

          .PINOCCHIO_ADD_PROPERTY(Self, joint1_id, "Index of first parent joint in the model tree.")
          .PINOCCHIO_ADD_PROPERTY(
            Self, joint2_id, "Index of second parent joint in the model tree.")
          .PINOCCHIO_ADD_PROPERTY(
            Self, joint1_placement, "Relative placement with respect to the frame of joint1.")
          .PINOCCHIO_ADD_PROPERTY(
            Self, joint2_placement, "Relative placement with respect to the frame of joint2.")

#ifndef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS
          .def(bp::self == bp::self)
          .def(bp::self != bp::self)
#endif
          ;
      }

      static void expose()
      {
        bp::class_<BilateralPointConstraintModel>(
          "BilateralPointConstraintModel",
          "Bilateral point constraint. Only the translation of the placements is constrained.\n",
          bp::no_init)
          .def(BilateralPointConstraintModelPythonVisitor())
          .def(CopyableVisitor<BilateralPointConstraintModel>());
      }
    };

  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_bilateral_box_constraint_hpp__
