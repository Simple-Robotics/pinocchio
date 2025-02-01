//
// Copyright (c) 2015-2025 CNRS INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_model_inheritance_hpp__
#define __pinocchio_python_algorithm_constraints_model_inheritance_hpp__

#include <boost/python.hpp>
#include <eigenpy/exception.hpp>
#include <eigenpy/eigen-to-python.hpp>

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/frame-constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/point-constraint-model-base.hpp"
#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/bindings/python/utils/macros.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    // Default inheritance Visitor Template
    template<class T, class TBase>
    struct ConstraintModelInheritancePythonVisitor
    : public bp::def_visitor<ConstraintModelInheritancePythonVisitor<T, TBase>>
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
    template<class T>
    struct ConstraintModelInheritancePythonVisitor<T, FrameConstraintModelBase<T>>
    : public bp::def_visitor<ConstraintModelInheritancePythonVisitor<T, FrameConstraintModelBase<T>>>
    {
      typedef typename T::Scalar Scalar;
      typedef ModelTpl<Scalar, T::Options, JointCollectionDefaultTpl> Model;
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl
          .def(bp::init<const Model &, JointIndex, const SE3 &, JointIndex, const SE3 &>(
            (bp::arg("self"), bp::arg("model"),
             bp::arg("joint1_id"), bp::arg("joint1_placement"),
             bp::arg("joint2_id"), bp::arg("joint2_placement")),
            "Contructor from given joint index and placement for the two joints "
            "implied in the constraint."))
          .def(bp::init<const Model &, JointIndex, const SE3 &>(
            (bp::arg("self"), bp::arg("model"),
             bp::arg("joint1_id"), bp::arg("joint1_placement")),
            "Contructor from given joint index and placement of the first joint "
            "implied in the constraint."))
          .def(bp::init<const Model &, JointIndex, const JointIndex &>(
            (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id"), bp::arg("joint2_id")),
            "Contructor from given joint index for the two joints "
            "implied in the constraint."))
          .def(bp::init<const Model &, JointIndex&>(
            (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id")),
            "Contructor from given joint index of the first joint "
            "implied in the constraint."))
          .PINOCCHIO_ADD_PROPERTY(T, joint1_id, "Index of the first joint in the model tree")
          .PINOCCHIO_ADD_PROPERTY(T, joint2_id, "Index of the second joint in the model tree")
          .PINOCCHIO_ADD_PROPERTY(T, joint1_placement, "Position of attached point with respect to the frame of joint1.")
          .PINOCCHIO_ADD_PROPERTY(T, joint2_placement, "Position of attached point with respect to the frame of joint2.")
          .PINOCCHIO_ADD_PROPERTY(T, desired_constraint_offset, "Desired constraint shift at position level.")
          .PINOCCHIO_ADD_PROPERTY(T, desired_constraint_velocity, "Desired constraint velocity at velocity level.")
          .PINOCCHIO_ADD_PROPERTY(T, desired_constraint_acceleration, "Desired constraint velocity at acceleration level.")
          .PINOCCHIO_ADD_PROPERTY(T, corrector_parameters, "Corrector parameters.")
          .PINOCCHIO_ADD_PROPERTY(T, colwise_joint1_sparsity, "Colwise sparsity pattern associated with joint 1.")
          .PINOCCHIO_ADD_PROPERTY(T, colwise_joint2_sparsity, "Colwise sparsity pattern associated with joint 2.")
          .PINOCCHIO_ADD_PROPERTY(T, joint1_span_indexes, " Jointwise span indexes associated with joint 1.")
          .PINOCCHIO_ADD_PROPERTY(T, joint2_span_indexes, "Jointwise span indexes associated with joint 2.")
          .PINOCCHIO_ADD_PROPERTY(T, loop_span_indexes, "Loop span indexes.")
          .PINOCCHIO_ADD_PROPERTY(T, colwise_sparsity, "Sparsity pattern associated to the constraint.")
          .PINOCCHIO_ADD_PROPERTY(T, colwise_span_indexes, "Indexes of the columns spanned by the constraints.")
          // .def("getRowSparsityPattern", ...)
          // .def("getRowActiveIndexes", ...)
          // .def("getA1", ...)
          // .def("getA2", ...)
          ;
      }
    };

    template<class T>
    struct ConstraintModelInheritancePythonVisitor<T, PointConstraintModelBase<T>>
    : public bp::def_visitor<ConstraintModelInheritancePythonVisitor<T, PointConstraintModelBase<T>>>
    {
      typedef typename T::Scalar Scalar;
      typedef ModelTpl<Scalar, T::Options, JointCollectionDefaultTpl> Model;
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl
          .def(bp::init<const Model &, JointIndex, const SE3 &, JointIndex, const SE3 &>(
            (bp::arg("self"), bp::arg("model"),
             bp::arg("joint1_id"), bp::arg("joint1_placement"),
             bp::arg("joint2_id"), bp::arg("joint2_placement")),
            "Contructor from given joint index and placement for the two joints "
            "implied in the constraint."))
          .def(bp::init<const Model &, JointIndex, const SE3 &>(
            (bp::arg("self"), bp::arg("model"),
             bp::arg("joint1_id"), bp::arg("joint1_placement")),
            "Contructor from given joint index and placement of the first joint "
            "implied in the constraint."))
          .def(bp::init<const Model &, JointIndex, const JointIndex &>(
            (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id"), bp::arg("joint2_id")),
            "Contructor from given joint index for the two joints "
            "implied in the constraint."))
          .def(bp::init<const Model &, JointIndex&>(
            (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id")),
            "Contructor from given joint index of the first joint "
            "implied in the constraint."))
          .PINOCCHIO_ADD_PROPERTY(T, joint1_id, "Index of the first joint in the model tree")
          .PINOCCHIO_ADD_PROPERTY(T, joint2_id, "Index of the second joint in the model tree")
          .PINOCCHIO_ADD_PROPERTY(T, joint1_placement, "Position of attached point with respect to the frame of joint1.")
          .PINOCCHIO_ADD_PROPERTY(T, joint2_placement, "Position of attached point with respect to the frame of joint2.")
          .PINOCCHIO_ADD_PROPERTY(T, desired_constraint_offset, "Desired constraint shift at position level.")
          .PINOCCHIO_ADD_PROPERTY(T, desired_constraint_velocity, "Desired constraint velocity at velocity level.")
          .PINOCCHIO_ADD_PROPERTY(T, desired_constraint_acceleration, "Desired constraint velocity at acceleration level.")
          .PINOCCHIO_ADD_PROPERTY(T, corrector_parameters, "Corrector parameters.")
          .PINOCCHIO_ADD_PROPERTY(T, colwise_joint1_sparsity, "Colwise sparsity pattern associated with joint 1.")
          .PINOCCHIO_ADD_PROPERTY(T, colwise_joint2_sparsity, "Colwise sparsity pattern associated with joint 2.")
          .PINOCCHIO_ADD_PROPERTY(T, joint1_span_indexes, " Jointwise span indexes associated with joint 1.")
          .PINOCCHIO_ADD_PROPERTY(T, joint2_span_indexes, "Jointwise span indexes associated with joint 2.")
          .PINOCCHIO_ADD_PROPERTY(T, loop_span_indexes, "Loop span indexes.")
          .PINOCCHIO_ADD_PROPERTY(T, colwise_sparsity, "Sparsity pattern associated to the constraint.")
          .PINOCCHIO_ADD_PROPERTY(T, colwise_span_indexes, "Indexes of the columns spanned by the constraints.")
          // .def("getRowSparsityPattern", &T::getRowSparsityPattern, args...)
          // .def("getRowActiveIndexes", &T::computeConstraintSpatialInertia , args...)
          // .def("getA1", &T::getA1, args...)
          // .def("getA2", &T::getA2, args...)
          // .def("computeConstraintSpatialInertia", &T::computeConstraintSpatialInertia, args...)
          // .def("appendConstraintDiagonalInertiaToJointInertias", &T::appendConstraintDiagonalInertiaToJointInertias, args...)
          // .def("mapConstraintForceToJointForces", &T::mapConstraintForceToJointForces, args...)
          ;
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_model_inheritance_hpp__
