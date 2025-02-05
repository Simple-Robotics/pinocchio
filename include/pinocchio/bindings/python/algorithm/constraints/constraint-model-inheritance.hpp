//
// Copyright (c) 2025 INRIA
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
      void visit(PyClass &) const
      {
      }
    };

    // Specialize
    template<class T>
    struct ConstraintModelInheritancePythonVisitor<T, FrameConstraintModelBase<T>>
    : public bp::def_visitor<
        ConstraintModelInheritancePythonVisitor<T, FrameConstraintModelBase<T>>>
    {
      typedef typename T::Scalar Scalar;
      typedef typename T::ConstraintSet ConstraintSet;
      typedef typename T::ConstraintData ConstraintData;
      typedef context::Model Model;
      typedef context::Data Data;
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def(bp::init<const Model &, JointIndex, const SE3 &, JointIndex, const SE3 &>(
            (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id"),
            bp::arg("joint1_placement"), bp::arg("joint2_id"), bp::arg("joint2_placement")),
            "Contructor from given joint index and placement for the two joints "
            "implied in the constraint."))
          .def(bp::init<const Model &, JointIndex, const SE3 &>(
            (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id"), bp::arg("joint1_placement")),
            "Contructor from given joint index and placement of the first joint "
            "implied in the constraint."))
          .def(bp::init<const Model &, JointIndex, JointIndex>(
            (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id"), bp::arg("joint2_id")),
            "Contructor from given joint index for the two joints "
            "implied in the constraint."))
          .def(bp::init<const Model &, JointIndex>(
            (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id")),
            "Contructor from given joint index of the first joint "
            "implied in the constraint."))
          .PINOCCHIO_ADD_PROPERTY(T, joint1_id, "Index of the first joint in the model tree.")
          .PINOCCHIO_ADD_PROPERTY(T, joint2_id, "Index of the second joint in the model tree.")
          .PINOCCHIO_ADD_PROPERTY(
            T, joint1_placement, "Position of attached point with respect to the frame of joint1.")
          .PINOCCHIO_ADD_PROPERTY(
            T, joint2_placement, "Position of attached point with respect to the frame of joint2.")
          .PINOCCHIO_ADD_PROPERTY(
            T, desired_constraint_offset, "Desired constraint shift at position level.")
          .PINOCCHIO_ADD_PROPERTY(
            T, desired_constraint_velocity, "Desired constraint velocity at velocity level.")
          .PINOCCHIO_ADD_PROPERTY(
            T, desired_constraint_acceleration,
            "Desired constraint velocity at acceleration level.")
          .PINOCCHIO_ADD_PROPERTY(T, corrector_parameters, "Corrector parameters.")
          .PINOCCHIO_ADD_PROPERTY(
            T, colwise_joint1_sparsity, "Colwise sparsity pattern associated with joint 1.")
          .PINOCCHIO_ADD_PROPERTY(
            T, colwise_joint2_sparsity, "Colwise sparsity pattern associated with joint 2.")
          .PINOCCHIO_ADD_PROPERTY(
            T, joint1_span_indexes, "Jointwise span indexes associated with joint 1.")
          .PINOCCHIO_ADD_PROPERTY(
            T, joint2_span_indexes, "Jointwise span indexes associated with joint 2.")
          .PINOCCHIO_ADD_PROPERTY(T, loop_span_indexes, "Loop span indexes.")
          .PINOCCHIO_ADD_PROPERTY(
            T, colwise_sparsity, "Sparsity pattern associated to the constraint.")
          .PINOCCHIO_ADD_PROPERTY(
            T, colwise_span_indexes, "Indexes of the columns spanned by the constraints.")
          .def("getA1", &getA1, bp::args("self", "constraint_data", "reference_frame"),
            "Returns the constraint projector associated with joint 1. "
            "This matrix transforms a spatial velocity expressed in a reference frame "
            "to the first component of the constraint associated with joint 1.")
          .def("getA2", &getA2, bp::args("self", "constraint_data", "reference_frame"),
            "Returns the constraint projector associated with joint 2. "
            "This matrix transforms a spatial velocity expressed in a reference frame "
            "to the first component of the constraint associated with joint 2.")
          ;
      }

      static context::Matrix6s getA1(const T & self, const ConstraintData & constraint_data, ReferenceFrame rf)
      {
        context::Matrix6s res;
        switch(rf) {
          case WORLD:
            res = self.getA1(constraint_data, WorldFrame());
          case LOCAL:
            res = self.getA1(constraint_data, LocalFrame());
          case LOCAL_WORLD_ALIGNED:
            res = self.getA1(constraint_data, LocalWorldAlignedFrame());
        }
        return res;
      }

      static context::Matrix6s getA2(const T & self, const ConstraintData & constraint_data, ReferenceFrame rf)
      {
        context::Matrix6s res;
        switch(rf) {
          case WORLD:
            res = self.getA2(constraint_data, WorldFrame());
          case LOCAL:
            res = self.getA2(constraint_data, LocalFrame());
          case LOCAL_WORLD_ALIGNED:
            res = self.getA2(constraint_data, LocalWorldAlignedFrame());
        }
        return res;
      }
    };


    template<class T>
    struct ConstraintModelInheritancePythonVisitor<T, PointConstraintModelBase<T>>
    : public bp::def_visitor<
        ConstraintModelInheritancePythonVisitor<T, PointConstraintModelBase<T>>>
    {
      typedef typename T::Scalar Scalar;
      typedef typename T::ConstraintSet ConstraintSet;
      typedef typename T::ConstraintData ConstraintData;
      typedef context::Model Model;
      typedef context::Data Data;
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def(bp::init<const Model &, JointIndex, const SE3 &, JointIndex, const SE3 &>(
                 (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id"),
                  bp::arg("joint1_placement"), bp::arg("joint2_id"), bp::arg("joint2_placement")),
                 "Contructor from given joint index and placement for the two joints "
                 "implied in the constraint."))
          .def(bp::init<const Model &, JointIndex, const SE3 &>(
            (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id"), bp::arg("joint1_placement")),
            "Contructor from the given joint index and the placement wrt the first joint "
            "implied in the constraint."))
          .def(bp::init<const Model &, JointIndex, JointIndex>(
            (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id"), bp::arg("joint2_id")),
            "Contructor from given joint indexes for the two joints "
            "implied in the constraint."))
          .def(bp::init<const Model &, JointIndex>(
            (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id")),
            "Contructor from given joint index of the first joint "
            "implied in the constraint."))
          .PINOCCHIO_ADD_PROPERTY(T, joint1_id, "Index of the first joint in the model tree.")
          .PINOCCHIO_ADD_PROPERTY(T, joint2_id, "Index of the second joint in the model tree.")
          .PINOCCHIO_ADD_PROPERTY(
            T, joint1_placement, "Position of attached point with respect to the frame of joint1.")
          .PINOCCHIO_ADD_PROPERTY(
            T, joint2_placement, "Position of attached point with respect to the frame of joint2.")
          .PINOCCHIO_ADD_PROPERTY(
            T, desired_constraint_offset, "Desired constraint shift at position level.")
          .PINOCCHIO_ADD_PROPERTY(
            T, desired_constraint_velocity, "Desired constraint velocity at velocity level.")
          .PINOCCHIO_ADD_PROPERTY(
            T, desired_constraint_acceleration,
            "Desired constraint velocity at acceleration level.")
          .PINOCCHIO_ADD_PROPERTY(T, corrector_parameters, "Corrector parameters.")
          .PINOCCHIO_ADD_PROPERTY(
            T, colwise_joint1_sparsity, "Colwise sparsity pattern associated with joint 1.")
          .PINOCCHIO_ADD_PROPERTY(
            T, colwise_joint2_sparsity, "Colwise sparsity pattern associated with joint 2.")
          .PINOCCHIO_ADD_PROPERTY(
            T, joint1_span_indexes, "Jointwise span indexes associated with joint 1.")
          .PINOCCHIO_ADD_PROPERTY(
            T, joint2_span_indexes, "Jointwise span indexes associated with joint 2.")
          .PINOCCHIO_ADD_PROPERTY(T, loop_span_indexes, "Loop span indexes.")
          .PINOCCHIO_ADD_PROPERTY(
            T, colwise_sparsity, "Sparsity pattern associated to the constraint.")
          .PINOCCHIO_ADD_PROPERTY(
            T, colwise_span_indexes, "Indexes of the columns spanned by the constraints.")
          .def("getA1", &getA1, bp::args("self", "constraint_data", "reference_frame"),
            "Returns the constraint projector associated with joint 1. "
            "This matrix transforms a spatial velocity expressed in a reference frame "
            "to the first component of the constraint associated with joint 1.")
          .def("getA2", &getA2, bp::args("self", "constraint_data", "reference_frame"),
            "Returns the constraint projector associated with joint 2. "
            "This matrix transforms a spatial velocity expressed in a reference frame "
            "to the first component of the constraint associated with joint 2.")
          .def("computeConstraintSpatialInertia", &computeConstraintSpatialInertia,
            bp::args("self", "placement", "diagonal_constraint_inertia"),
            "This function computes the spatial inertia associated with the constraint."
          )
          // The two following methods are not exposed as they rely on allocators.
          // .def("appendConstraintDiagonalInertiaToJointInertias",
          //   &appendConstraintDiagonalInertiaToJointInertias,
          //   bp::args("self", "model", "data", "constraint_data", "diagonal_constraint_inertia", "inertias"),
          //   "Append the constraint diagonal inertia to the joint inertias."
          // )
          // .def("mapConstraintForceToJointForces",
          //   &mapConstraintForceToJointForces,
          //   bp::args("self", "model", "data", "constraint_data", "constraint_forces", "joint_forces"),
          //   "Map the constraint forces (aka constraint Lagrange multipliers) to the forces supported by the joints."
          // )
          ;
      }

      static context::Matrix36s getA1(const T & self, const ConstraintData & constraint_data, ReferenceFrame rf)
      {
        context::Matrix36s res;
        switch(rf) {
          case WORLD:
            res = self.getA1(constraint_data, WorldFrame());
          case LOCAL:
            res = self.getA1(constraint_data, LocalFrame());
          case LOCAL_WORLD_ALIGNED:
            res = self.getA1(constraint_data, LocalWorldAlignedFrame());
        }
        return res;
      }

      static context::Matrix36s getA2(const T & self, const ConstraintData & constraint_data, ReferenceFrame rf)
      {
        context::Matrix36s res;
        switch(rf) {
          case WORLD:
            res = self.getA2(constraint_data, WorldFrame());
          case LOCAL:
            res = self.getA2(constraint_data, LocalFrame());
          case LOCAL_WORLD_ALIGNED:
            res = self.getA2(constraint_data, LocalWorldAlignedFrame());
        }
        return res;
      }

      static context::Matrix6s computeConstraintSpatialInertia(const T & self, const context::SE3 & placement, const context::Vector3s & diagonal_constraint_inertia)
      {
        return self.computeConstraintSpatialInertia(placement, diagonal_constraint_inertia);
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_model_inheritance_hpp__
