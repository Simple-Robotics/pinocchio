//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_baumgarte_corrector_vector_parameters_hpp__
#define __pinocchio_python_algorithm_constraints_baumgarte_corrector_vector_parameters_hpp__

#include "pinocchio/algorithm/constraints/baumgarte-corrector-vector-parameters.hpp"

#include "pinocchio/bindings/python/utils/cast.hpp"
#include "pinocchio/bindings/python/utils/copyable.hpp"
#include "pinocchio/bindings/python/utils/comparable.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<typename BaumgarteCorrectorVectorParameters>
    struct BaumgarteCorrectorVectorParametersPythonVisitor
    : public boost::python::def_visitor<
        BaumgarteCorrectorVectorParametersPythonVisitor<BaumgarteCorrectorVectorParameters>>
    {
      typedef typename BaumgarteCorrectorVectorParameters::Scalar Scalar;
      typedef typename BaumgarteCorrectorVectorParameters::VectorType VectorType;
      typedef BaumgarteCorrectorVectorParameters Self;

    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl
          // .def(bp::init<int>(bp::args("self", "size"), "Default constructor."))

          .def_readwrite("Kp", &Self::Kp, "Proportional corrector value.")
          .def_readwrite("Kd", &Self::Kd, "Damping corrector value.")

          // .def(CastVisitor<Self>())
          // .def(ExposeConstructorByCastVisitor<
          //      Self,
          //      ::pinocchio::context::RigidConstraintModel::BaumgarteCorrectorVectorParameters>())
          .def(ComparableVisitor<Self, pinocchio::is_floating_point<Scalar>::value>());
      }

      static void expose(const std::string & classname)
      {
        // eigenpy::enableEigenPySpecific<VectorType>();
        if (eigenpy::check_registration<BaumgarteCorrectorVectorParameters>())
          return;
        bp::class_<BaumgarteCorrectorVectorParameters>(
          classname.c_str(), "Parameters of the Baumgarte Corrector.", bp::no_init)
          .def(BaumgarteCorrectorVectorParametersPythonVisitor());
      }
    };

  } // namespace python
} // namespace pinocchio

#endif // ifndef
       // __pinocchio_python_algorithm_constraints_baumgarte_corrector_vector_parameters_hpp__
