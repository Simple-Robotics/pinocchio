//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_baumgarte_corrector_parameters_hpp__
#define __pinocchio_python_algorithm_constraints_baumgarte_corrector_parameters_hpp__

#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hpp"

#include "pinocchio/bindings/python/utils/cast.hpp"
#include "pinocchio/bindings/python/utils/copyable.hpp"
#include "pinocchio/bindings/python/utils/comparable.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<typename BaumgarteCorrectorParameters>
    struct BaumgarteCorrectorParametersPythonVisitor
    : public boost::python::def_visitor<
        BaumgarteCorrectorParametersPythonVisitor<BaumgarteCorrectorParameters>>
    {
      typedef typename BaumgarteCorrectorParameters::Scalar Scalar;
      typedef BaumgarteCorrectorParameters Self;

    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def(bp::init<Scalar, Scalar>(bp::args("self", "Kp", "Kd"), "Default constructor."))
          .def_readwrite("Kp", &Self::Kp, "Proportional corrector value.")
          .def_readwrite("Kd", &Self::Kd, "Damping corrector value.")
          .def(ComparableVisitor<Self, pinocchio::is_floating_point<Scalar>::value>());
      }

      static void expose()
      {
        if (eigenpy::check_registration<BaumgarteCorrectorParameters>())
          return;
        bp::class_<BaumgarteCorrectorParameters>(
          "BaumgarteCorrectorParameters", "Parameters of the Baumgarte Corrector.", bp::no_init)
          .def(BaumgarteCorrectorParametersPythonVisitor());
      }
    };

  } // namespace python
} // namespace pinocchio

#endif // ifndef
       // __pinocchio_python_algorithm_constraints_baumgarte_corrector_parameters_hpp__
