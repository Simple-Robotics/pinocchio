//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_model_hpp__
#define __pinocchio_python_algorithm_constraints_model_hpp__

#include "pinocchio/algorithm/constraints/constraint-model-generic.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/constraint-model-base.hpp"
#include "pinocchio/bindings/python/utils/printable.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<typename ConstraintModel>
    struct ExtractConstraintModelVariantTypeVisitor
    {
      typedef typename ConstraintModel::ConstraintCollection ConstraintCollection;
      typedef bp::object result_type;

      template<typename ConstraintModelDerived>
      result_type operator()(const ConstraintModelBase<ConstraintModelDerived> & cmodel) const
      {
        bp::object obj(boost::ref(cmodel.derived()));
        return obj;
      }

      result_type operator()(boost::blank) const
      {
        bp::object obj;
        return obj;
      }

      static result_type extract(const ConstraintModel & cmodel)
      {
        return boost::apply_visitor(ExtractConstraintModelVariantTypeVisitor(), cmodel);
      }
    };

    template<typename ConstraintModel>
    struct ConstraintModelPythonVisitor
    {
      static void expose()
      {
        bp::class_<ConstraintModel>("ConstraintModel", "Generic Constraint Model", bp::no_init)
          // .def(bp::init<>(bp::arg("self"), "Default constructor"))
          .def(bp::init<const ConstraintModel &>(bp::args("self", "other"), "Copy constructor"))
          .def(ConstraintModelBasePythonVisitor<ConstraintModel>())
          .def(PrintableVisitor<ConstraintModel>())
          .def(
            "extract", ExtractConstraintModelVariantTypeVisitor<ConstraintModel>::extract,
            bp::arg("self"),
            "Returns a reference of the internal joint managed by the ConstraintModel",
            bp::with_custodian_and_ward_postcall<0, 1>());
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_model_hpp__
