//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_data_hpp__
#define __pinocchio_python_algorithm_constraints_data_hpp__

#include "pinocchio/algorithm/constraints/constraint-data-generic.hpp"
#include "pinocchio/serialization/constraints-data.hpp"

#include "pinocchio/bindings/python/algorithm/constraints/constraint-data-base.hpp"
#include "pinocchio/bindings/python/utils/printable.hpp"
#include "pinocchio/bindings/python/serialization/serializable.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<typename ConstraintData>
    struct ExtractConstraintDataVariantTypeVisitor
    {
      typedef typename ConstraintData::ConstraintCollection ConstraintCollection;
      typedef bp::object result_type;

      template<typename ConstraintDataDerived>
      result_type operator()(const ConstraintDataBase<ConstraintDataDerived> & cdata) const
      {
        bp::object obj(boost::ref(cdata.derived()));
        return obj;
      }

      result_type operator()(boost::blank) const
      {
        bp::object obj;
        return obj;
      }

      static result_type extract(const ConstraintData & cdata)
      {
        return boost::apply_visitor(ExtractConstraintDataVariantTypeVisitor(), cdata);
      }
    };

    template<typename ConstraintData>
    struct ConstraintDataPythonVisitor
    {
      static void expose()
      {
        bp::class_<ConstraintData>("ConstraintData", "Generic Constraint Data.", bp::no_init)
          .def(bp::init<>(bp::arg("self"), "Default constructor."))
          .def(bp::init<const ConstraintData &>(
            bp::args("self", "constraint_data"), "Copy constructor."))
          .def(ConstraintDataBasePythonVisitor<ConstraintData>())
          .def(PrintableVisitor<ConstraintData>())
          .def(SerializableVisitor<ConstraintData>())
          .def(
            "extract", ExtractConstraintDataVariantTypeVisitor<ConstraintData>::extract,
            bp::arg("self"),
            "Returns a reference of the internal constraint managed by the ConstraintData.",
            bp::with_custodian_and_ward_postcall<0, 1>());
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_data_hpp__
