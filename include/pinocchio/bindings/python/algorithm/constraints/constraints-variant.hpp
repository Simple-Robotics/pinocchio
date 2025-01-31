//
// Copyright (c) 2015-2025 CNRS INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_variant_hpp__
#define __pinocchio_python_algorithm_constraints_variant_hpp__

#include <boost/algorithm/string/replace.hpp>

#include "pinocchio/bindings/python/algorithm/constraints/constraint-model.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/constraint-data.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/constraints-models.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/constraints-datas.hpp"
#include "pinocchio/bindings/python/utils/std-aligned-vector.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<typename T>
    std::string sanitizedClassname()
    {
      std::string className = boost::replace_all_copy(T::classname(), "<", "_");
      boost::replace_all(className, ">", "");
      return className;
    }

    template<typename VariantType>
    struct ConstraintVariantVisitor : boost::static_visitor<PyObject *>
    {
      static result_type convert(VariantType const & jv)
      {
        return apply_visitor(ConstraintVariantVisitor<VariantType>(), jv);
      }

      template<typename T>
      result_type operator()(T const & t) const
      {
        return bp::incref(bp::object(t).ptr());
      }
    };

    struct ConstraintDataExposer
    {
      template<class T>
      void operator()(T)
      {
        expose_constraint_data<T>(
          expose_constraint_data_inheritance<T>(
            bp::class_<T>(
              sanitizedClassname<T>().c_str(), sanitizedClassname<T>().c_str(), bp::init<>())
              .def(ConstraintDataBasePythonVisitor<T>())
              .def(PrintableVisitor<T>())));
        bp::implicitly_convertible<T, context::ConstraintData>();
      }

      void operator()(boost::blank)
      {}
    };

    struct ConstraintModelExposer
    {
      template<class T>
      void operator()(T)
      {
        expose_constraint_model<T>(
          expose_constraint_model_inheritance<T>(
            bp::class_<T>(
              sanitizedClassname<T>().c_str(), sanitizedClassname<T>().c_str(), bp::no_init)
              .def(ConstraintModelBasePythonVisitor<T>())
              .def(PrintableVisitor<T>())));
        bp::implicitly_convertible<T, context::ConstraintModel>();
      }

      void operator()(boost::blank)
      {}
    };

    struct ConstraintStdExposer
    {
      template<class T>
      void operator()(T)
      {
        StdAlignedVectorPythonVisitor<T>::expose(
          "StdVec_" + sanitizedClassname<T>().c_str()
        );
      }

      void operator()(boost::blank)
      {}
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_variant_hpp__
