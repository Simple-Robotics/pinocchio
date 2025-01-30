//
// Copyright (c) 2015-2025 CNRS INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_datas_hpp__
#define __pinocchio_python_algorithm_constraints_datas_hpp__

#include "pinocchio/algorithm/constraints/constraint-data-generic.hpp"
#include "pinocchio/algorithm/constraints/constraint-collection-default.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/constraint-data-inheritance.hpp"
#include "pinocchio/bindings/python/utils/printable.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    // Add the inheritance
    template<class T>
    inline bp::class_<T> & expose_constraint_data_inheritance(bp::class_<T> & cl)
    {
      return ConstraintDataInheritanceHelper<T, T::Base>::expose_inheritance();
    }

    // generic expose_constraint_data : do nothing special
    template<class T>
    inline bp::class_<T> & expose_constraint_data(bp::class_<T> & cl)
    {
      return cl
      ;
    }

    // specialization for ConstraintDataRevoluteUnaligned
    // template<>
    // inline bp::class_<ConstraintDataRevoluteUnaligned> &
    // expose_constraint_data<ConstraintDataRevoluteUnaligned>(bp::class_<ConstraintDataRevoluteUnaligned> & cl)
    // {
    //   return cl
    //     .def(bp::init<const Eigen::Vector3d &>(
    //       bp::args("axis"), "Init ConstraintDataRevoluteUnaligned from an axis with x-y-z components"));
    // };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_datas_hpp__
