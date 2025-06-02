//
// Copyright (c) 2020 INRIA
//

#ifndef __pinocchio_python_utils_deprecation_hpp__
#define __pinocchio_python_utils_deprecation_hpp__

#include <eigenpy/eigenpy.hpp>
#include <eigenpy/deprecation-policy.hpp>
#include "pinocchio/deprecated.hpp"

PINOCCHIO_DEPRECATED_HEADER("Include <eigenpy/deprecation-policy.hpp> instead.")

namespace pinocchio
{
  namespace python
  {
    template<class P = boost::python::default_call_policies>
    using deprecated_warning_policy =
      eigenpy::deprecation_warning_policy<eigenpy::DeprecationType::DEPRECATION, P>;
    using eigenpy::deprecated_function;
    using eigenpy::deprecated_member;
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_utils_deprecation_hpp__
