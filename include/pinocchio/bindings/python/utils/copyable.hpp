//
// Copyright (c) 2016-2023 CNRS INRIA
//

#ifndef __pinocchio_python_utils_copyable_hpp__
#define __pinocchio_python_utils_copyable_hpp__

#include <eigenpy/copyable.hpp>

#include "pinocchio/deprecated.hpp"

PINOCCHIO_DEPRECATED_HEADER("Directly include <eigenpy/copyable.hpp>")

namespace pinocchio
{
  namespace python
  {
    using eigenpy::CopyableVisitor;
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_utils_copyable_hpp__
