//
// Copyright (c) 2026 INRIA
//

#pragma once
// IWYU pragma: begin_keep
#include "pinocchio/macros.hpp"
#include "pinocchio/deprecated.hpp"

#include "pinocchio/spatial.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/collision/coal-pinocchio-conversions.hpp"
// IWYU pragma: end_exports

// clang-format off
PINOCCHIO_MOVED_HEADER_PINOCCHIO4(pinocchio/collision/fcl-pinocchio-conversions.hpp, pinocchio/collision/coal-pinocchio-conversions.hpp)
// clang-format on

namespace pinocchio
{

  template<typename Scalar>
  PINOCCHIO_DEPRECATED inline coal::Transform3s toFclTransform3f(const SE3Tpl<Scalar> & m)
  {
    return toCoalTransform3s(m);
  }

} // namespace pinocchio
