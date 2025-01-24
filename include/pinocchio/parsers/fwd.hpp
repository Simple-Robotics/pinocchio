//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_parsers_fwd_hpp__
#define __pinocchio_parsers_fwd_hpp__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/geometry.hpp"
#include "pinocchio/algorithm/constraints/point-bilateral-constraint.hpp"

namespace pinocchio
{
  namespace parsers
  {

    // Parsers will work on an object of type ::pinocchio::parsers::Model, which is a
    // specialization of the ModelTpl template of Pinocchio.
    // Once constructed, this model can then be cast to another scalar or to another joint
    // collection.
    typedef ::pinocchio::ModelTpl<double, 0, JointCollectionDefaultTpl> Model;
    typedef Model::JointModel JointModel;

  } // namespace parsers
} // namespace pinocchio

#endif // __pinocchio_parsers_fwd_hpp__
