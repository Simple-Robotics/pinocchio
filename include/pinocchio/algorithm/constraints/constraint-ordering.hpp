//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_constraint_ordering_hpp__
#define __pinocchio_algorithm_constraints_constraint_ordering_hpp__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"

#include "pinocchio/algorithm/check.hpp"
#include "pinocchio/algorithm/constraints/constraints.hpp"

namespace pinocchio
{

  ///
  /// \brief Init the data according to the contact information contained in constraint_models.
  ///
  /// \tparam JointCollection Collection of Joint types.
  /// \tparam Allocator Allocator class for the std::vector.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_models Vector of contact information related to the problem.
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator>
  inline void computeJointMinimalOrdering(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models);

} // namespace pinocchio

#include "pinocchio/algorithm/constraints/constraint-ordering.hxx"

#endif // ifndef __pinocchio_algorithm_constraints_constraint_ordering_hpp__
