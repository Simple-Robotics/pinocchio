//
// Copyright (c) 2019-2025 INRIA
//

#ifndef __pinocchio_multibody_joint_joint_common_operations_hpp__
#define __pinocchio_multibody_joint_joint_common_operations_hpp__

#include "pinocchio/macros.hpp"
#include "pinocchio/math/matrix.hpp"
#include "pinocchio/math/matrix-inverse.hpp"

namespace pinocchio
{

  ///
  /// \brief Linear affine transformation of the configuration vector.
  ///        Valide for most common joints which are evolving on a vector space.
  ///
  struct LinearAffineTransform
  {
    template<typename ConfigVectorIn, typename Scalar, typename ConfigVectorOut>
    static void run(
      const Eigen::MatrixBase<ConfigVectorIn> & qIn,
      const Scalar & scaling,
      const Scalar & offset,
      const Eigen::MatrixBase<ConfigVectorOut> & qOut)
    {
      assert(qIn.size() == qOut.size());
      qOut.const_cast_derived().noalias() =
        scaling * qIn + ConfigVectorOut::Constant(qOut.size(), offset);
    }
  };

  struct UnboundedRevoluteAffineTransform
  {
    template<typename ConfigVectorIn, typename Scalar, typename ConfigVectorOut>
    static void run(
      const Eigen::MatrixBase<ConfigVectorIn> & qIn,
      const Scalar & scaling,
      const Scalar & offset,
      const Eigen::MatrixBase<ConfigVectorOut> & qOut)
    {
      assert(qIn.size() == 2);
      assert(qOut.size() == 2);

      const typename ConfigVectorIn::Scalar & ca = qIn(0);
      const typename ConfigVectorIn::Scalar & sa = qIn(1);

      const typename ConfigVectorIn::Scalar & theta = math::atan2(sa, ca);
      const typename ConfigVectorIn::Scalar & theta_transform = scaling * theta + offset;

      auto & dest_ = qOut.const_cast_derived();
      SINCOS(theta_transform, &dest_.coeffRef(1), &dest_.coeffRef(0));
    }
  };

  struct NoAffineTransform
  {
    template<typename ConfigVectorIn, typename Scalar, typename ConfigVectorOut>
    static void run(
      const Eigen::MatrixBase<ConfigVectorIn> &,
      const Scalar &,
      const Scalar &,
      const Eigen::MatrixBase<ConfigVectorOut> &)
    {
      assert(false && "Joint cannot be used with JointMimic.");
    }
  };

  ///
  /// \brief Assign the correct configuration vector space affine transformation according to the
  /// joint type. Must be specialized for every joint type.
  ///
  template<typename Joint>
  struct ConfigVectorAffineTransform
  {
    typedef NoAffineTransform Type;
  };

} // namespace pinocchio

#endif // ifndef __pinocchio_multibody_joint_joint_common_operations_hpp__
