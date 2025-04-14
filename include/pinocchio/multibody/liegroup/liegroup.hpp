//
// Copyright (c) 2016-2019 CNRS INRIA
//

#ifndef __pinocchio_lie_group_hpp__
#define __pinocchio_lie_group_hpp__

#include "pinocchio/multibody/liegroup/vector-space.hpp"
#include "pinocchio/multibody/liegroup/cartesian-product.hpp"
#include "pinocchio/multibody/liegroup/cartesian-product-variant.hpp"
#include "pinocchio/multibody/liegroup/special-orthogonal.hpp"
#include "pinocchio/multibody/liegroup/special-euclidean.hpp"
#include "pinocchio/multibody/liegroup/liegroup-collection.hpp"

#include "pinocchio/multibody/joint/fwd.hpp"
// #include "pinocchio/multibody/joint/joint-mimic.hpp"

namespace pinocchio
{
  struct LieGroupMap
  {

    // Default LieGroup operation is euclidean
    template<typename JointModel>
    struct operation
    {
      typedef VectorSpaceOperationTpl<
        JointModel::NQ,
        typename JointModel::Scalar,
        JointModel::Options>
        type;
    };

    // LieGroupCollectionDefaultTpl is implicitely chosen by the LieGroupMap
    // inside the aggregation type
    template<typename Scalar, int Options>
    struct product_variant
    {
      typedef CartesianProductOperationVariantTpl<Scalar, Options, LieGroupCollectionDefaultTpl>
        type;
    }
  };

  // template<typename JointModel>
  // struct LieGroup
  // {
  //   typedef typename LieGroupMap::operation<JointModel>::type type;
  // };

  // The type for generic, composite and mimic is directly the variant type
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  struct LieGroupMap::operation<JointModelTpl<Scalar, Options, JointCollectionTpl>>
  {
    typedef typename LieGroupMap::product_variant<Scalar, Options>::type type;
  };

  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  struct LieGroupMap::operation<JointModelCompositeTpl<Scalar, Options, JointCollectionTpl>>
  {
    typedef typename LieGroupMap::product_variant<Scalar, Options>::type type;
  };

  // TODO: Fix after the rebase on Mimic
  // template<typename JointModelRef>
  // struct LieGroupMap::operation<JointModelMimic<JointModelRef>>
  // {
  //   typedef typename LieGroupMap::product_variant<typename JointModelRef::Scalar,
  //   JointModelRef::Options>::type type;
  // };

  // Atomic joint that are not euclidean
  template<typename Scalar, int Options>
  struct LieGroupMap::operation<JointModelSphericalTpl<Scalar, Options>>
  {
    typedef SpecialOrthogonalOperationTpl<3, Scalar, Options> type;
  };

  template<typename Scalar, int Options>
  struct LieGroupMap::operation<JointModelFreeFlyerTpl<Scalar, Options>>
  {
    typedef SpecialEuclideanOperationTpl<3, Scalar, Options> type;
  };

  template<typename Scalar, int Options>
  struct LieGroupMap::operation<JointModelPlanarTpl<Scalar, Options>>
  {
    typedef SpecialEuclideanOperationTpl<2, Scalar, Options> type;
  };

  template<typename Scalar, int Options, int axis>
  struct LieGroupMap::operation<JointModelRevoluteUnboundedTpl<Scalar, Options, axis>>
  {
    typedef SpecialOrthogonalOperationTpl<2, Scalar, Options> type;
  };

  template<typename Scalar, int Options>
  struct LieGroupMap::operation<JointModelRevoluteUnboundedUnalignedTpl<Scalar, Options>>
  {
    typedef SpecialOrthogonalOperationTpl<2, Scalar, Options> type;
  };

} // namespace pinocchio

#endif // ifndef __pinocchio_lie_group_hpp__
