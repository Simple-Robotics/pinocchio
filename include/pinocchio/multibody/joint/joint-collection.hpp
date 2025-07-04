//
// Copyright (c) 2018 CNRS
// Copyright (c) 2018-2025 INRIA
//

#ifndef __pinocchio_multibody_joint_collection_hpp__
#define __pinocchio_multibody_joint_collection_hpp__

#include "pinocchio/multibody/joint/fwd.hpp"
#include "pinocchio/multibody/joint/joints.hpp"

#include <boost/variant.hpp>
#include <boost/variant/recursive_wrapper.hpp>

namespace pinocchio
{

  template<typename _Scalar, int _Options>
  struct JointCollectionDefaultTpl
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    // Joint Revolute
    typedef JointModelRevoluteTpl<Scalar, Options, 0> JointModelRX;
    typedef JointModelRevoluteTpl<Scalar, Options, 1> JointModelRY;
    typedef JointModelRevoluteTpl<Scalar, Options, 2> JointModelRZ;

    // Joint Revolute Unaligned
    typedef JointModelRevoluteUnalignedTpl<Scalar, Options> JointModelRevoluteUnaligned;

    // Joint Revolute UBounded
    typedef JointModelRevoluteUnboundedTpl<Scalar, Options, 0> JointModelRUBX;
    typedef JointModelRevoluteUnboundedTpl<Scalar, Options, 1> JointModelRUBY;
    typedef JointModelRevoluteUnboundedTpl<Scalar, Options, 2> JointModelRUBZ;

    // Joint Revolute Unbounded Unaligned
    typedef JointModelRevoluteUnboundedUnalignedTpl<Scalar, Options>
      JointModelRevoluteUnboundedUnaligned;

    // Joint Prismatic
    typedef JointModelPrismaticTpl<Scalar, Options, 0> JointModelPX;
    typedef JointModelPrismaticTpl<Scalar, Options, 1> JointModelPY;
    typedef JointModelPrismaticTpl<Scalar, Options, 2> JointModelPZ;

    // Joint Prismatic Unaligned
    typedef JointModelPrismaticUnalignedTpl<Scalar, Options> JointModelPrismaticUnaligned;

    // Joint Spherical
    typedef JointModelSphericalTpl<Scalar, Options> JointModelSpherical;

    // Joint Spherical ZYX
    typedef JointModelSphericalZYXTpl<Scalar, Options> JointModelSphericalZYX;

    // Joint Translation
    typedef JointModelTranslationTpl<Scalar, Options> JointModelTranslation;

    // Joint FreeFlyer
    typedef JointModelFreeFlyerTpl<Scalar, Options> JointModelFreeFlyer;

    // Joint Planar
    typedef JointModelPlanarTpl<Scalar, Options> JointModelPlanar;

    // Joint Composite
    typedef JointModelCompositeTpl<Scalar, Options, ::pinocchio::JointCollectionDefaultTpl>
      JointModelComposite;

    // Joint Mimic
    typedef JointModelMimicTpl<Scalar, Options, ::pinocchio::JointCollectionDefaultTpl>
      JointModelMimic;

    // Joint Helical
    typedef JointModelHelicalTpl<Scalar, Options, 0> JointModelHx;
    typedef JointModelHelicalTpl<Scalar, Options, 1> JointModelHy;
    typedef JointModelHelicalTpl<Scalar, Options, 2> JointModelHz;

    // Joint Helical Unaligned
    typedef JointModelHelicalUnalignedTpl<Scalar, Options> JointModelHelicalUnaligned;

    // Joint Universal
    typedef JointModelUniversalTpl<Scalar, Options> JointModelUniversal;

    typedef boost::variant<
      //    JointModelVoid,
      JointModelRX,
      JointModelRY,
      JointModelRZ,
      JointModelFreeFlyer,
      JointModelPlanar,
      JointModelRevoluteUnaligned,
      JointModelSpherical,
      JointModelSphericalZYX,
      JointModelPX,
      JointModelPY,
      JointModelPZ,
      JointModelPrismaticUnaligned,
      JointModelTranslation,
      JointModelRUBX,
      JointModelRUBY,
      JointModelRUBZ,
      JointModelRevoluteUnboundedUnaligned,
      JointModelHx,
      JointModelHy,
      JointModelHz,
      JointModelHelicalUnaligned,
      JointModelUniversal,
      boost::recursive_wrapper<JointModelComposite>,
      boost::recursive_wrapper<JointModelMimic>>
      JointModelVariant;

    // Joint Revolute
    typedef JointDataRevoluteTpl<Scalar, Options, 0> JointDataRX;
    typedef JointDataRevoluteTpl<Scalar, Options, 1> JointDataRY;
    typedef JointDataRevoluteTpl<Scalar, Options, 2> JointDataRZ;

    // Joint Revolute Unaligned
    typedef JointDataRevoluteUnalignedTpl<Scalar, Options> JointDataRevoluteUnaligned;

    // Joint Revolute Unaligned
    typedef JointDataRevoluteUnboundedUnalignedTpl<Scalar, Options>
      JointDataRevoluteUnboundedUnaligned;

    // Joint Revolute UBounded
    typedef JointDataRevoluteUnboundedTpl<Scalar, Options, 0> JointDataRUBX;
    typedef JointDataRevoluteUnboundedTpl<Scalar, Options, 1> JointDataRUBY;
    typedef JointDataRevoluteUnboundedTpl<Scalar, Options, 2> JointDataRUBZ;

    // Joint Prismatic
    typedef JointDataPrismaticTpl<Scalar, Options, 0> JointDataPX;
    typedef JointDataPrismaticTpl<Scalar, Options, 1> JointDataPY;
    typedef JointDataPrismaticTpl<Scalar, Options, 2> JointDataPZ;

    // Joint Prismatic Unaligned
    typedef JointDataPrismaticUnalignedTpl<Scalar, Options> JointDataPrismaticUnaligned;

    // Joint Spherical
    typedef JointDataSphericalTpl<Scalar, Options> JointDataSpherical;

    // Joint Spherical ZYX
    typedef JointDataSphericalZYXTpl<Scalar, Options> JointDataSphericalZYX;

    // Joint Translation
    typedef JointDataTranslationTpl<Scalar, Options> JointDataTranslation;

    // Joint FreeFlyer
    typedef JointDataFreeFlyerTpl<Scalar, Options> JointDataFreeFlyer;

    // Joint Planar
    typedef JointDataPlanarTpl<Scalar, Options> JointDataPlanar;

    // Joint Composite
    typedef JointDataCompositeTpl<Scalar, Options, ::pinocchio::JointCollectionDefaultTpl>
      JointDataComposite;

    // Joint Mimic
    typedef JointDataMimicTpl<Scalar, Options, ::pinocchio::JointCollectionDefaultTpl>
      JointDataMimic;

    // Joint Helical
    typedef JointDataHelicalTpl<Scalar, Options, 0> JointDataHx;
    typedef JointDataHelicalTpl<Scalar, Options, 1> JointDataHy;
    typedef JointDataHelicalTpl<Scalar, Options, 2> JointDataHz;

    // Joint Helical Unaligned
    typedef JointDataHelicalUnalignedTpl<Scalar, Options> JointDataHelicalUnaligned;

    // Joint Universal
    typedef JointDataUniversalTpl<Scalar, Options> JointDataUniversal;

    typedef boost::variant<
      //    JointDataVoid
      JointDataRX,
      JointDataRY,
      JointDataRZ,
      JointDataFreeFlyer,
      JointDataPlanar,
      JointDataRevoluteUnaligned,
      JointDataSpherical,
      JointDataSphericalZYX,
      JointDataPX,
      JointDataPY,
      JointDataPZ,
      JointDataPrismaticUnaligned,
      JointDataTranslation,
      JointDataRUBX,
      JointDataRUBY,
      JointDataRUBZ,
      JointDataRevoluteUnboundedUnaligned,
      JointDataHx,
      JointDataHy,
      JointDataHz,
      JointDataHelicalUnaligned,
      JointDataUniversal,
      boost::recursive_wrapper<JointDataComposite>,
      boost::recursive_wrapper<JointDataMimic>>
      JointDataVariant;
  };

  typedef JointCollectionDefault::JointModelVariant JointModelVariant;
  typedef JointCollectionDefault::JointDataVariant JointDataVariant;

} // namespace pinocchio

#endif // ifndef __pinocchio_multibody_joint_collection_hpp__
