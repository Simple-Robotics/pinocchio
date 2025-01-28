//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_multibody_joint_generic_hxx__
#define __pinocchio_multibody_joint_generic_hxx__

namespace pinocchio
{

  namespace details
  {
    template<typename JointModelType, typename JointCollection>
    struct IsContainedInJointCollection
    {
      typedef typename JointCollection::JointModelVariant ModelVariant;
      static constexpr bool value =
        boost::mpl::contains<typename ModelVariant::types, JointModelType>();
    };

    template<bool CopyToJointCollection, typename JointModelDerived>
    struct CopyJointToJointCollection
    {
      template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
      void operator()(
        const JointModelBase<JointModelDerived> & joint_in,
        JointModelTpl<Scalar, Options, JointCollectionTpl> &)
      {
        std::stringstream ss;
        ss << joint_in.classname();
        ss << " not contained in new joint collection.\n";
        PINOCCHIO_THROW(std::invalid_argument, ss.str());
      }
    };

    template<typename JointModelDerived>
    struct CopyJointToJointCollection<true, JointModelDerived>
    {
      template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
      void operator()(
        const JointModelBase<JointModelDerived> & joint_in,
        JointModelTpl<Scalar, Options, JointCollectionTpl> & joint_out)
      {
        joint_out = joint_in;
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      template<typename, int> class OtherJointCollectionTpl>
    struct CopyJointFromOtherJointCollection
    : fusion::JointUnaryVisitorBase<CopyJointFromOtherJointCollection<
        Scalar,
        Options,
        JointCollectionTpl,
        OtherJointCollectionTpl>>
    {
      typedef JointCollectionTpl<Scalar, Options> JointCollection;
      typedef OtherJointCollectionTpl<Scalar, Options> OtherJointCollection;
      typedef boost::fusion::vector<JointModelTpl<Scalar, Options, JointCollectionTpl> &> ArgsType;

      template<typename JointModelDerived>
      static void algo(
        const JointModelDerived & joint_in,
        JointModelTpl<Scalar, Options, JointCollectionTpl> & joint_out)
      {
        CopyJointToJointCollection<
          IsContainedInJointCollection<JointModelDerived, JointCollection>::value,
          JointModelDerived>{}(joint_in, joint_out);
      }
    };

  } // namespace details

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  template<template<typename, int> class OtherJointCollectionTpl>
  JointModelTpl<Scalar, Options, JointCollectionTpl> &
  JointModelTpl<Scalar, Options, JointCollectionTpl>::operator=(
    const JointModelTpl<Scalar, Options, OtherJointCollectionTpl> & other)
  {
    typedef details::CopyJointFromOtherJointCollection<
      Scalar, Options, JointCollectionTpl, OtherJointCollectionTpl>
      Algo;
    typename Algo::ArgsType args(*this);
    Algo::run(other, args);
    return *this;
  }

} // namespace pinocchio

#endif // __pinocchio_multibody_joint_generic_hxx__
