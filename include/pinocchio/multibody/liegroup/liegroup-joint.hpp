//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_multibody_liegroup_liegroup_joint_hpp__
#define __pinocchio_multibody_liegroup_liegroup_joint_hpp__

#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/multibody/joint/joint-composite.hpp"
#include "pinocchio/multibody/joint/joint-mimic.hpp"
#include "pinocchio/multibody/liegroup/liegroup-algo.hpp"

namespace pinocchio
{
  // Overload the composite using the dispatch of the Lie group algo
  template<typename _Scalar, int _Options, template<typename S, int O> class JointCollectionTpl>
  template<typename LieGroupMap>
  typename LieGroupMap::template operation<
    JointModelCompositeTpl<_Scalar, _Options, JointCollectionTpl>>::type
  JointModelTpl<_Scalar, _Options, JointCollectionTpl>::lie_group_impl() const
  {
    typedef LieGroupInstanceStep<_Scalar, _Options, LieGroupMap> Algo;
    typedef typename Algo::LgType LgType;
    typedef typename Algo::Args Args;

    LgType res;
    ::pinocchio::details::Dispatch<Algo>::run(*this, ArgsType(res))

      return res;
  }

  // Write a Lie group for the generic overload
  template<typename _Scalar, int _Options, typename LieGroupMap>
  struct JointModelLieGroupVisitor
  : boost::static_visitor<typename LieGroupMap::template product_variant<_Scalar, _Options>>
  {
    typedef typename LieGroupMap::template product_variant<_Scalar, _Options> LgType;

    // Default behavior, instantiate the type in an agregation
    template<typename JointModelDerived>
    LgType operator()(const JointModelBase<JointModelDerived> & jmodel) const
    {
      return LgType(jmodel.template lie_group<LieGroupMap>());
    }

    // Composite and Mimic are already agregated lie group type
    template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
    LgType
    operator()(const JointModelCompositeTpl<Scalar, Options, JointCollectionTpl> & jmodel) const
    {
      return jmodel.template lie_group<LieGroupMap>();
    }

    template<typename JointModelRef>
    LgType operator()(const JointModelMimic<JointModelRef> & jmodel) const
    {
      return jmodel.template lie_group<LieGroupMap>();
    }

    template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
    static std::string run(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel)
    {
      return boost::apply_visitor(JointModelLieGroupVisitor(), jmodel);
    }
  };

  // Overload the generic with a simple visiting
  template<typename _Scalar, int _Options, template<typename S, int O> class JointCollectionTpl>
  template<typename LieGroupMap>
  typename LieGroupMap::template operation<
    JointModelTpl<_Scalar, _Options, JointCollectionTpl>>::type
  JointModelTpl<_Scalar, _Options, JointCollectionTpl>::lie_group_impl() const
  {
    typedef JointModelLieGroupVisitor<_Scalar, _Options, LieGroupMap> Algo;
    return Algo::run(*this);
  }
} // namespace pinocchio

#endif // #ifndef __pinocchio_multibody_liegroup_liegroup_joint_hpp__
