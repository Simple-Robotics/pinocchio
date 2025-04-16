//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_multibody_liegroup_liegroup_joint_hpp__
#define __pinocchio_multibody_liegroup_liegroup_joint_hpp__

#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/multibody/joint/joint-generic.hpp"
#include "pinocchio/multibody/liegroup/liegroup-algo.hpp"
#include "pinocchio/multibody/liegroup/liegroup-generic.hpp"

namespace pinocchio
{
  // Overload the composite using the dispatch of the Lie group algo
  template<typename _Scalar, int _Options, template<typename S, int O> class JointCollectionTpl>
  template<typename LieGroup_t>
  typename LieGroup_t::template operation<
    JointModelCompositeTpl<_Scalar, _Options, JointCollectionTpl>>::type
  JointModelCompositeTpl<_Scalar, _Options, JointCollectionTpl>::lie_group_impl() const
  {
    typedef LieGroupInstanceStep<LieGroup_t, _Scalar, _Options> Algo;
    typedef typename Algo::LgType LgType;
    typedef typename Algo::ArgsType ArgsType;

    LgType res;
    Algo::run(*this, ArgsType(res));
    return res;
  }

  // Write a Lie group for the generic overload
  template<typename LieGroup_t, typename _Scalar, int _Options>
  struct JointModelLieGroupVisitor
  : boost::static_visitor<typename LieGroup_t::template product_variant<_Scalar, _Options>::type>
  {
    typedef typename LieGroup_t::template product_variant<_Scalar, _Options>::type LgType;

    // Default behavior, instantiate the atomic type in the agregation type
    template<typename JointModelDerived>
    LgType operator()(const JointModelBase<JointModelDerived> & jmodel) const
    {
      return LgType(jmodel.template lie_group<LieGroup_t>());
    }

    // Composite and Mimic are already agregated lie group type
    template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
    LgType
    operator()(const JointModelCompositeTpl<Scalar, Options, JointCollectionTpl> & jmodel) const
    {
      return jmodel.template lie_group<LieGroup_t>();
    }

    template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
    LgType operator()(const JointModelMimicTpl<Scalar, Options, JointCollectionTpl> & jmodel) const
    {
      return jmodel.template lie_group<LieGroup_t>();
    }

    template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
    static LgType run(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel)
    {
      return boost::apply_visitor(JointModelLieGroupVisitor(), jmodel);
    }
  };

  // Overload the generic with a simple visiting
  template<typename _Scalar, int _Options, template<typename S, int O> class JointCollectionTpl>
  template<typename LieGroup_t>
  typename LieGroup_t::template operation<
    JointModelTpl<_Scalar, _Options, JointCollectionTpl>>::type
  JointModelTpl<_Scalar, _Options, JointCollectionTpl>::lie_group_impl() const
  {
    typedef JointModelLieGroupVisitor<LieGroup_t, _Scalar, _Options> Algo;
    return Algo::run(*this);
  }
} // namespace pinocchio

#endif // #ifndef __pinocchio_multibody_liegroup_liegroup_joint_hpp__
