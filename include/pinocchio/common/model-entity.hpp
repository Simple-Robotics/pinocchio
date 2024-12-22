//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_common_model_entity_hpp__
#define __pinocchio_common_model_entity_hpp__

#include "pinocchio/common/fwd.hpp"

namespace pinocchio
{

  template<typename Derived>
  struct ModelEntity
  {
    typedef typename traits<Derived>::Scalar Scalar;
    typedef typename traits<Derived>::Data Data;

    Derived & derived()
    {
      return static_cast<Derived &>(*this);
    }

    const Derived & derived() const
    {
      return static_cast<const Derived &>(*this);
    }

    Data createData() const
    {
      derived().createData();
    }
  };

  template<typename Derived>
  typename traits<Derived>::Data createData(const ModelEntity<Derived> & entity)
  {
    return entity.createData();
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_common_model_entity_hpp__
