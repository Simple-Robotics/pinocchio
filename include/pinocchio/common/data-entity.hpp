//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_common_data_entity_hpp__
#define __pinocchio_common_data_entity_hpp__

#include "pinocchio/common/fwd.hpp"

namespace pinocchio
{

  template<typename Derived>
  struct DataEntity
  {
    typedef typename traits<Derived>::Scalar Scalar;
    typedef typename traits<Derived>::Model Model;

    Derived & derived()
    {
      return static_cast<Derived &>(*this);
    }

    const Derived & derived() const
    {
      return static_cast<const Derived &>(*this);
    }
  };

} // namespace pinocchio

#endif // ifndef __pinocchio_common_data_entity_hpp__
