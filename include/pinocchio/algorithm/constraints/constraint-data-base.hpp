//
// Copyright (c) 2023-2024 INRIA
//

#ifndef __pinocchio_algorithm_constraint_data_base_hpp__
#define __pinocchio_algorithm_constraint_data_base_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/common/data-entity.hpp"

namespace pinocchio
{

  template<typename Derived>
  struct ConstraintDataBase
  : NumericalBase<Derived>
  , DataEntity<Derived>
  {
    typedef typename traits<Derived>::Scalar Scalar;
    typedef typename traits<Derived>::ConstraintModel ConstraintModel;
    typedef DataEntity<Derived> Base;

    Derived & derived()
    {
      return static_cast<Derived &>(*this);
    }
    const Derived & derived() const
    {
      return static_cast<const Derived &>(*this);
    }

    std::string shortname() const
    {
      return derived().shortname();
    }
    static std::string classname()
    {
      return Derived::classname();
    }

    void disp(std::ostream & os) const
    {
      using namespace std;
      os << shortname() << endl;
    }

    friend std::ostream & operator<<(std::ostream & os, const ConstraintDataBase<Derived> & constraint)
    {
      constraint.disp(os);
      return os;
    }

    template<typename OtherDerived>
    bool operator==(const ConstraintDataBase<OtherDerived> &) const
    {
      return true;
    }

    template<typename OtherDerived>
    bool operator!=(const ConstraintDataBase<OtherDerived> & other) const
    {
      return !(*this == other);
    }

  protected:
    /// \brief Default constructor
    ConstraintDataBase()
    {
    }
  };
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraint_data_base_hpp__
