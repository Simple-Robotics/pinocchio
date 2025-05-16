//
// Copyright (c) 2020-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_baumgarte_corrector_parameters_hpp__
#define __pinocchio_algorithm_constraints_baumgarte_corrector_parameters_hpp__

#include "pinocchio/fwd.hpp"

namespace pinocchio
{

  template<typename Scalar>
  struct BaumgarteCorrectorParametersTpl;

  template<typename NewScalar, typename Scalar>
  struct CastType<NewScalar, BaumgarteCorrectorParametersTpl<Scalar>>
  {
    typedef BaumgarteCorrectorParametersTpl<NewScalar> type;
  };

  template<typename _Scalar>
  struct traits<BaumgarteCorrectorParametersTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
  };

  template<typename _Scalar>
  struct BaumgarteCorrectorParametersTpl : NumericalBase<BaumgarteCorrectorParametersTpl<_Scalar>>
  {
    typedef _Scalar Scalar;

    /// \brief Default constructor initializes Kp and Kd to 0 (no correction).
    /// It is needed for constraints that don't have baumgarte correction.
    BaumgarteCorrectorParametersTpl()
    : Kp(Scalar(0))
    , Kd(Scalar(0))
    {
    }

    /// \brief Constructor from Kp and Kd.
    BaumgarteCorrectorParametersTpl(const Scalar Kp, const Scalar Kd)
    : Kp(Kp)
    , Kd(Kd)
    {
    }

    bool operator==(const BaumgarteCorrectorParametersTpl & other) const
    {
      if (this == &other)
        return true;
      return Kp == other.Kp && Kd == other.Kd;
    }

    bool operator!=(const BaumgarteCorrectorParametersTpl & other) const
    {
      return !(*this == other);
    }

    // parameters
    /// \brief Proportional corrector values.
    Scalar Kp;

    /// \brief Damping corrector values.
    Scalar Kd;

    template<typename NewScalar>
    typename CastType<NewScalar, BaumgarteCorrectorParametersTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, BaumgarteCorrectorParametersTpl>::type ReturnType;
      ReturnType res;
      res.Kp = NewScalar(Kp);
      res.Kd = NewScalar(Kd);
      return res;
    }

  }; // struct BaumgarteCorrectorParametersTpl
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_baumgarte_corrector_parameters_hpp__
