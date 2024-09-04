//
// Copyright (c) 2020-2024 INRIA
//

#ifndef __pinocchio_algorithm_baumgarte_corrector_parameters_hpp__
#define __pinocchio_algorithm_baumgarte_corrector_parameters_hpp__

#include "pinocchio/algorithm/fwd.hpp"

namespace pinocchio
{

  template<typename _Scalar>
  struct BaumgarteCorrectorParametersTpl;

  template<typename _Scalar>
  struct traits<BaumgarteCorrectorParametersTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
  };

  template<typename _Scalar>
  struct BaumgarteCorrectorParametersTpl : NumericalBase<BaumgarteCorrectorParametersTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, -1, 1, Eigen::ColMajor, 6> Vector6Max;

    explicit BaumgarteCorrectorParametersTpl(int size = 6)
    : Kp(size)
    , Kd(size)
    {
      Kp.setZero();
      Kd.setZero();
    }

    bool operator==(const BaumgarteCorrectorParametersTpl & other) const
    {
      return Kp == other.Kp && Kd == other.Kd;
    }

    bool operator!=(const BaumgarteCorrectorParametersTpl & other) const
    {
      return !(*this == other);
    }

    // parameters
    /// \brief Proportional corrector value.
    Vector6Max Kp;

    /// \brief Damping corrector value.
    Vector6Max Kd;

    template<typename NewScalar>
    BaumgarteCorrectorParametersTpl<NewScalar> cast() const
    {
      typedef BaumgarteCorrectorParametersTpl<NewScalar> ReturnType;
      ReturnType res;
      res.Kp = Kp.template cast<NewScalar>();
      res.Kd = Kd.template cast<NewScalar>();
      return res;
    }
  }; // struct BaumgarteCorrectorParametersTpl
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_baumgarte_corrector_parameters_hpp__
