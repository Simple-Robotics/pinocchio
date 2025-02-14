//
// Copyright (c) 2020-2024 INRIA
//

#ifndef __pinocchio_algorithm_baumgarte_corrector_parameters_hpp__
#define __pinocchio_algorithm_baumgarte_corrector_parameters_hpp__

#include "pinocchio/algorithm/fwd.hpp"

namespace pinocchio
{

  template<typename _BaumgarteVector>
  struct BaumgarteCorrectorParametersTpl;

  template<typename NewScalar, typename BaumgarteVector>
  struct CastType<NewScalar, BaumgarteCorrectorParametersTpl<BaumgarteVector>>
  {
    enum
    {
      RowsAtCompileTime = BaumgarteVector::RowsAtCompileTime,
      ColsAtCompileTime = BaumgarteVector::ColsAtCompileTime
    };

    typedef Eigen::Matrix<NewScalar, RowsAtCompileTime, ColsAtCompileTime> NewBaumgarteVector;
    typedef BaumgarteCorrectorParametersTpl<NewBaumgarteVector> type;
  };

  template<typename _BaumgarteVector>
  struct traits<BaumgarteCorrectorParametersTpl<_BaumgarteVector>>
  {
    typedef _BaumgarteVector BaumgarteVector;
    typedef typename BaumgarteVector::Scalar Scalar;
  };

  template<typename _BaumgarteVector>
  struct BaumgarteCorrectorParametersTpl
  : NumericalBase<BaumgarteCorrectorParametersTpl<_BaumgarteVector>>
  {
    typedef _BaumgarteVector BaumgarteVector;
    typedef typename BaumgarteVector::Scalar Scalar;

    /// \brief Default constructor which does not set Kp/Kd.
    /// It is needed for constraints that don't have baumgarte correction.
    BaumgarteCorrectorParametersTpl()
    {
    }

    explicit BaumgarteCorrectorParametersTpl(int size)
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

    template<typename OtherBaumgarteVector>
    BaumgarteCorrectorParametersTpl &
    operator=(const BaumgarteCorrectorParametersTpl<OtherBaumgarteVector> & other)
    {
      Kp = other.Kp;
      Kd = other.Kd;
      return *this;
    }

    // parameters
    /// \brief Proportional corrector value.
    BaumgarteVector Kp;

    /// \brief Damping corrector value.
    BaumgarteVector Kd;

    template<typename NewScalar>
    typename CastType<NewScalar, BaumgarteCorrectorParametersTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, BaumgarteCorrectorParametersTpl>::type ReturnType;
      ReturnType res;
      res.Kp = Kp.template cast<NewScalar>();
      res.Kd = Kd.template cast<NewScalar>();
      return res;
    }

  }; // struct BaumgarteCorrectorParametersTpl
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_baumgarte_corrector_parameters_hpp__
