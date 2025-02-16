//
// Copyright (c) 2020-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_baumgarte_corrector_parameters_hpp__
#define __pinocchio_algorithm_constraints_baumgarte_corrector_parameters_hpp__

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
      ColsAtCompileTime = BaumgarteVector::ColsAtCompileTime,
      Options = BaumgarteVector::Options
    };

    typedef Eigen::Matrix<NewScalar, RowsAtCompileTime, ColsAtCompileTime, Options>
      NewBaumgarteVector;
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

    template<typename OtherBaumgarteVector>
    friend struct BaumgarteCorrectorParametersTpl;

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

    /// \brief Get reference to baumgarte parameters.
    /// It is needed for the generic constraint model.
    template<typename OtherBaumgarteVector>
    BaumgarteCorrectorParametersTpl<Eigen::Ref<OtherBaumgarteVector>> get_ref()
    {
      typedef BaumgarteCorrectorParametersTpl<Eigen::Ref<OtherBaumgarteVector>> ReturnType;
      ReturnType res(::pinocchio::make_ref(Kp), ::pinocchio::make_ref(Kd));
      return res;
    }

    /// \brief Get const reference to baumgarte parameters.
    /// It is needed for the generic constraint model.
    template<typename OtherBaumgarteVector>
    BaumgarteCorrectorParametersTpl<Eigen::Ref<const OtherBaumgarteVector>> get_const_ref() const
    {
      typedef BaumgarteCorrectorParametersTpl<Eigen::Ref<const OtherBaumgarteVector>> ReturnType;
      ReturnType res(::pinocchio::make_const_ref(Kp), ::pinocchio::make_const_ref(Kd));
      return res;
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

  protected:
    /// \brief Constructor from BaumgarteVector.
    /// It is needed for the generic constraint model.
    BaumgarteCorrectorParametersTpl(const BaumgarteVector & Kp, const BaumgarteVector & Kd)
    : Kp(Kp)
    , Kd(Kd)
    {
    }

  }; // struct BaumgarteCorrectorParametersTpl
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_baumgarte_corrector_parameters_hpp__
