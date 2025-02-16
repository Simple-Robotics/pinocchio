//
// Copyright (c) 2020-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_baumgarte_corrector_parameters_hpp__
#define __pinocchio_algorithm_constraints_baumgarte_corrector_parameters_hpp__

#include "pinocchio/algorithm/fwd.hpp"

namespace pinocchio
{

  template<typename VectorType>
  struct BaumgarteCorrectorParametersTpl;

  template<typename NewScalar, typename VectorType>
  struct CastType<NewScalar, BaumgarteCorrectorParametersTpl<VectorType>>
  {
    enum
    {
      RowsAtCompileTime = VectorType::RowsAtCompileTime,
      ColsAtCompileTime = VectorType::ColsAtCompileTime,
      Options = VectorType::Options
    };

    typedef Eigen::Matrix<NewScalar, RowsAtCompileTime, ColsAtCompileTime, Options> NewVectorType;
    typedef BaumgarteCorrectorParametersTpl<NewVectorType> type;
  };

  template<typename _VectorType>
  struct traits<BaumgarteCorrectorParametersTpl<_VectorType>>
  {
    typedef _VectorType VectorType;
    typedef typename VectorType::Scalar Scalar;
  };

  template<typename _VectorType>
  struct BaumgarteCorrectorParametersTpl
  : NumericalBase<BaumgarteCorrectorParametersTpl<_VectorType>>
  {
    typedef _VectorType VectorType;
    typedef typename VectorType::Scalar Scalar;

    template<typename OtherVectorType>
    friend struct BaumgarteCorrectorParametersTpl;

    /// \brief Default constructor with 0-size Kp and Kd.
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

    /// \brief Constructor from VectorType.
    /// It is needed for the generic constraint model.
    template<typename Vector1Like, typename Vector2Like>
    BaumgarteCorrectorParametersTpl(
      const Eigen::MatrixBase<Vector1Like> & Kp, const Eigen::MatrixBase<Vector2Like> & Kd)
    : Kp(Kp)
    , Kd(Kd)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        (this->Kp.array() >= Scalar(0)).all(), "Kp should only contain non negative quantities.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        (this->Kd.array() >= Scalar(0)).all(), "Kp should only contain non negative quantities.");
    }

    /// \brief Get reference to baumgarte parameters.
    /// It is needed for the generic constraint model.
    template<typename OtherVectorType>
    BaumgarteCorrectorParametersTpl<Eigen::Ref<OtherVectorType>> get_ref()
    {
      typedef BaumgarteCorrectorParametersTpl<Eigen::Ref<OtherVectorType>> ReturnType;
      ReturnType res(::pinocchio::make_ref(Kp), ::pinocchio::make_ref(Kd));
      return res;
    }

    /// \brief Get const reference to baumgarte parameters.
    /// It is needed for the generic constraint model.
    template<typename OtherVectorType>
    BaumgarteCorrectorParametersTpl<Eigen::Ref<const OtherVectorType>> get_const_ref() const
    {
      typedef BaumgarteCorrectorParametersTpl<Eigen::Ref<const OtherVectorType>> ReturnType;
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

    template<typename OtherVectorType>
    BaumgarteCorrectorParametersTpl &
    operator=(const BaumgarteCorrectorParametersTpl<OtherVectorType> & other)
    {
      Kp = other.Kp;
      Kd = other.Kd;
      return *this;
    }

    // parameters
    /// \brief Proportional corrector values.
    VectorType Kp;

    /// \brief Damping corrector values.
    VectorType Kd;

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

#endif // ifndef __pinocchio_algorithm_constraints_baumgarte_corrector_parameters_hpp__
