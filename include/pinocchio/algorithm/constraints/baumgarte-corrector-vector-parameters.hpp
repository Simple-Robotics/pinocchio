//
// Copyright (c) 2020-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_baumgarte_corrector_vector_parameters_hpp__
#define __pinocchio_algorithm_constraints_baumgarte_corrector_vector_parameters_hpp__

#include "pinocchio/algorithm/fwd.hpp"

namespace pinocchio
{

  template<typename VectorType>
  struct BaumgarteCorrectorVectorParametersTpl;

  template<typename NewScalar, typename VectorType>
  struct CastType<NewScalar, BaumgarteCorrectorVectorParametersTpl<VectorType>>
  {
    enum
    {
      RowsAtCompileTime = VectorType::RowsAtCompileTime,
      ColsAtCompileTime = VectorType::ColsAtCompileTime,
      Options = VectorType::Options
    };

    typedef Eigen::Matrix<NewScalar, RowsAtCompileTime, ColsAtCompileTime, Options> NewVectorType;
    typedef BaumgarteCorrectorVectorParametersTpl<NewVectorType> type;
  };

  template<typename _VectorType>
  struct traits<BaumgarteCorrectorVectorParametersTpl<_VectorType>>
  {
    typedef _VectorType VectorType;
    typedef typename VectorType::Scalar Scalar;
  };

  template<typename _VectorType>
  struct BaumgarteCorrectorVectorParametersTpl
  : NumericalBase<BaumgarteCorrectorVectorParametersTpl<_VectorType>>
  {
    typedef _VectorType VectorType;
    typedef typename VectorType::Scalar Scalar;

    template<typename OtherVectorType>
    friend struct BaumgarteCorrectorVectorParametersTpl;

    /// \brief Default constructor with 0-size Kp and Kd.
    /// It is needed for constraints that don't have baumgarte correction.
    BaumgarteCorrectorVectorParametersTpl()
    {
    }

    explicit BaumgarteCorrectorVectorParametersTpl(int size)
    : Kp(size)
    , Kd(size)
    {
      Kp.setZero();
      Kd.setZero();
    }

    /// \brief Constructor from VectorType.
    /// It is needed for the generic constraint model.
    template<typename Vector1Like, typename Vector2Like>
    BaumgarteCorrectorVectorParametersTpl(
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
    BaumgarteCorrectorVectorParametersTpl<Eigen::Ref<OtherVectorType>> ref()
    {
      typedef BaumgarteCorrectorVectorParametersTpl<Eigen::Ref<OtherVectorType>> ReturnType;
      ReturnType res(::pinocchio::make_ref(Kp), ::pinocchio::make_ref(Kd));
      return res;
    }

    /// \brief Get const reference to baumgarte parameters.
    /// It is needed for the generic constraint model.
    template<typename OtherVectorType>
    BaumgarteCorrectorVectorParametersTpl<Eigen::Ref<const OtherVectorType>> ref() const
    {
      typedef BaumgarteCorrectorVectorParametersTpl<Eigen::Ref<const OtherVectorType>> ReturnType;
      ReturnType res(::pinocchio::make_const_ref(Kp), ::pinocchio::make_const_ref(Kd));
      return res;
    }

    bool operator==(const BaumgarteCorrectorVectorParametersTpl & other) const
    {
      if (this == &other)
        return true;
      return Kp == other.Kp && Kd == other.Kd;
    }

    bool operator!=(const BaumgarteCorrectorVectorParametersTpl & other) const
    {
      return !(*this == other);
    }

    template<typename OtherVectorType>
    BaumgarteCorrectorVectorParametersTpl &
    operator=(const BaumgarteCorrectorVectorParametersTpl<OtherVectorType> & other)
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
    typename CastType<NewScalar, BaumgarteCorrectorVectorParametersTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, BaumgarteCorrectorVectorParametersTpl>::type ReturnType;
      ReturnType res;
      res.Kp = Kp.template cast<NewScalar>();
      res.Kd = Kd.template cast<NewScalar>();
      return res;
    }

  }; // struct BaumgarteCorrectorVectorParametersTpl
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_baumgarte_corrector_vector_parameters_hpp__
