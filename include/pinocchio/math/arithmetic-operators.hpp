//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_math_arithmetic_operators_hpp__
#define __pinocchio_math_arithmetic_operators_hpp__

#include "pinocchio/math/matrix.hpp"

namespace pinocchio
{
  template<typename LhsType, typename RhsType>
  struct MultiplicationOperatorReturnType;
  
  template<typename LhsMatrixDerived, typename RhsMatrixDerived>
  struct MultiplicationOperatorReturnType<Eigen::MatrixBase<LhsMatrixDerived>,Eigen::MatrixBase<RhsMatrixDerived>>
  : MatrixMatrixProduct<LhsMatrixDerived,RhsMatrixDerived>
  {
    typedef MatrixMatrixProduct<LhsMatrixDerived,RhsMatrixDerived> Base;
    typedef typename Base::type type;
  };
}

#endif //#ifndef __pinocchio_math_arithmetic_operators_hpp__
