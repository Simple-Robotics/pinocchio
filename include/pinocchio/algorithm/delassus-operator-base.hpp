//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_base_hpp__
#define __pinocchio_algorithm_delassus_operator_base_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/math/eigenvalues.hpp"
#include "pinocchio/math/arithmetic-operators.hpp"

namespace pinocchio {
  
  template<typename DelassusOperatorDerived> struct DelassusOperatorBase;
  
  template<typename DelassusOperatorDerived, typename MatrixDerived>
  struct DelassusOperatorApplyOnTheRightReturnType;
  
  template<typename DelassusOperatorDerived, typename MatrixDerived>
  struct traits<DelassusOperatorApplyOnTheRightReturnType<DelassusOperatorDerived,MatrixDerived> >
  {
    typedef typename traits<DelassusOperatorDerived>::Scalar Scalar;
    typedef typename traits<DelassusOperatorDerived>::Matrix Matrix;
  };
  
  template<typename DelassusOperatorDerived, typename MatrixDerived>
  struct MultiplicationOperatorReturnType<DelassusOperatorBase<DelassusOperatorDerived>,Eigen::MatrixBase<MatrixDerived>>
  {
    typedef DelassusOperatorApplyOnTheRightReturnType<DelassusOperatorDerived,MatrixDerived> type;
  };
  
}

namespace Eigen {
  namespace internal {
    
    template<typename DelassusOperatorDerived, typename MatrixDerived>
    struct traits<pinocchio::DelassusOperatorApplyOnTheRightReturnType<DelassusOperatorDerived,MatrixDerived> >
    {
      typedef typename ::pinocchio::traits<DelassusOperatorDerived>::Matrix ReturnType;
      enum { Flags = 0 };
    };
  
    template< typename DstXprType, typename DelassusOperatorDerived, typename MatrixDerived, typename Functor>
    struct Assignment<DstXprType, Eigen::ReturnByValue<pinocchio::DelassusOperatorApplyOnTheRightReturnType<DelassusOperatorDerived,MatrixDerived>>, Functor, Dense2Dense, void>
    {
      typedef Eigen::ReturnByValue<pinocchio::DelassusOperatorApplyOnTheRightReturnType<DelassusOperatorDerived,MatrixDerived>> SrcXprType;
      
      EIGEN_DEVICE_FUNC
      static EIGEN_STRONG_INLINE void run(DstXprType &dst, const SrcXprType &src, const Functor &func)
      {
        Index dstRows = src.rows();
        Index dstCols = src.cols();
        if((dst.rows()!=dstRows) || (dst.cols()!=dstCols))
          dst.resize(dstRows, dstCols);
        
        eigen_assert(dst.rows() == src.rows() && dst.cols() == src.cols());
        src.evalTo(dst);
      }
    };
   
  } // namespace internal
} // namespace Eigen


namespace pinocchio {
  
  template<typename DelassusOperatorDerived, typename MatrixDerived>
  struct DelassusOperatorApplyOnTheRightReturnType
  : public Eigen::ReturnByValue<DelassusOperatorApplyOnTheRightReturnType<DelassusOperatorDerived,MatrixDerived> >
  {
    typedef DelassusOperatorApplyOnTheRightReturnType Self;
    
    DelassusOperatorApplyOnTheRightReturnType(const DelassusOperatorDerived & lhs,
                                              const MatrixDerived & rhs)
    : m_lhs(lhs)
    , m_rhs(rhs)
    {}
    
    template <typename ResultType>
    inline void evalTo(ResultType& result) const
    {
      m_lhs.applyOnTheRight(m_rhs.derived(),result);
    }
    
    EIGEN_CONSTEXPR Eigen::Index rows() const EIGEN_NOEXCEPT { return m_lhs.rows(); }
    EIGEN_CONSTEXPR Eigen::Index cols() const EIGEN_NOEXCEPT { return m_rhs.cols(); }
    
  protected:
    
    const DelassusOperatorDerived & m_lhs;
    const MatrixDerived & m_rhs;
  };

template<typename DelassusOperatorDerived>
struct DelassusOperatorBase
{
  typedef typename traits<DelassusOperatorDerived>::Scalar Scalar;
  typedef typename traits<DelassusOperatorDerived>::Vector Vector;
  typedef PowerIterationAlgoTpl<Vector> PowerIterationAlgo;

  DelassusOperatorDerived & derived() { return static_cast<DelassusOperatorDerived&>(*this); }
  const DelassusOperatorDerived & derived() const { return static_cast<const DelassusOperatorDerived&>(*this); }

  explicit DelassusOperatorBase()
  {}

  template<typename VectorLike>
  void updateDamping(const Eigen::MatrixBase<VectorLike> & vec)
  {
    derived().updateDamping(vec.derived());
  }

  void updateDamping(const Scalar mu)
  {
    derived().updateDamping(mu);
  }

  template<typename MatrixLike>
  void solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat) const
  {
    derived().solveInPlace(mat.const_cast_derived());
  }

  template<typename MatrixLike>
  typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLike)
  solve(const Eigen::MatrixBase<MatrixLike> & mat) const
  {
    return derived().solve(mat);
  }

  template<typename MatrixDerivedIn, typename MatrixDerivedOut>
  void solve(const Eigen::MatrixBase<MatrixDerivedIn> & x,
             const Eigen::MatrixBase<MatrixDerivedOut> & res) const
  {
    derived().solve(x.derived(), res.const_cast_derived());
  }

  template<typename MatrixIn, typename MatrixOut>
  void applyOnTheRight(const Eigen::MatrixBase<MatrixIn> & x,
                       const Eigen::MatrixBase<MatrixOut> & res) const
  {
    derived().applyOnTheRight(x.derived(), res.const_cast_derived());
  }

  template<typename MatrixDerived>
  typename MultiplicationOperatorReturnType<DelassusOperatorBase,Eigen::MatrixBase<MatrixDerived>>::type
  operator*(const Eigen::MatrixBase<MatrixDerived> & x) const
  {
    typedef typename MultiplicationOperatorReturnType<DelassusOperatorBase,Eigen::MatrixBase<MatrixDerived>>::type ReturnType;
    return ReturnType(derived(),x.derived());
  }

  Eigen::DenseIndex size() const { return derived().size(); }
  Eigen::DenseIndex rows() const { return derived().rows(); }
  Eigen::DenseIndex cols() const { return derived().cols(); }

}; // struct DelassusOperatorBase

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_delassus_operator_base_hpp__
