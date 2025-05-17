//
// Copyright (c) 2025 INRIA
//

#include "pinocchio/fwd.hpp"
#include "pinocchio/multibody/joint/joint-common-operations.hpp"

#include <benchmark/benchmark.h>

using namespace pinocchio;

template<int Size, int Options = 0>
using Matrix = Eigen::Matrix<double, Size, Size, Options>;
using DynamicMatrix = Matrix<Eigen::Dynamic>;

#define DEFINE_MATRIX(size)                                                                        \
  using Matrix##size = Matrix<size>;                                                               \
  using RowMatrix##size = Matrix<size, Eigen::RowMajor>;

DEFINE_MATRIX(1)
DEFINE_MATRIX(2)
DEFINE_MATRIX(3)
DEFINE_MATRIX(4)
DEFINE_MATRIX(5)
DEFINE_MATRIX(6)
DEFINE_MATRIX(7)
DEFINE_MATRIX(8)
DEFINE_MATRIX(9)
DEFINE_MATRIX(10)
DEFINE_MATRIX(11)
DEFINE_MATRIX(12)

static void CustomArgumentsStaticMatrix(benchmark::internal::Benchmark * b)
{
  b->MinWarmUpTime(3.);
}

struct MatrixInverseEigen
{
  template<typename M1, typename M2>
  PINOCCHIO_DONT_INLINE static void
  run(const Eigen::MatrixBase<M1> & mat, const Eigen::MatrixBase<M2> & mat_inv)
  {
    mat_inv.const_cast_derived().noalias() = mat.inverse();
  }
};

struct MatrixInversePartialPivLU
{
  template<typename M1, typename M2>
  PINOCCHIO_DONT_INLINE static void
  run(const Eigen::MatrixBase<M1> & mat, const Eigen::MatrixBase<M2> & mat_inv)
  {
    mat_inv.const_cast_derived().noalias() = mat.partialPivLu().inverse();
  }
};

struct MatrixInverseLLT
{
  template<typename M1, typename M2>
  PINOCCHIO_DONT_INLINE static void
  run(const Eigen::MatrixBase<M1> & mat, const Eigen::MatrixBase<M2> & mat_inv)
  {
    mat_inv.const_cast_derived().setIdentity();
    mat.llt().solveInPlace(mat_inv.const_cast_derived());
  }
};

struct MatrixInverseLDLT
{
  template<typename M1, typename M2>
  PINOCCHIO_DONT_INLINE static void
  run(const Eigen::MatrixBase<M1> & mat, const Eigen::MatrixBase<M2> & mat_inv)
  {
    mat_inv.const_cast_derived().setIdentity();
    mat.ldlt().solveInPlace(mat_inv.const_cast_derived());
  }
};

struct MatrixInversePinocchio
{
  template<typename M1, typename M2>
  PINOCCHIO_DONT_INLINE static void
  run(const Eigen::MatrixBase<M1> & mat, const Eigen::MatrixBase<M2> & mat_inv)
  {
    ::pinocchio::matrix_inversion(mat, mat_inv.const_cast_derived());
  }
};

struct MatrixInverseCodeGenerated
{
  template<typename M1, typename M2>
  PINOCCHIO_DONT_INLINE static void
  run(const Eigen::MatrixBase<M1> & mat, const Eigen::MatrixBase<M2> & mat_inv)
  {
    ::pinocchio::matrix_inversion_code_generated(mat, mat_inv.const_cast_derived());
  }
};

template<typename InputMatrix, typename OutputMatrix, class MatrixInverseFunctor>
static void static_matrix_inversion_call(benchmark::State & st)
{
  const InputMatrix input_matrix = InputMatrix::Identity();
  OutputMatrix res = OutputMatrix::Zero(input_matrix.rows(), input_matrix.cols());
  for (auto _ : st)
  {
    MatrixInverseFunctor::run(input_matrix, res);
    benchmark::DoNotOptimize(res);
  }
}

static void CustomArgumentsDynamicMatrix(benchmark::internal::Benchmark * b)
{
  b->MinWarmUpTime(3.);
  for (int size = 1; size <= 12; ++size)
    b->Arg(size);

  b->Arg(20)->Arg(50)->Arg(100);
}

template<typename InputMatrix, typename OutputMatrix, class MatrixInverseFunctor>
static void dynamic_matrix_inversion_call(benchmark::State & st)
{
  const auto size = st.range(0);
  const InputMatrix input_matrix = InputMatrix::Identity(size, size);
  OutputMatrix res = OutputMatrix::Zero(input_matrix.rows(), input_matrix.cols());
  for (auto _ : st)
  {
    MatrixInverseFunctor::run(input_matrix, res);
    // pinocchio::internal::MatrixInversionDynamicMatrixImpl::run(input_matrix, res);
    benchmark::DoNotOptimize(res);
  }
}

template<typename Scalar>
PINOCCHIO_DONT_INLINE void scalar_inversion_op(const Scalar & input_scalar, Scalar & output)
{
  output = Scalar(1) / input_scalar;
}

void scalar_inversion(benchmark::State & st)
{
  const double input_scalar = Matrix1::Random().coeff(0, 0);
  double scalar_inv = 0.;
  for (auto _ : st)
  {
    scalar_inversion_op(input_scalar, scalar_inv);
    benchmark::DoNotOptimize(scalar_inv);
  }
}

template<typename Scalar>
PINOCCHIO_DONT_INLINE void scalar_sqrt_op(const Scalar & input_scalar, Scalar & output)
{
  output = math::sqrt(input_scalar);
}

void scalar_sqrt(benchmark::State & st)
{
  const double input_scalar = Matrix1::Random().coeff(0, 0);
  double res_scalar = 0.;
  for (auto _ : st)
  {
    scalar_sqrt_op(input_scalar, res_scalar);
    benchmark::DoNotOptimize(res_scalar);
  }
}

template<typename Scalar>
PINOCCHIO_DONT_INLINE void scalar_multiplication_op(const Scalar & input_scalar, Scalar & output)
{
  output = input_scalar * input_scalar;
}

void scalar_multiplication(benchmark::State & st)
{
  const double input_scalar = Matrix1::Random().coeff(0, 0);
  double res_scalar = 0.;
  for (auto _ : st)
  {
    scalar_multiplication_op(input_scalar, res_scalar);
    benchmark::DoNotOptimize(res_scalar);
  }
}

#define BENCH_MATRIX_INVERSION(Call, Type, MatrixInverseFunctor, Arg)                              \
  BENCHMARK(Call<Type, Type, MatrixInverseFunctor>)->Apply(Arg);                                   \
  //BENCHMARK(Call<Row##Type,Type,MatrixInverseEigen>); \
//BENCHMARK(Call<Type,Row##Type,MatrixInverseEigen>); \
//BENCHMARK(Call<Row##Type,Row##Type,MatrixInverseEigen>);

#define BENCH_STATIC_MATRIX_INVERSION(Type, MatrixInverseFunctor)                                  \
  BENCH_MATRIX_INVERSION(                                                                          \
    static_matrix_inversion_call, Type, MatrixInverseFunctor, CustomArgumentsStaticMatrix)

#define BENCH_STATIC_MATRIX_INVERSION_ALL(MatrixInverseFunctor)                                    \
  BENCH_STATIC_MATRIX_INVERSION(Matrix1, MatrixInverseFunctor)                                     \
  BENCH_STATIC_MATRIX_INVERSION(Matrix2, MatrixInverseFunctor)                                     \
  BENCH_STATIC_MATRIX_INVERSION(Matrix3, MatrixInverseFunctor)                                     \
  BENCH_STATIC_MATRIX_INVERSION(Matrix4, MatrixInverseFunctor)                                     \
  BENCH_STATIC_MATRIX_INVERSION(Matrix5, MatrixInverseFunctor)                                     \
  BENCH_STATIC_MATRIX_INVERSION(Matrix6, MatrixInverseFunctor)                                     \
  BENCH_STATIC_MATRIX_INVERSION(Matrix7, MatrixInverseFunctor)                                     \
  BENCH_STATIC_MATRIX_INVERSION(Matrix8, MatrixInverseFunctor)                                     \
  BENCH_STATIC_MATRIX_INVERSION(Matrix9, MatrixInverseFunctor)                                     \
  BENCH_STATIC_MATRIX_INVERSION(Matrix10, MatrixInverseFunctor)                                    \
  BENCH_STATIC_MATRIX_INVERSION(Matrix11, MatrixInverseFunctor)                                    \
  BENCH_STATIC_MATRIX_INVERSION(Matrix12, MatrixInverseFunctor)

#define BENCH_DYNAMIC_MATRIX_INVERSION_ALL(MatrixInverseFunctor)                                   \
  BENCH_MATRIX_INVERSION(                                                                          \
    dynamic_matrix_inversion_call, DynamicMatrix, MatrixInverseFunctor,                            \
    CustomArgumentsDynamicMatrix)

BENCHMARK(scalar_inversion)->Apply(CustomArgumentsStaticMatrix);
BENCHMARK(scalar_sqrt)->Apply(CustomArgumentsStaticMatrix);
BENCHMARK(scalar_multiplication)->Apply(CustomArgumentsStaticMatrix);

BENCH_STATIC_MATRIX_INVERSION_ALL(MatrixInverseEigen)
BENCH_STATIC_MATRIX_INVERSION_ALL(MatrixInversePartialPivLU)
BENCH_STATIC_MATRIX_INVERSION_ALL(MatrixInverseLLT)
BENCH_STATIC_MATRIX_INVERSION_ALL(MatrixInverseLDLT)
BENCH_STATIC_MATRIX_INVERSION_ALL(MatrixInverseCodeGenerated)
BENCH_STATIC_MATRIX_INVERSION_ALL(MatrixInversePinocchio)

BENCH_DYNAMIC_MATRIX_INVERSION_ALL(MatrixInverseEigen)
BENCH_DYNAMIC_MATRIX_INVERSION_ALL(MatrixInversePinocchio)

BENCHMARK_MAIN();
