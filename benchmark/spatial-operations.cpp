//
// Copyright (c) 2025 INRIA
//

#include "pinocchio/spatial/inertia.hpp"
#include "pinocchio/spatial/se3.hpp"
#include "pinocchio/spatial/motion.hpp"
#include "pinocchio/spatial/force.hpp"

#include "pinocchio/algorithm/aba.hpp"

#include <benchmark/benchmark.h>

using namespace pinocchio;

static void CustomArguments(benchmark::internal::Benchmark * b)
{
  b->MinWarmUpTime(30.)->Iterations(int(1e8));
}

void homgeneous_multiplication(benchmark::State & st)
{
  using Matrix4 = SE3::Matrix4;
  Matrix4 M1 = Matrix4::Random();
  Matrix4 M2 = Matrix4::Random();
  for (auto _ : st)
  {
    Matrix4 M3 = M1 * M2;
    benchmark::DoNotOptimize(M3);
  }
}

template<typename SpatialType>
void se3_action(benchmark::State & st)
{
  SE3 M = SE3::Random();
  SpatialType obj = SpatialType::Random();
  for (auto _ : st)
  {
    auto res = M.act(obj);
    benchmark::DoNotOptimize(res);
  }
}

void cross_operation(benchmark::State & st)
{
  SE3::Vector3 t = SE3::Vector3::Random();
  for (auto _ : st)
  {
    auto t_skew = skew(t);
    benchmark::DoNotOptimize(t_skew);
  }
}

template<int Size>
void matrix_matrix_add(benchmark::State & st)
{
  using Matrix = Eigen::Matrix<double, Size, Size>;
  Matrix M1 = Matrix::Random();
  Matrix M2 = Matrix::Random();
  for (auto _ : st)
  {
    Matrix M3 = M1 + M2;
    benchmark::DoNotOptimize(M3);
  }
}

template<int Size>
void matrix_matrix_plus_equal(benchmark::State & st)
{
  using Matrix = Eigen::Matrix<double, Size, Size>;
  Matrix M1 = Matrix::Random();
  Matrix M2 = Matrix::Random();
  for (auto _ : st)
  {
    M2 += M1;
    benchmark::DoNotOptimize(M2);
  }
}

void inertia_matrix_transformation_ouput_arg(benchmark::State & st)
{
  const Inertia spatial_inertia = Inertia::Random();
  Inertia::Matrix6 matrix;
  for (auto _ : st)
  {
    spatial_inertia.matrix(matrix);
    benchmark::DoNotOptimize(matrix);
  }
}

void inertia_matrix_transformation_assignment(benchmark::State & st)
{
  const Inertia spatial_inertia = Inertia::Random();
  for (auto _ : st)
  {
    auto matrix = spatial_inertia.matrix();
    benchmark::DoNotOptimize(matrix);
  }
}

void inertia_se3_action(benchmark::State & st)
{
  Inertia spatial_inertia = Inertia::Random();
  SE3 placement = SE3::Random();
  for (auto _ : st)
  {
    auto inertia = placement.act(spatial_inertia);
    benchmark::DoNotOptimize(inertia);
  }
}

void dense_inertia_se3_action(benchmark::State & st)
{
  using Matrix6 = Inertia::Matrix6;
  using Scalar = Inertia::Scalar;
  Inertia spatial_inertia = Inertia::Random();
  Matrix6 dense_inertia_matrix = spatial_inertia.matrix();
  SE3 placement = SE3::Random();
  for (auto _ : st)
  {
    Matrix6 res = impl::internal::SE3actOn<Scalar>::run(placement, dense_inertia_matrix);
    benchmark::DoNotOptimize(res);
  }
}

void inertia_set_zero(benchmark::State & st)
{
  for (auto _ : st)
  {
    Inertia spatial_inertia;
    spatial_inertia.setZero();
    benchmark::DoNotOptimize(spatial_inertia);
  }
}

void inertia_init_from_zero(benchmark::State & st)
{
  for (auto _ : st)
  {
    Inertia spatial_inertia = Inertia::Zero();
    benchmark::DoNotOptimize(spatial_inertia);
  }
}

void inertia_motion_product(benchmark::State & st)
{
  Inertia spatial_inertia = Inertia::Random();
  Motion spatial_motion = Motion::Random();
  for (auto _ : st)
  {
    Force spatial_force = spatial_inertia * spatial_motion;
    benchmark::DoNotOptimize(spatial_force);
  }
}

void force_init_from_zero(benchmark::State & st)
{
  for (auto _ : st)
  {
    Force spatial_force = Force::Zero();
    benchmark::DoNotOptimize(spatial_force);
  }
}

void force_set_zero(benchmark::State & st)
{
  for (auto _ : st)
  {
    Force spatial_force;
    spatial_force.setZero();
    benchmark::DoNotOptimize(spatial_force);
  }
}

void force_no_init(benchmark::State & st)
{
  for (auto _ : st)
  {
    Force spatial_force;
    benchmark::DoNotOptimize(spatial_force);
  }
}

void matrix_times_vector6(benchmark::State & st)
{
  Inertia::Matrix6 mat6 = Inertia::Matrix6::Random();
  Motion::Vector6 vec6 = Motion::Vector6::Random();
  Force res = Force::Zero();
  for (auto _ : st)
  {
    res.toVector().noalias() += mat6 * vec6;
    benchmark::DoNotOptimize(res);
  }
}

void matrix_times_vector_static_dispatch(benchmark::State & st)
{
  typedef Eigen::Matrix<double, 6, 6> Matrix6;
  typedef Eigen::Matrix<double, 6, 1> Vector6;

  Matrix6 mat = Matrix6::Random();
  Vector6 vec = Vector6::Random();
  Force res = Force::Zero();
  for (auto _ : st)
  {
    auto size = st.range(0);
#define CASE_OP(n)                                                                                 \
  case n:                                                                                          \
    res.toVector().noalias() += mat.leftCols<n>() * vec.head<n>();                                 \
    break;

#define CASE_OP_DYN(n)                                                                             \
  case n:                                                                                          \
    res.toVector().noalias() += mat.leftCols(size) * vec.head(size);                               \
    break;

    switch (size)
    {
      CASE_OP(1)
      CASE_OP(2)
      CASE_OP(3)
      CASE_OP(4)
      CASE_OP(5)
      CASE_OP(6)
    default:
      break;
    }

#undef CASE_OP
    benchmark::DoNotOptimize(size);
    benchmark::DoNotOptimize(res);
  }
}

void matrix_times_vector_dynamic_dispatch(benchmark::State & st)
{
  typedef Eigen::Matrix<double, 6, 6> Matrix6;
  typedef Eigen::Matrix<double, 6, 1> Vector6;

  Matrix6 mat = Matrix6::Random();
  Vector6 vec = Vector6::Random();
  Force res = Force::Zero();
  for (auto _ : st)
  {
    auto size = st.range(0);

    res.toVector().noalias() += mat.leftCols(size) * vec.head(size);
    benchmark::DoNotOptimize(size);
    benchmark::DoNotOptimize(res);
  }
}

// Matrix operations
BENCHMARK(matrix_matrix_add<2>)->Apply(CustomArguments);
BENCHMARK(matrix_matrix_add<3>)->Apply(CustomArguments);
BENCHMARK(matrix_matrix_add<4>)->Apply(CustomArguments);
BENCHMARK(matrix_matrix_add<5>)->Apply(CustomArguments);
BENCHMARK(matrix_matrix_add<6>)->Apply(CustomArguments);

BENCHMARK(matrix_matrix_plus_equal<2>)->Apply(CustomArguments);
BENCHMARK(matrix_matrix_plus_equal<3>)->Apply(CustomArguments);
BENCHMARK(matrix_matrix_plus_equal<4>)->Apply(CustomArguments);
BENCHMARK(matrix_matrix_plus_equal<5>)->Apply(CustomArguments);
BENCHMARK(matrix_matrix_plus_equal<6>)->Apply(CustomArguments);

// SE3
BENCHMARK(homgeneous_multiplication)->Apply(CustomArguments);
BENCHMARK(se3_action<SE3>)->Apply(CustomArguments);
BENCHMARK(se3_action<Force>)->Apply(CustomArguments);
BENCHMARK(se3_action<Motion>)->Apply(CustomArguments);
BENCHMARK(se3_action<Inertia>)->Apply(CustomArguments);
BENCHMARK(cross_operation)->Apply(CustomArguments);

// Inertias
BENCHMARK(inertia_init_from_zero)->Apply(CustomArguments);
BENCHMARK(inertia_set_zero)->Apply(CustomArguments);
BENCHMARK(inertia_matrix_transformation_ouput_arg)->Apply(CustomArguments);
BENCHMARK(inertia_matrix_transformation_assignment)->Apply(CustomArguments);
BENCHMARK(inertia_se3_action)->Apply(CustomArguments);
BENCHMARK(inertia_motion_product)->Apply(CustomArguments);
BENCHMARK(dense_inertia_se3_action)->Apply(CustomArguments);

// Forces
BENCHMARK(force_no_init)->Apply(CustomArguments);
BENCHMARK(force_init_from_zero)->Apply(CustomArguments);
BENCHMARK(force_set_zero)->Apply(CustomArguments);

// Others
BENCHMARK(matrix_times_vector6)->Apply(CustomArguments);
BENCHMARK(matrix_times_vector_static_dispatch)
  ->Apply(CustomArguments)
  ->Arg(1)
  ->Arg(2)
  ->Arg(3)
  ->Arg(4)
  ->Arg(5)
  ->Arg(6);
BENCHMARK(matrix_times_vector_dynamic_dispatch)
  ->Apply(CustomArguments)
  ->Arg(1)
  ->Arg(2)
  ->Arg(3)
  ->Arg(4)
  ->Arg(5)
  ->Arg(6);
BENCHMARK_MAIN();
