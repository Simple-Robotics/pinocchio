#
# Copyright (c) 2015-2018 CNRS
# Copyright (c) 2018-2025 INRIA
#

# ----------------------------------------------------
# --- BENCHMARK --------------------------------------
# ----------------------------------------------------
add_custom_target(${PROJECT_NAME}-benchmarks)

add_project_private_dependency(benchmark REQUIRED)

function(ADD_PINOCCHIO_BENCHMARK name)
  set(options
      PARSERS
      COLLISION
      PARALLEL
      COLLISION_PARALLEL
      CPPADCG
      CPPAD
      PARSERS_OPTIONAL
      COLLISION_OPTIONAL
      PARALLEL_OPTIONAL
      COLLISION_PARALLEL_OPTIONAL
      CPPADCG_OPTIONAL
      CPPAD_OPTIONAL)
  set(oneValueArgs)
  set(multiValueArgs)
  cmake_parse_arguments(bench "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(bench_name "pinocchio-benchmark-${name}")
  add_executable(${bench_name} ${name}.cpp)

  target_compile_definitions(${bench_name} PRIVATE PINOCCHIO_MODEL_DIR="${PINOCCHIO_MODEL_DIR}")

  target_link_libraries(${bench_name} PRIVATE pinocchio_default)
  target_link_libraries(${bench_name} PRIVATE benchmark::benchmark)

  set(cppad_special_define OFF)

  if(bench_PARSERS OR (bench_PARSERS_OPTIONAL AND TARGET pinocchio_parsers))
    target_link_libraries(${bench_name} PRIVATE pinocchio_parsers)
  endif()
  if(bench_COLLISION OR (bench_COLLISION_OPTIONAL AND TARGET pinocchio_collision))
    target_link_libraries(${bench_name} PRIVATE pinocchio_collision)
  endif()
  if(bench_PARALLEL OR (bench_PARALLEL_OPTIONAL AND TARGET pinocchio_parallel))
    target_link_libraries(${bench_name} PRIVATE pinocchio_parallel)
  endif()
  if(bench_COLLISION_PARALLEL
     OR (bench_COLLISION_PARALLEL_OPTIONAL AND TARGET pinocchio_collision_parallel))
    target_link_libraries(${bench_name} PRIVATE pinocchio_collision_parallel)
  endif()
  if(bench_CPPADCG OR (bench_CPPADCG_OPTIONAL AND TARGET pinocchio_cppadcg))
    target_link_libraries(${bench_name} PRIVATE pinocchio_cppadcg)
    set(cppad_special_define ON)
  endif()
  if(bench_CPPAD OR (bench_CPPAD_OPTIONAL AND TARGET pinocchio_cppad))
    target_link_libraries(${bench_name} PRIVATE pinocchio_cppad)
    set(cppad_special_define ON)
  endif()

  if(cppad_special_define)
    target_compile_definitions(${bench_name}
                               PRIVATE PINOCCHIO_CXX_COMPILER=\"${CMAKE_CXX_COMPILER}\")
    target_link_libraries(${bench_name} PRIVATE ${CMAKE_DL_LIBS})
  endif()

  add_dependencies(${PROJECT_NAME}-benchmarks ${bench_name})
endfunction()

add_pinocchio_benchmark(spatial-operations)

# timings
#
add_pinocchio_benchmark(timings PARSERS)
if(BUILD_WITH_CODEGEN_SUPPORT)
  add_pinocchio_benchmark(timings-cg PARSERS CPPADCG)
endif()

if(BUILD_WITH_OPENMP_SUPPORT)
  add_pinocchio_benchmark(timings-parallel PARSERS PARALLEL COLLISION_PARALLEL_OPTIONAL)
endif()

# timings cholesky
#
add_pinocchio_benchmark(timings-cholesky PARSERS)
add_pinocchio_benchmark(timings-loop-constrained-aba PARSERS)

# timings derivatives
#
add_pinocchio_benchmark(timings-derivatives PARSERS)
if(BUILD_WITH_AUTODIFF_SUPPORT)
  # timings-cppad-jit
  add_pinocchio_benchmark(timings-cppad-jit CPPAD)
endif()

# timings-eigen
add_pinocchio_benchmark(timings-eigen)

# timings-linalg
add_pinocchio_benchmark(timings-linalg-inverse)

# timings-geometry
#
if(BUILD_WITH_URDF_SUPPORT AND BUILD_WITH_COLLISION_SUPPORT)
  add_pinocchio_benchmark(timings-geometry PARSERS COLLISION)
endif()

# timings-jacobian
#
add_pinocchio_benchmark(timings-jacobian PARSERS)

# timings-contact-dynamics
#
add_pinocchio_benchmark(timings-delassus-operations PARSERS)
add_pinocchio_benchmark(timings-contact-dynamics PARSERS)
add_pinocchio_benchmark(timings-constrained-dynamics-derivatives PARSERS)

# timings-impulse-dynamics
#
add_pinocchio_benchmark(timings-impulse-dynamics PARSERS)
add_pinocchio_benchmark(timings-impulse-dynamics-derivatives PARSERS)
