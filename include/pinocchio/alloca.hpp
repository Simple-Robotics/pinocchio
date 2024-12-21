//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_alloca_hpp__
#define __pinocchio_alloca_hpp__

#ifdef WIN32
  #include <malloc.h>
#else
  #include <alloca.h>
#endif

#define PINOCCHIO_ALLOCA EIGEN_ALLOCA
#define PINOCCHIO_ALIGNED_PTR(ptr, align)                                                          \
  reinterpret_cast<void *>(((intptr_t)ptr + (align - 1)) & ~(align - 1))
#define PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, rows, cols)                                             \
  PINOCCHIO_EIGEN_MAP_ALLOCA_ALIGNED(Scalar, rows, cols, EIGEN_DEFAULT_ALIGN_BYTES)
#define PINOCCHIO_EIGEN_MAP_ALLOCA_ALIGNED(Scalar, rows, cols, align)                              \
  static_cast<Scalar *>(PINOCCHIO_ALIGNED_PTR(                                                     \
    PINOCCHIO_ALLOCA(size_t(rows * cols) * sizeof(Scalar) + (align > 0 ? (align - 1) : 1)),        \
    align)),                                                                                       \
    rows, cols

#endif // ifndef __pinocchio_alloca_hpp__
