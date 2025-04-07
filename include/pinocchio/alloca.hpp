//
// Copyright (c) 2024-2025 INRIA
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
#define PINOCCHIO_EIGEN_MAP_ALLOCA(S, rows, cols)                                                  \
  PINOCCHIO_EIGEN_MAP_ALLOCA_ALIGNED(S, rows, cols, EIGEN_DEFAULT_ALIGN_BYTES)
#define PINOCCHIO_EIGEN_MAP_ALLOCA_ALIGNED(S, rows, cols, align)                                   \
  static_cast<S *>(PINOCCHIO_ALIGNED_PTR(                                                          \
    PINOCCHIO_ALLOCA(size_t(rows * cols) * sizeof(S) + (align > 0 ? (align - 1) : 0)), align)),    \
    rows, cols

#endif // ifndef __pinocchio_alloca_hpp__
