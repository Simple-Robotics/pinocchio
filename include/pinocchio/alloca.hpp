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

#define PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, rows, cols)                                             \
  static_cast<Scalar *>(alloca(size_t(rows * cols) * sizeof(Scalar))), rows, cols
#define PINOCCHIO_ALLOCA EIGEN_ALLOCA

#endif // ifndef __pinocchio_alloca_hpp__
