//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_alloca_hpp__
#define __pinocchio_alloca_hpp__

#ifdef WIN32
#else
  #include <alloca.h>
#endif

#define PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, rows, cols)                                             \
  static_cast<Scalar *>(alloca(size_t(rows * cols) * sizeof(Scalar))), rows, cols

#endif // ifndef __pinocchio_alloca_hpp__
