//
// Copyright (c) 2015-2024 CNRS INRIA
// Copyright (c) 2015 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#ifndef __pinocchio_python_fwd_hpp__
#define __pinocchio_python_fwd_hpp__

#include "pinocchio/bindings/python/context.hpp"
#include <eigenpy/eigenpy.hpp>

#include <memory>
#define PINOCCHIO_SHARED_PTR_HOLDER_TYPE(T) ::std::shared_ptr<T>

namespace pinocchio
{
  namespace python
  {
    using eigenpy::register_symbolic_link_to_registered_type;

    // Expose spatial classes
    void exposeSE3();
    void exposeForce();
    void exposeMotion();
    void exposeInertia();
    void exposeSymmetric3();
    void exposeExplog();
    void exposeSkew();
    void exposeLieGroups();

    // Expose math module
    void exposeRpy();
    void exposeEigenTypes();
    void exposeConversions();
    void exposeLinalg();
    void exposeTridiagonalMatrix();
    void exposeLanczosDecomposition();
    void exposeGramSchmidtOrthonormalisation();

    // Expose multibody classes
    void exposeJoints();
    void exposeModel();
    void exposeFrame();
    void exposeData();
    void exposeSampleModels();

    // Expose geometry module
    void exposeGeometry();

    // Expose parsers
    void exposeParsers();

    // Expose algorithms
    void exposeAlgorithms();
    void exposeExtras();
    void exposeConstraints();

#ifdef PINOCCHIO_PYTHON_INTERFACE_WITH_HPP_FCL_PYTHON_BINDINGS
    void exposeFCL();
    void exposeCollision();
#endif // PINOCCHIO_PYTHON_INTERFACE_WITH_HPP_FCL_PYTHON_BINDINGS

#ifdef PINOCCHIO_PYTHON_INTERFACE_WITH_OPENMP
    void exposePool();
    void exposeParallelAlgorithms();
#endif

#if defined(PINOCCHIO_PYTHON_INTERFACE_WITH_HPP_FCL_PYTHON_BINDINGS)                               \
  && defined(PINOCCHIO_PYTHON_INTERFACE_WITH_OPENMP)
    void exposeParallelCollision();
    void exposePoolCollision();
#endif

  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_fwd_hpp__
