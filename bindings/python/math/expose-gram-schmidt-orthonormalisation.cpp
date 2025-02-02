//
// Copyright (c) 2025 INRIA
//

#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/math/gram-schmidt-orthonormalisation.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    void exposeGramSchmidtOrthonormalisation()
    {
#ifndef PINOCCHIO_PYTHON_SKIP_CASADI_UNSUPPORTED
      const context::Scalar prec = Eigen::NumTraits<context::Scalar>::dummy_precision();

      bp::def(
        "orthogonalization",
        +[](
           const context::MatrixXs & mat, context::RefVectorXs vec,
           const context::Scalar threshold) -> void { orthogonalization(mat, vec, threshold); },
        (bp::args("basis", "vec"), bp::arg("threshold") = context::Scalar(0)),
        "Perform the Gram-Schmidt orthogonalization on the input/output vector for a given input "
        "basis.");

      bp::def(
        "orthonormalization",
        +[](context::MatrixXs & mat, const context::Scalar threshold) -> void {
          orthonormalization(mat, threshold);
        },
        (bp::arg("matrix"), bp::arg("threshold") = context::Scalar(0)),
        "Perform the orthonormalization of the input matrix via the Gram-Schmidt procedure.");

      bp::def(
        "isOrthonormal",
        +[](context::MatrixXs & matrix, const context::Scalar prec) -> bool {
          return isOrthonormal(matrix, prec);
        },
        (bp::arg("matrix"), bp::arg("prec") = prec),
        "Check whether the input matrix is orthonormal up to a given input precision.");
#endif // PINOCCHIO_PYTHON_SKIP_CASADI_UNSUPPORTED
    }

  } // namespace python
} // namespace pinocchio
