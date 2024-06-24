//
// Copyright (c) 2015-2022 CNRS INRIA
//

#ifndef __pinocchio_multibody_geometry_object_hxx__
#define __pinocchio_multibody_geometry_object_hxx__

#include <limits>

namespace pinocchio
{

  inline std::ostream & operator<<(std::ostream & os, const GeometryObject & geom_object)
  {
    os << "Name: \t \n"
       << geom_object.name << "\n"
       << "Parent frame ID: \t \n"
       << geom_object.parentFrame << "\n"
       << "Parent joint ID: \t \n"
       << geom_object.parentJoint << "\n"
       << "Position in parent frame: \t \n"
       << geom_object.placement << "\n"
       << "Absolute path to mesh file: \t \n"
       << geom_object.meshPath << "\n"
       << "Scale for transformation of the mesh: \t \n"
       << geom_object.meshScale.transpose() << "\n"
       << "Disable collision: \t \n"
       << geom_object.disableCollision << "\n"
       << std::endl;
    return os;
  }

  inline FrictionCoefficientMatrix::FrictionCoefficientMatrix()
  {
    // Initialize the matrix with zeros, in case we forget to set some coefficients.
    for (Index i = 0; i < friction_coefficient_matrix.rows(); ++i)
    {
      for (Index j = 0; j < friction_coefficient_matrix.cols(); ++j)
      {
        friction_coefficient_matrix.coeffRef(i, j) = 0.0;
      }
    }

    // Sources: https://en.wikipedia.org/wiki/Friction
    //          https://www.engineeringtoolbox.com/friction-coefficients-d_778.html
    // These are really rough estimates, and should be replaced by more accurate values if
    // available.
    friction_coefficient_matrix.coeffRef(METAL, METAL) = 0.75;
    friction_coefficient_matrix.coeffRef(METAL, WOOD) = 0.5;
    friction_coefficient_matrix.coeffRef(METAL, PLASTIC) = 0.2;
    friction_coefficient_matrix.coeffRef(METAL, ICE) = 0.03;
    friction_coefficient_matrix.coeffRef(METAL, CONCRETE) = 0.85;

    friction_coefficient_matrix.coeffRef(WOOD, WOOD) = 0.4;
    friction_coefficient_matrix.coeffRef(WOOD, PLASTIC) = 0.3;
    friction_coefficient_matrix.coeffRef(WOOD, ICE) = 0.03;
    friction_coefficient_matrix.coeffRef(WOOD, CONCRETE) = 0.65;

    friction_coefficient_matrix.coeffRef(PLASTIC, PLASTIC) = 0.2;
    friction_coefficient_matrix.coeffRef(PLASTIC, ICE) = 0.02;
    friction_coefficient_matrix.coeffRef(PLASTIC, CONCRETE) = 0.55;

    friction_coefficient_matrix.coeffRef(ICE, ICE) = 0.01;
    friction_coefficient_matrix.coeffRef(ICE, CONCRETE) = 0.25;

    friction_coefficient_matrix.coeffRef(CONCRETE, CONCRETE) = 1.0;

    // Symmetrize the matrix
    friction_coefficient_matrix.triangularView<Eigen::StrictlyLower>() =
      friction_coefficient_matrix.transpose();
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_multibody_geometry_object_hxx__
