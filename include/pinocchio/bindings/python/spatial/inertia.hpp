//
// Copyright (c) 2015-2024 CNRS INRIA
// Copyright (c) 2016 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#ifndef __pinocchio_python_spatial_inertia_hpp__
#define __pinocchio_python_spatial_inertia_hpp__

#include <eigenpy/exception.hpp>
#include <eigenpy/eigenpy.hpp>
#include <eigenpy/memory.hpp>
#include <boost/python/tuple.hpp>

#include "pinocchio/spatial/inertia.hpp"

#include "pinocchio/bindings/python/utils/cast.hpp"
#include "pinocchio/bindings/python/utils/copyable.hpp"
#include "pinocchio/bindings/python/utils/printable.hpp"

#if EIGENPY_VERSION_AT_MOST(2, 8, 1)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(pinocchio::Inertia)
#endif

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<typename Inertia>
    struct InertiaPythonVisitor : public boost::python::def_visitor<InertiaPythonVisitor<Inertia>>
    {
      enum
      {
        Options = Inertia::Options
      };
      typedef typename Inertia::Scalar Scalar;
      typedef typename Inertia::Vector3 Vector3;
      typedef typename Inertia::Matrix3 Matrix3;
      typedef typename Inertia::Vector6 Vector6;
      typedef typename Inertia::Matrix6 Matrix6;
      typedef typename Inertia::Matrix4 Matrix4;

      typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
      typedef MotionTpl<Scalar, Options> Motion;
      typedef ForceTpl<Scalar, Options> Force;

    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        static const Scalar dummy_precision = Eigen::NumTraits<Scalar>::dummy_precision();
        PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
        PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_SELF_ASSIGN_OVERLOADED
        cl.def(
            "__init__",
            bp::make_constructor(
              &InertiaPythonVisitor::makeFromMCI, bp::default_call_policies(),
              bp::args("mass", "lever", "inertia")),
            "Initialize from mass, lever and 3d inertia.")
          .def(bp::init<>(bp::arg("self"), "Default constructor."))
          .def(bp::init<const Inertia &>((bp::arg("self"), bp::arg("clone")), "Copy constructor"))

          .add_property(
            "mass", &InertiaPythonVisitor::getMass, &InertiaPythonVisitor::setMass,
            "Mass of the Spatial Inertia.")
          .add_property(
            "lever",
            bp::make_function(
              (typename Inertia::Vector3 & (Inertia::*)()) & Inertia::lever,
              bp::return_internal_reference<>()),
            &InertiaPythonVisitor::setLever,
            "Center of mass location of the Spatial Inertia. It corresponds to the location of the "
            "center of mass regarding to the frame where the Spatial Inertia is expressed.")
          .add_property(
            "inertia", &InertiaPythonVisitor::getInertia, &InertiaPythonVisitor::setInertia,
            "Rotational part of the Spatial Inertia, i.e. a symmetric matrix "
            "representing the rotational inertia around the center of mass.")

          .def("matrix", (Matrix6(Inertia::*)() const) & Inertia::matrix, bp::arg("self"))
          .def("inverse", (Matrix6(Inertia::*)() const) & Inertia::inverse, bp::arg("self"))
          .def(
            "se3Action", &Inertia::template se3Action<Scalar, Options>, bp::args("self", "M"),
            "Returns the result of the action of M on *this.")
          .def(
            "se3ActionInverse", &Inertia::template se3ActionInverse<Scalar, Options>,
            bp::args("self", "M"), "Returns the result of the action of the inverse of M on *this.")

          .def(
            "setIdentity", &Inertia::setIdentity, bp::arg("self"),
            "Set *this to be the Identity inertia.")
          .def(
            "setZero", &Inertia::setZero, bp::arg("self"),
            "Set all the components of *this to zero.")
          .def(
            "setRandom", &Inertia::setRandom, bp::arg("self"),
            "Set all the components of *this to random values.")

          .def(bp::self + bp::self)
          .def(bp::self += bp::self)
          .def(bp::self - bp::self)
          .def(bp::self -= bp::self)
          .def(bp::self * bp::other<Motion>())

          .add_property("np", (Matrix6(Inertia::*)() const) & Inertia::matrix)
          .def(
            "vxiv", &Inertia::template vxiv<Motion>, bp::args("self", "v"),
            "Returns the result of v x Iv.")
          .def(
            "vtiv", &Inertia::template vtiv<Motion>, bp::args("self", "v"),
            "Returns the result of v.T * Iv.")
          .def(
            "vxi",
            (Matrix6(Inertia::*)(const MotionDense<Motion> &) const)
              & Inertia::template vxi<Motion>,
            bp::args("self", "v"), "Returns the result of v x* I, a 6x6 matrix.")
          .def(
            "ivx",
            (Matrix6(Inertia::*)(const MotionDense<Motion> &) const)
              & Inertia::template ivx<Motion>,
            bp::args("self", "v"), "Returns the result of I vx, a 6x6 matrix.")
          .def(
            "variation",
            (Matrix6(Inertia::*)(const MotionDense<Motion> &) const)
              & Inertia::template variation<Motion>,
            bp::args("self", "v"), "Returns the time derivative of the inertia.")

#ifndef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS
          .def(bp::self == bp::self)
          .def(bp::self != bp::self)
#endif

#ifndef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS
          .def(
            "isApprox", &Inertia::isApprox,
            (bp::arg("self"), bp::arg("other"), bp::arg("prec") = dummy_precision),
            "Returns true if *this is approximately equal to other, within the precision given "
            "by prec.")

          .def(
            "isZero", &Inertia::isZero, (bp::arg("self"), bp::arg("prec") = dummy_precision),
            "Returns true if *this is approximately equal to the zero Inertia, within the "
            "precision given by prec.")
#endif

          .def("Identity", &Inertia::Identity, "Returns the identity Inertia.")
          .staticmethod("Identity")
          .def("Zero", &Inertia::Zero, "Returns the zero Inertia.")
          .staticmethod("Zero")
          .def("Random", &Inertia::Random, "Returns a random Inertia.")
          .staticmethod("Random")

          .def(
            "toDynamicParameters", &InertiaPythonVisitor::toDynamicParameters_proxy,
            bp::arg("self"),
            "Returns the representation of the matrix as a vector of dynamic parameters."
            "\nThe parameters are given as v = [m, mc_x, mc_y, mc_z, I_{xx}, I_{xy}, I_{yy}, "
            "I_{xz}, I_{yz}, I_{zz}]^T "
            "where I = I_C + mS^T(c)S(c) and I_C has its origin at the barycenter")
          .def(
            "FromDynamicParameters", &Inertia::template FromDynamicParameters<VectorXs>,
            bp::args("dynamic_parameters"),
            "Builds and inertia matrix from a vector of dynamic parameters."
            "\nThe parameters are given as dynamic_parameters = [m, mc_x, mc_y, mc_z, I_{xx}, "
            "I_{xy}, I_{yy}, I_{xz}, I_{yz}, I_{zz}]^T "
            "where I = I_C + mS^T(c)S(c) and I_C has its origin at the barycenter.")
          .staticmethod("FromDynamicParameters")
          .def(
            "toPseudoInertia", &InertiaPythonVisitor::toPseudoInertia_proxy, 
            bp::arg("self"),
            "Converts the inertia to a pseudo inertia matrix."
            "\nThe returned 4x4 pseudo inertia matrix has the form:"
            "\n[[ -0.5*I_xx + 0.5*I_yy + 0.5*I_zz, -I_xy, -I_xz, mr_x],"
            "\n [ -I_xy, 0.5*I_xx - 0.5*I_yy + 0.5*I_zz, -I_yz, mr_y],"
            "\n [ -I_xz, -I_yz, 0.5*I_xx + 0.5*I_yy - 0.5*I_zz, mr_z],"
            "\n [ mr_x, mr_y, mr_z, m ]].")
          .def(
            "FromPseudoInertia", &Inertia::FromPseudoInertia, 
            bp::args("pseudo_inertia"),
            "Builds an inertia matrix from a 4x4 pseudo inertia matrix."
            "\nThe parameters are given as"
            "\npseudo_inertia = [[ -0.5*I_xx + 0.5*I_yy + 0.5*I_zz, -I_xy, -I_xz, mr_x],"
            "\n [ -I_xy, 0.5*I_xx - 0.5*I_yy + 0.5*I_zz, -I_yz, mr_y],"
            "\n [ -I_xz, -I_yz, 0.5*I_xx + 0.5*I_yy - 0.5*I_zz, mr_z],"
            "\n [ mr_x, mr_y, mr_z, m ]].")
          .staticmethod("FromPseudoInertia")
          .def(
            "LogcholToDynamicParameters", &Inertia::template LogcholToDynamicParameters<VectorXs>,
            bp::args("log_cholesky"),
            "Converts logarithmic Cholesky parameters directly to theta parameters."
            "\nThe parameters are given as log_cholesky = [alpha, d_1, d_2, d_3, s_{12}, s_{23}, "
            "s_{13}, t_1, t_2, t_3]."
            "\nThe returned vector contains the dynamic parameters.")
          .staticmethod("LogcholToDynamicParameters")
          .def(
            "FromLogCholeskyParameters", &Inertia::template FromLogCholeskyParameters<VectorXs>,
            bp::args("log_cholesky"),
            "Builds an InertiaTpl from log Cholesky parameters."
            "\nThe parameters are given as log_cholesky = [alpha, d_1, d_2, d_3, s_{12}, s_{23}, "
            "s_{13}, t_1, t_2, t_3] "
            "\nThe returned InertiaTpl object is constructed from the provided log Cholesky "
            "parameters.")
          .staticmethod("FromLogCholeskyParameters")
          .def(
            "FromSphere", &Inertia::FromSphere, bp::args("mass", "radius"),
            "Returns the Inertia of a sphere defined by a given mass and radius.")
          .staticmethod("FromSphere")
          .def(
            "FromEllipsoid", &Inertia::FromEllipsoid,
            bp::args("mass", "length_x", "length_y", "length_z"),
            "Returns the Inertia of an ellipsoid shape defined by a mass and given dimensions "
            "the semi-axis of values length_{x,y,z}.")
          .staticmethod("FromEllipsoid")
          .def(
            "FromCylinder", &Inertia::FromCylinder, bp::args("mass", "radius", "length"),
            "Returns the Inertia of a cylinder defined by its mass, radius and length along the "
            "Z axis.")
          .staticmethod("FromCylinder")
          .def(
            "FromBox", &Inertia::FromBox, bp::args("mass", "length_x", "length_y", "length_z"),
            "Returns the Inertia of a box shape with a mass and of dimension the semi axis of "
            "length_{x,y,z}.")
          .staticmethod("FromBox")
          .def(
            "FromCapsule", &Inertia::FromCapsule, bp::args("mass", "radius", "height"),
            "Computes the Inertia of a capsule defined by its mass, radius and length along the "
            "Z axis. Assumes a uniform density.")
          .staticmethod("FromCapsule")

          .def("__array__", (Matrix6(Inertia::*)() const) & Inertia::matrix)
          .def(
            "__array__", &__array__,
            (bp::arg("self"), bp::arg("dtype") = bp::object(), bp::arg("copy") = bp::object()))
#ifndef PINOCCHIO_PYTHON_NO_SERIALIZATION
          .def_pickle(Pickle())
#endif
          ;
        PINOCCHIO_COMPILER_DIAGNOSTIC_POP
      }

      static Scalar getMass(const Inertia & self)
      {
        return self.mass();
      }
      static void setMass(Inertia & self, Scalar mass)
      {
        self.mass() = mass;
      }

      static void setLever(Inertia & self, const Vector3 & lever)
      {
        self.lever() = lever;
      }

      static Matrix3 getInertia(const Inertia & self)
      {
        return self.inertia().matrix();
      }
      //      static void setInertia(Inertia & self, const Vector6 & minimal_inertia) {
      //      self.inertia().data() = minimal_inertia; }
      static void setInertia(Inertia & self, const Matrix3 & symmetric_inertia)
      {
        if (!check_expression_if_real<Scalar>(
              isZero(symmetric_inertia - symmetric_inertia.transpose())))
          throw eigenpy::Exception("The 3d inertia should be symmetric.");
        self.inertia().data() << symmetric_inertia(0, 0), symmetric_inertia(1, 0),
          symmetric_inertia(1, 1), symmetric_inertia(0, 2), symmetric_inertia(1, 2),
          symmetric_inertia(2, 2);
      }

      static VectorXs toDynamicParameters_proxy(const Inertia & self)
      {
        return self.toDynamicParameters();
      }

      static Matrix4 toPseudoInertia_proxy(const Inertia & self)
      {
        return self.toPseudoInertia();
      }

      static Inertia *
      makeFromMCI(const Scalar & mass, const Vector3 & lever, const Matrix3 & inertia)
      {
#ifndef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS
        if (!inertia.isApprox(inertia.transpose()))
          throw eigenpy::Exception("The 3d inertia should be symmetric.");
        if (
          (Vector3::UnitX().transpose() * inertia * Vector3::UnitX() < 0)
          || (Vector3::UnitY().transpose() * inertia * Vector3::UnitY() < 0)
          || (Vector3::UnitZ().transpose() * inertia * Vector3::UnitZ() < 0))
          throw eigenpy::Exception("The 3d inertia should be positive.");
#endif
        return new Inertia(mass, lever, inertia);
      }

      static void expose()
      {
#if PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION == 6 && EIGENPY_VERSION_AT_LEAST(2, 9, 0)
        typedef PINOCCHIO_SHARED_PTR_HOLDER_TYPE(Inertia) HolderType;
#else
        typedef ::boost::python::detail::not_specified HolderType;
#endif
        bp::class_<Inertia, HolderType>(
          "Inertia",
          "This class represenses a sparse version of a Spatial Inertia and its is defined by its "
          "mass, its center of mass location and the rotational inertia expressed around this "
          "center of mass.\n\n"
          "Supported operations ...",
          bp::no_init)
          .def(InertiaPythonVisitor<Inertia>())
          .def(CastVisitor<Inertia>())
          .def(ExposeConstructorByCastVisitor<Inertia, ::pinocchio::Inertia>())
          .def(CopyableVisitor<Inertia>())
          .def(PrintableVisitor<Inertia>());
      }

    private:
      static Matrix6 __array__(const Inertia & self, bp::object, bp::object)
      {
        return self.matrix();
      }

      struct Pickle : bp::pickle_suite
      {
        static boost::python::tuple getinitargs(const Inertia & I)
        {
          return bp::make_tuple(I.mass(), (Vector3)I.lever(), I.inertia().matrix());
        }

        static bool getstate_manages_dict()
        {
          return true;
        }
      };

    }; // struct InertiaPythonVisitor

  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_spatial_inertia_hpp__
