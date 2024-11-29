//
// Copyright (c) 2024 INRIA
//

#include "pinocchio/parsers/mjcf.hpp"
#include "pinocchio/bindings/python/parsers/mjcf.hpp"
#include "pinocchio/bindings/python/utils/path.hpp"

#include <boost/python.hpp>
#include <boost/filesystem/fstream.hpp>

namespace pinocchio
{
  namespace python
  {

    namespace bp = boost::python;

    Model buildModelFromMJCF(const bp::object & filename)
    {
      Model model;
      ::pinocchio::mjcf::buildModel(path(filename), model);
      return model;
    }

    Model buildModelFromMJCF(const bp::object & filename, Model & model)
    {
      return ::pinocchio::mjcf::buildModel(path(filename), model);
    }

    bp::tuple buildModelFromMJCF(
      const bp::object & filename,
      Model & model,
      PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(BilateralPointConstraintModel) & constraint_models)
    {
      ::pinocchio::mjcf::buildModel(path(filename), model, constraint_models);
      return bp::make_tuple(model, constraint_models);
    }

    Model buildModelFromMJCF(const bp::object & filename, const JointModel & root_joint)
    {
      Model model;
      ::pinocchio::mjcf::buildModel(path(filename), root_joint, model);
      return model;
    }

    Model
    buildModelFromMJCF(const bp::object & filename, const JointModel & root_joint, Model & model)
    {
      return ::pinocchio::mjcf::buildModel(path(filename), root_joint, model);
    }

    bp::tuple buildModelFromMJCF(
      const bp::object & filename,
      const JointModel & root_joint,
      Model & model,
      PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(BilateralPointConstraintModel) & constraint_models)
    {
      ::pinocchio::mjcf::buildModel(path(filename), root_joint, model, constraint_models);
      return bp::make_tuple(model, constraint_models);
    }

    bp::tuple buildModelFromMJCF(
      const bp::object & filename,
      const JointModel & root_joint,
      const std::string & root_joint_name)
    {
      Model model;
      PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(BilateralPointConstraintModel) constraint_models;
      ::pinocchio::mjcf::buildModel(
        path(filename), root_joint, root_joint_name, model, constraint_models);
      return bp::make_tuple(model, constraint_models);
    }

    bp::tuple buildModelFromMJCF(
      const bp::object & filename,
      const JointModel & root_joint,
      const std::string & root_joint_name,
      Model & model)
    {
      PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(BilateralPointConstraintModel) constraint_models;
      ::pinocchio::mjcf::buildModel(
        path(filename), root_joint, root_joint_name, model, constraint_models);
      return bp::make_tuple(model, constraint_models);
    }

    bp::tuple buildModelFromMJCF(
      const bp::object & filename,
      const JointModel & root_joint,
      const std::string & root_joint_name,
      Model & model,
      PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(BilateralPointConstraintModel) & constraint_models)
    {
      ::pinocchio::mjcf::buildModel(
        path(filename), root_joint, root_joint_name, model, constraint_models);
      return bp::make_tuple(model, constraint_models);
    }

    Model buildModelFromMJCFString(const std::string & xml_string)
    {
      Model model;
      ::pinocchio::mjcf::buildModelFromXML(xml_string, model);
      return model;
    }

    Model buildModelFromMJCFString(const std::string & xml_string, const JointModel & root_joint)
    {
      Model model;
      ::pinocchio::mjcf::buildModelFromXML(xml_string, root_joint, model);
      return model;
    }

    bp::tuple buildModelFromMJCFString(
      const std::string & xml_string,
      const JointModel & root_joint,
      const std::string & root_joint_name)
    {
      Model model;
      PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(BilateralPointConstraintModel) constraint_models;
      ::pinocchio::mjcf::buildModelFromXML(
        xml_string, root_joint, root_joint_name, model, constraint_models);
      return bp::make_tuple(model, constraint_models);
    }

    void exposeMJCFModel()
    {
      bp::def(
        "buildModelFromMJCF",
        static_cast<Model (*)(const bp::object &)>(pinocchio::python::buildModelFromMJCF),
        bp::args("mjcf_filename"),
        "Parse the MJCF file given in input and return a pinocchio Model.");

      bp::def(
        "buildModelFromMJCF",
        static_cast<Model (*)(const bp::object &, Model &)>(pinocchio::python::buildModelFromMJCF),
        bp::args("mjcf_filename", "model"),
        "Parse the MJCF file given in input and return a pinocchio Model.");

      bp::def(
        "buildModelFromMJCF",
        static_cast<bp::tuple (*)(
          const bp::object &, Model &,
          PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(BilateralPointConstraintModel) &)>(
          pinocchio::python::buildModelFromMJCF),
        bp::args("mjcf_filename", "model", "constraint_models"),
        "Parse the MJCF file given in input and return a pinocchio Model.");

      bp::def(
        "buildModelFromMJCF",
        static_cast<Model (*)(const bp::object &, const JointModel &)>(
          pinocchio::python::buildModelFromMJCF),
        bp::args("mjcf_filename", "root_joint"),
        "Parse the MJCF file and return a pinocchio Model with the given root Joint.");

      bp::def(
        "buildModelFromMJCF",
        static_cast<Model (*)(const bp::object &, const JointModel &, Model &)>(
          pinocchio::python::buildModelFromMJCF),
        bp::args("mjcf_filename", "root_joint", "model"),
        "Parse the MJCF file and return a pinocchio Model with the given root Joint.");

      bp::def(
        "buildModelFromMJCF",
        static_cast<bp::tuple (*)(
          const bp::object &, const JointModel &, Model &,
          PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(BilateralPointConstraintModel) &)>(
          pinocchio::python::buildModelFromMJCF),
        bp::args("mjcf_filename", "root_joint", "model", "constraint_models"),
        "Parse the MJCF file and return a pinocchio Model with the given root Joint and the "
        "constraint models.");

      bp::def(
        "buildModelFromMJCF",
        static_cast<bp::tuple (*)(const bp::object &, const JointModel &, const std::string &)>(
          pinocchio::python::buildModelFromMJCF),
        bp::args("mjcf_filename", "root_joint", "root_joint_name"),
        "Parse the MJCF file and return a pinocchio Model with the given root Joint and its "
        "specified name as well as a constraint list if some are present in the MJCF file.");

      bp::def(
        "buildModelFromMJCF",
        static_cast<bp::tuple (*)(
          const bp::object &, const JointModel &, const std::string &, Model &)>(
          pinocchio::python::buildModelFromMJCF),
        bp::args("mjcf_filename", "root_joint", "root_joint_name", "model"),
        "Parse the MJCF file and return a pinocchio Model with the given root Joint and its "
        "specified name as well as a constraint list if some are present in the MJCF file.");

      bp::def(
        "buildModelFromMJCF",
        static_cast<bp::tuple (*)(
          const bp::object &, const JointModel &, const std::string &, Model &,
          PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(BilateralPointConstraintModel) &)>(
          pinocchio::python::buildModelFromMJCF),
        bp::args("mjcf_filename", "root_joint", "root_joint_name", "model", "constraint_models"),
        "Parse the MJCF file and return a pinocchio Model with the given root Joint and its "
        "specified name as well as a constraint list if some are present in the MJCF file.");

      bp::def(
        "buildModelFromMJCFString",
        static_cast<Model (*)(const std::string &)>(pinocchio::python::buildModelFromMJCFString),
        bp::args("xml_string"),
        "Parse the MJCF string given in input and return a pinocchio Model.");

      bp::def(
        "buildModelFromMJCFString",
        static_cast<Model (*)(const std::string &, const JointModel &)>(
          pinocchio::python::buildModelFromMJCFString),
        bp::args("mjcf_filename", "root_joint"),
        "Parse the MJCF string and return a pinocchio Model with the given root Joint.");

      bp::def(
        "buildModelFromMJCFString",
        static_cast<bp::tuple (*)(const std::string &, const JointModel &, const std::string &)>(
          pinocchio::python::buildModelFromMJCFString),
        bp::args("mjcf_filename", "root_joint", "root_joint_name"),
        "Parse the MJCF string and return a pinocchio Model with the given root Joint and its "
        "specified name as well as a constraint list if some are present in the MJCF file.");
    }
  } // namespace python
} // namespace pinocchio
