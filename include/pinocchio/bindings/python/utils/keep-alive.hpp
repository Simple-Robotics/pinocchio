//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_python_utils_keep_alive_hpp__
#define __pinocchio_python_utils_keep_alive_hpp__

#include <eigenpy/eigenpy.hpp>

namespace pinocchio
{
  namespace python
  {
    template<
      typename ReturnType,
      std::size_t custodian,
      std::size_t ward,
      class BasePolicy_ = boost::python::default_call_policies>
    struct keep_alive;

    template<std::size_t custodian, std::size_t ward, class BasePolicy_>
    struct keep_alive<boost::python::tuple, custodian, ward, BasePolicy_> : BasePolicy_
    {
      BOOST_STATIC_ASSERT(custodian != ward);

      template<class ArgumentPackage>
      static PyObject * postcall(ArgumentPackage const & args_, PyObject * result)
      {
        std::size_t arity_ = boost::python::detail::arity(args_);
        // check if either custodian or ward exceeds the arity
        // (this weird formulation avoids "always false" warnings
        // for arity_ = 0)
        if ((std::max)(custodian, ward) > arity_)
        {
          PyErr_SetString(PyExc_IndexError, "keep_alive_with_tuple: argument index out of range");
          return 0;
        }

        if (!PyTuple_Check(result))
        {
          PyErr_SetString(PyExc_RuntimeError, "keep_alive_with_tuple: result is not a tuple");
          return 0;
        }

        int tuple_size = PyTuple_Size(result);
        PINOCCHIO_THROW_PRETTY_IF(
          custodian < std::size_t(tuple_size), std::runtime_error,
          "custodian is greater than tuple_size.");

        // keep_alive<Nurse, Patient> indicates that the argument with index Patient should be kept
        // alive at least until the argument with index Nurse is freed by the garbage collector.
        PyObject * patient = boost::python::detail::get_prev<ward>::execute(args_, result);
        PyObject * nurse = PyTuple_GetItem(result, custodian);

        if (nurse == 0)
          return 0;

        nurse = BasePolicy_::postcall(args_, nurse);
        if (nurse == 0)
          return 0;

        if (boost::python::objects::make_nurse_and_patient(nurse, patient) == 0)
        {
          Py_XDECREF(nurse);
          return 0;
        }
        return result;
      }
    };

  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_utils_keep_alive_hpp__
