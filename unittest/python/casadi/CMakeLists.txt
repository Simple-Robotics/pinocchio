#
# Copyright (c) 2020-2021 INRIA
#

set(${PROJECT_NAME}_PYTHON_CASADI_TESTS bindings_main_algo)

foreach(test ${${PROJECT_NAME}_PYTHON_CASADI_TESTS})
  add_python_unit_test("test-py-casadi-${test}" "unittest/python/casadi/${test}.py"
                       "bindings/python")
endforeach(test ${${PROJECT_NAME}_PYTHON_CASADI_TESTS})
