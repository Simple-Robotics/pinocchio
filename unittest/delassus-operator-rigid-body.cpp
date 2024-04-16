//
// Copyright (c) 2024 INRIA
//

#include <iostream>
#include <memory>

#include "pinocchio/algorithm/cholesky.hpp"
#include "pinocchio/algorithm/contact-info.hpp"
#include "pinocchio/algorithm/contact-dynamics.hpp"
#include "pinocchio/algorithm/contact-cholesky.hxx"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/context.hpp"
#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/algorithm/delassus-operator-rigid-body.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;
using namespace Eigen;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(default_constructor_shared_ptr)
{
  typedef DelassusOperatorRigidBodyTpl<double,0,JointCollectionDefaultTpl,std::shared_ptr> DelassusOperatorRigidBodySharedPtr;
  typedef typename DelassusOperatorRigidBodySharedPtr::ConstraintModelVector ConstraintModelVector;
  typedef typename DelassusOperatorRigidBodySharedPtr::ConstraintDataVector ConstraintDataVector;

  std::shared_ptr<Model> model_shared_ptr = std::make_shared<Model>();
  Model * model_ptr = model_shared_ptr.get();
  Model & model = *model_ptr;
  buildModels::humanoidRandom(model,true);
  
  std::shared_ptr<Data> data_shared_ptr = std::make_shared<Data>(model);
  Data * data_ptr = data_shared_ptr.get();
//  Data & data = *data_ptr;
  
  std::shared_ptr<ConstraintModelVector> constraint_models_shared_ptr = std::make_shared<ConstraintModelVector>();
  ConstraintModelVector * constraint_models_ptr = constraint_models_shared_ptr.get();
  std::shared_ptr<ConstraintDataVector> constraint_datas_shared_ptr = std::make_shared<ConstraintDataVector>();
  ConstraintDataVector * constraint_datas_ptr = constraint_datas_shared_ptr.get();

  DelassusOperatorRigidBodySharedPtr delassus_operator(model_shared_ptr,
                                                              data_shared_ptr,
                                                              constraint_models_shared_ptr,
                                                              constraint_datas_shared_ptr);
  
  BOOST_CHECK(delassus_operator.size() == 0);
  BOOST_CHECK(&delassus_operator.model() == model_ptr);
  BOOST_CHECK(&delassus_operator.data() == data_ptr);
  BOOST_CHECK(&delassus_operator.constraint_models() == constraint_models_ptr);
  BOOST_CHECK(&delassus_operator.constraint_datas() == constraint_datas_ptr);
}

BOOST_AUTO_TEST_CASE(default_constructor_reference_wrapper)
{
  typedef DelassusOperatorRigidBodyTpl<double,0,JointCollectionDefaultTpl,std::reference_wrapper> DelassusOperatorRigidBodyReferenceWrapper;
  typedef typename DelassusOperatorRigidBodyReferenceWrapper::ConstraintModelVector ConstraintModelVector;
  typedef typename DelassusOperatorRigidBodyReferenceWrapper::ConstraintDataVector ConstraintDataVector;
  
  Model model;
  std::reference_wrapper<Model> model_ref = model;
  buildModels::humanoidRandom(model,true);
  
  Data data(model);
  std::reference_wrapper<Data> data_ref = data;
  
  ConstraintModelVector constraint_models;
  std::reference_wrapper<ConstraintModelVector> constraint_models_ref = constraint_models;
  ConstraintDataVector constraint_datas;
  std::reference_wrapper<ConstraintDataVector> constraint_datas_ref = constraint_datas;

  DelassusOperatorRigidBodyReferenceWrapper delassus_operator(model_ref, data_ref,
                                                                     constraint_models_ref,
                                                                     constraint_datas_ref);
  
  BOOST_CHECK(delassus_operator.size() == 0);
  BOOST_CHECK(&delassus_operator.model() == &model);
  BOOST_CHECK(&delassus_operator.data() == &data);
  BOOST_CHECK(&delassus_operator.constraint_models() == &constraint_models);
  BOOST_CHECK(&delassus_operator.constraint_datas() == &constraint_datas);
}

BOOST_AUTO_TEST_SUITE_END()
