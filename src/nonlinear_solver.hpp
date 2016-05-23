#ifndef NOSH_NONLINEAR_SOLVER_HPP
#define NOSH_NONLINEAR_SOLVER_HPP

#include <LOCA_Thyra_SaveDataStrategy.H>
#include <Thyra_ModelEvaluatorDefaultBase.hpp>
#include <Tpetra_Vector.hpp>

#include <boost/any.hpp>

#include <map>
#include <memory>

using list = std::map<std::string, boost::any>;
namespace nosh {
  std::shared_ptr<const Tpetra::Vector<double,int,int>>
  nonlinear_solve(
      const std::shared_ptr<Thyra::ModelEvaluatorDefaultBase<double>> & model,
      std::map<std::string, boost::any> solver_params
      );

  void
  parameter_continuation(
      const std::shared_ptr<Thyra::ModelEvaluatorDefaultBase<double>> & model,
      const std::shared_ptr<LOCA::Thyra::SaveDataStrategy> & data_saver,
      std::map<std::string, boost::any> solver_params
      );
}
#endif // NOSH_NONLINEAR_SOLVER_HPP
