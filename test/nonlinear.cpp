#include <catch.hpp>

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_ModelEvaluatorDefaultBase.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <mikado.hpp>

using dict = std::map<std::string, boost::any>;
// ===========================================================================
class jac_sqrt_two: public Tpetra::Operator<double,int,int>
{
  public:
    explicit
    jac_sqrt_two(
        const Tpetra::Vector<double,int,int> & x0
        ):
      x0_(x0)
    {
    }

    virtual
    ~jac_sqrt_two()
    {
    }

    virtual
    Teuchos::RCP<const Tpetra::Map<int,int>>
    getDomainMap () const
    {
      return x0_.getMap();
    }

    virtual
    Teuchos::RCP<const Tpetra::Map<int,int>>
    getRangeMap () const
    {
      return x0_.getMap();
    }

    void
    apply(
        const Tpetra::MultiVector<double,int,int> & X,
        Tpetra::MultiVector<double,int,int> & Y,
        Teuchos::ETransp mode=Teuchos::NO_TRANS,
        double alpha=Teuchos::ScalarTraits<double>::one(),
        double beta=Teuchos::ScalarTraits<double>::zero()
        ) const
    {
      for (size_t k = 0; k < Y.getNumVectors(); k++) {
        const auto x_data = X.getData(k);
        const auto x0_data = x0_.getData();
        auto y_data = Y.getDataNonConst(k);
        for (size_t i = 0; i < y_data.size(); i++) {
          y_data[i] = 2 * x0_data[i] * x_data[i];
        }
      }
      return;
    }

    void set_x0(const Tpetra::Vector<double,int,int> & x0)
    {
      x0_ = x0;
      return;
    }

  private:
    Tpetra::Vector<double,int,int> x0_;
};
// ===========================================================================
class sqrt_two : public Thyra::ModelEvaluatorDefaultBase<double>
{
public:
  explicit
  sqrt_two (
      const Teuchos::RCP<const Tpetra::Vector<double,int,int>> & init_x
    ):
    comm_(init_x->getMap()->getComm()),
    jac_(Teuchos::rcp(new jac_sqrt_two(*init_x))),
    space_(Thyra::createVectorSpace<double>(init_x->getMap())),
    p_map_(init_p_map_()),
    nominal_values_(init_nominal_values_(init_x))
  {
  }

  virtual
  ~sqrt_two() {};

  virtual
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
  get_x_space() const
  {
#ifndef NDEBUG
    TEUCHOS_ASSERT(!space_.is_null());
#endif
    return space_;
  }

  virtual
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
  get_f_space() const
  {
#ifndef NDEBUG
    TEUCHOS_ASSERT(!space_.is_null());
#endif
    return space_;
  }

  virtual
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
  get_p_space(int l) const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        l != 0,
        "Nosh can only deal with one parameter vector."
        );
    return Thyra::createVectorSpace<double>(p_map_);
  }

  virtual
  Teuchos::RCP<const Teuchos::Array<std::string> >
  get_p_names(int l) const
  {
    (void) l;
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
    return Teuchos::null;
  }

  virtual
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
  get_g_space(int l) const
  {
    (void) l;
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
    return Teuchos::null;
  }

  virtual
  Teuchos::ArrayView<const std::string>
  get_g_names(int j) const
  {
    (void) j;
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
    return Teuchos::null;
  }

  virtual
  Thyra::ModelEvaluatorBase::InArgs<double>
  getNominalValues() const
  {
    return nominal_values_;
  }

  virtual
  Thyra::ModelEvaluatorBase::InArgs<double>
  getLowerBounds() const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
    return this->createInArgs();
  }

  virtual
  Thyra::ModelEvaluatorBase::InArgs<double>
  getUpperBounds() const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
    return this->createInArgs();
  }

  virtual
  Teuchos::RCP<Thyra::LinearOpBase<double>>
  create_W_op() const
  {
    Teuchos::RCP<Tpetra::Operator<double,int,int>> op = jac_;
    return Thyra::createLinearOp(op, space_, space_);
  }

  virtual
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double>>
  get_W_factory() const
  {
    Stratimikos::DefaultLinearSolverBuilder builder;

    const std::map<std::string, boost::any> linear_solver_params = {
      {"package", std::string("Belos")},
      {"method", std::string("Pseudo Block GMRES")},
      {"parameters", dict{
        {"Output Frequency", 1},
        {"Output Style", 1},
        {"Verbosity", 33}
      }}
      };
    auto p = Teuchos::rcp(new Teuchos::ParameterList());
    mikado::std_map_to_teuchos_list(
        mikado::convert_to_belos_parameters(linear_solver_params),
        *p
        );
    builder.setParameterList(p);

    auto lowsFactory = builder.createLinearSolveStrategy("");

    lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

    return lowsFactory;
  }

  virtual
  Teuchos::RCP<Thyra::PreconditionerBase<double>>
  create_W_prec() const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
    return Teuchos::null;
  }

  virtual
  Thyra::ModelEvaluatorBase::InArgs<double>
  createInArgs() const
  {
    Thyra::ModelEvaluatorBase::InArgsSetup<double> in_args;

    in_args.setModelEvalDescription("sqrt(2)");

    in_args.set_Np(1);
    in_args.setSupports(IN_ARG_x, true);

    return in_args;
  }

  virtual
  void
  reportFinalPoint(
      const Thyra::ModelEvaluatorBase::InArgs<double> &finalPoint,
      const bool wasSolved
      )
  {
    (void) finalPoint;
    (void) wasSolved;
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
  }


protected:

  virtual
  Thyra::ModelEvaluatorBase::OutArgs<double>
  createOutArgsImpl() const
  {
    Thyra::ModelEvaluatorBase::OutArgsSetup<double> out_args;
    out_args.setModelEvalDescription("sqrt(2)");
    out_args.set_Np_Ng(1, 0);
    out_args.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f);
    out_args.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op);
    return out_args;
  }

  virtual
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<double> &in_args,
      const Thyra::ModelEvaluatorBase::OutArgs<double> &out_args
      ) const
  {
    const auto & x_in = in_args.get_x();
#ifndef NDEBUG
    TEUCHOS_ASSERT(!x_in.is_null());
#endif
    // create corresponding tpetra vector
    auto x_in_tpetra =
      Thyra::TpetraOperatorVectorExtraction<double,int,int>::getConstTpetraVector(
          x_in
          );

    // compute F
    const auto & f_out = out_args.get_f();
    if (!f_out.is_null()) {
      auto f_out_tpetra =
        Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraVector(
            f_out
            );
      const auto x_data = x_in_tpetra->getData();
      auto f_data = f_out_tpetra->getDataNonConst();
      for (size_t i = 0; i < f_data.size(); i++) {
        f_data[i] = x_data[i] * x_data[i] - 2.0;
      }
    }

    // Fill Jacobian.
    const auto & W_out = out_args.get_W_op();
    if(!W_out.is_null()) {
      auto W_outT =
        Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraOperator(
            W_out
            );
      const auto & jac = Teuchos::rcp_dynamic_cast<jac_sqrt_two>(W_outT, true);
      jac->set_x0(*x_in_tpetra);
    }

    return;
  }

private:
  Teuchos::RCP<Tpetra::Map<int,int>>
  init_p_map_() const
  {
    return Teuchos::rcp(new Tpetra::Map<int,int>(
          1,  // number of params
          0,
          comm_
          ));
  }

  Thyra::ModelEvaluatorBase::InArgs<double>
  init_nominal_values_(
      const Teuchos::RCP<const Tpetra::Vector<double,int,int>> & x
      )
  {
    auto nominal_values = this->createInArgs();
    const auto xxx = Thyra::createConstVector(x, space_);
    nominal_values.set_x(xxx);
    return nominal_values;
  }

private:
  const Teuchos::RCP<const Tpetra::Comm<int>> comm_;
  const Teuchos::RCP<jac_sqrt_two> jac_;
  const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> space_;
  Teuchos::RCP<const Tpetra::Map<int,int>> p_map_;
  Thyra::ModelEvaluatorBase::InArgs<double> nominal_values_;
};
// ===========================================================================
TEST_CASE("nonlinear solver", "[nonlinear]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const int n = 100;
  const auto map = Teuchos::rcp(new Tpetra::Map<int,int>(n, 0, comm));

  auto x = Teuchos::rcp(new Tpetra::Vector<double,int,int>(map));
  x->putScalar(2.0);

  const auto model = std::make_shared<sqrt_two>(x);

  const auto sol = mikado::nonlinear_solve(
      model,
      {
        {"method", std::string("Newton")}
      }
      );

  const auto sol_data = sol->getData();
  REQUIRE(sol_data[0] == Approx(std::sqrt(2.0)));
  REQUIRE(sol_data[1] == Approx(std::sqrt(2.0)));
  REQUIRE(sol_data[2] == Approx(std::sqrt(2.0)));
}
// ===========================================================================
