#include <catch.hpp>

#include <LOCA_Thyra_SaveDataStrategy.H>
#include <NOX_Thyra_Vector.H>
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_ModelEvaluatorDefaultBase.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <mikado.hpp>

using dict = std::map<std::string, boost::any>;
// ===========================================================================
class jac_sqrt_alpha: public Tpetra::Operator<double,int,int>
{
  public:
    explicit
    jac_sqrt_alpha(
        const Tpetra::Vector<double,int,int> & x0
        ):
      x0_(x0)
    {
    }

    virtual
    ~jac_sqrt_alpha()
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
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        double alpha = Teuchos::ScalarTraits<double>::one(),
        double beta = Teuchos::ScalarTraits<double>::zero()
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
class sqrt_alpha : public Thyra::ModelEvaluatorDefaultBase<double>
{
public:
  explicit
  sqrt_alpha (
      const Teuchos::RCP<const Tpetra::Vector<double,int,int>> & init_x,
      const double alpha0
    ):
    comm_(init_x->getMap()->getComm()),
    jac_(Teuchos::rcp(new jac_sqrt_alpha(*init_x))),
    space_(Thyra::createVectorSpace<double>(init_x->getMap())),
    p_map_(init_p_map_()),
    p_names_(init_p_names_()),
    nominal_values_(init_nominal_values_(init_x, alpha0))
  {
  }

  virtual
  ~sqrt_alpha() {};

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
        "Mikado can only deal with one parameter vector."
        );
    return Thyra::createVectorSpace<double>(p_map_);
  }

  virtual
  Teuchos::RCP<const Teuchos::Array<std::string> >
  get_p_names(int l) const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        l != 0,
        "Mikado can only deal with one parameter vector."
        );
    return p_names_;
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
      {"method", std::string("Pseudo Block CG")},
      {"parameters", dict{
        {"Output Frequency", 1},
        {"Verbosity", 0}
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

    in_args.setModelEvalDescription("sqrt(alpha)");

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
    out_args.setModelEvalDescription("sqrt(alpha)");
    out_args.set_Np_Ng(1, 0);
    out_args.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f);
    out_args.setSupports(
        Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,
        0,
        DerivativeSupport(DERIV_MV_BY_COL)
        );
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
      // Dissect in_args.get_p(0) into parameter sublists.
      const auto & p_in = in_args.get_p(0);
#ifndef NDEBUG
      TEUCHOS_ASSERT(!p_in.is_null());
      // Make sure the parameters aren't NaNs.
      for (int k = 0; k < p_in->space()->dim(); k++) {
        TEUCHOS_ASSERT(!std::isnan(Thyra::get_ele(*p_in, k)));
      }
#endif
      // Fill the parameters into a std::map.
      const auto param_names = this->get_p_names(0);
      const double alph = Thyra::get_ele(*p_in, 0);

      auto f_out_tpetra =
        Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraVector(
            f_out
            );
      const auto x_data = x_in_tpetra->getData();
      auto f_data = f_out_tpetra->getDataNonConst();
      for (size_t i = 0; i < f_data.size(); i++) {
        f_data[i] = x_data[i] * x_data[i] - alph;
      }
    }

    // Compute df/dp.
    const auto & derivMv = out_args.get_DfDp(0).getDerivativeMultiVector();
    const auto & dfdp_out = derivMv.getMultiVector();
    if (!dfdp_out.is_null()) {
      auto dfdp_out_tpetra =
        Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraMultiVector(
            dfdp_out
            );

      TEUCHOS_ASSERT_EQUALITY(dfdp_out_tpetra->getNumVectors(), 1);
      auto out = dfdp_out_tpetra->getVectorNonConst(0);
      auto out_data = out->getDataNonConst();
      for (size_t k = 0; k < out_data.size(); k++) {
        out_data[k] = -1.0;
      }
    }

    // Fill Jacobian.
    const auto & W_out = out_args.get_W_op();
    if(!W_out.is_null()) {
      auto W_outT =
        Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraOperator(
            W_out
            );
      const auto & jac = Teuchos::rcp_dynamic_cast<jac_sqrt_alpha>(W_outT, true);
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

  Teuchos::RCP<Teuchos::Array<std::string>>
  init_p_names_()
  {
    auto p_names = Teuchos::rcp(new Teuchos::Array<std::string>(1));
    (*p_names_)[0] = "alpha";
    return p_names;
  }

  Thyra::ModelEvaluatorBase::InArgs<double>
  init_nominal_values_(
      const Teuchos::RCP<const Tpetra::Vector<double,int,int>> & x,
      const double alpha0
      )
  {
    auto p_init = Thyra::createMember(this->get_p_space(0));
    Thyra::set_ele(0, alpha0, p_init());

    const auto xxx = Thyra::createConstVector(x, space_);

    auto nominal_values = this->createInArgs();
    nominal_values.set_p(0, p_init);
    nominal_values.set_x(xxx);

    return nominal_values;
  }

private:
  const Teuchos::RCP<const Tpetra::Comm<int>> comm_;
  const Teuchos::RCP<jac_sqrt_alpha> jac_;
  const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> space_;
  Teuchos::RCP<const Tpetra::Map<int,int>> p_map_;
  Teuchos::RCP<Teuchos::Array<std::string> > p_names_;
  Thyra::ModelEvaluatorBase::InArgs<double> nominal_values_;
};
// ===========================================================================
class continuation_data_saver: public LOCA::Thyra::SaveDataStrategy
{
  public:
  continuation_data_saver():
    latest_x(Teuchos::null)
  {};

  virtual ~continuation_data_saver() {};

  virtual
  void
  saveSolution(
      const NOX::Abstract::Vector &x,
      double p
      )
  {
    (void) p;

    // extract Tpetra vector
    const auto x_nox_thyra = dynamic_cast<const NOX::Thyra::Vector*>(&x);
    TEUCHOS_ASSERT(x_nox_thyra != nullptr);
    const auto x_thyra = x_nox_thyra->getThyraRCPVector();
    auto x_tpetra =
      Thyra::TpetraOperatorVectorExtraction<double,int,int>::getConstTpetraVector(
          x_thyra
          );

    // copy over x
    latest_x = x_tpetra;
  }

  public:
  Teuchos::RCP<const Tpetra::Vector<double,int,int>> latest_x;
};
// ===========================================================================
TEST_CASE("nonlinear solver", "[nonlinear]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const int n = 100;
  const auto map = Teuchos::rcp(new Tpetra::Map<int,int>(n, 0, comm));

  auto x = Teuchos::rcp(new Tpetra::Vector<double,int,int>(map));
  x->putScalar(1.0);

  const auto model = std::make_shared<sqrt_alpha>(x, 2.0);

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
TEST_CASE("parameter continuation", "[param cont]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const int n = 100;
  const auto map = Teuchos::rcp(new Tpetra::Map<int,int>(n, 0, comm));

  auto x = Teuchos::rcp(new Tpetra::Vector<double,int,int>(map));
  x->putScalar(1.0);

  const auto model = std::make_shared<sqrt_alpha>(x, 1.0);

  const auto saver = std::make_shared<continuation_data_saver>();

  mikado::parameter_continuation(
      model, saver,
      {
        {"NOX", dict{
          {"Status Tests", dict{
            {"Test Type", std::string("NormF")},
            {"Norm Type", std::string("Two Norm")},
            {"Tolerance", 1.0e-8}
          }},
          {"Printing", dict{
           {"Output Information", dict{
             {"Details", true},
             {"Outer Iteration", true},
             {"Outer Iteration Status Test", true},
             {"Inner Iteration", true},
             {"Linear Solver Details", true},
             {"Parameters", true},
             {"Warning", true},
             {"Debug", true},
             {"Test Details", true},
             {"Error", true},
             {"Stepper Iteration", true},
             {"Stepper Details", true},
             {"Stepper Parameters", true}
           }}
          }}
        }},
        {"LOCA", dict{
          {"Predictor", dict{
            {"Method", std::string("Tangent")}
          }},
          {"Stepper", dict{
            {"Continuation Method", std::string("Arc Length")},
            {"Continuation Parameter", std::string("alpha")},
            {"Initial Value", 1.0},
            {"Min Value", 0.0},
            {"Max Value", 3.0},
            {"Max Nonlinear Iterations", 5},
          }},
          {"Step Size", dict{
            {"Initial Step Size", 1.0e-1},
            {"Min Step Size", 1.0e-5},
            {"Max Step Size", 5.0e-1},
            {"Aggressiveness", 0.1}
          }}
        }}
      }
      );

  const auto sol_data = saver->latest_x->getData();
  REQUIRE(sol_data[0] == Approx(std::sqrt(3.0)));
  // REQUIRE(sol_data[1] == Approx(std::sqrt(3.0)));
  // REQUIRE(sol_data[2] == Approx(std::sqrt(3.0)));
}
// ===========================================================================
