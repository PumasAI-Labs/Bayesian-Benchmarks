#include <stan/math/torsten/test/unit/test_fixture_onecpt.hpp>
#include <stan/math/torsten/test/unit/test_functors.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <boost/mp11.hpp>

TYPED_TEST_SUITE_P(test_onecpt);

TYPED_TEST_P(test_onecpt, single_bolus_with_tlag) {
  this -> reset_events(2);      // only need two events
  this -> amt[0] = 1000;
  this -> evid[0] = 1;
  this -> cmt[0] = 1;
  this -> ii[0] = 0;
  this -> addl[0] = 0;
  this -> time[0] = 0.0;
  this -> tlag[0][0] = 1.5;
  this -> time[1] = 2.5;

  this -> compare_param_overload(1.e-10);
}

TYPED_TEST_P(test_onecpt, multiple_bolus) {
  this -> compare_overload(1.e-10);
}

TYPED_TEST_P(test_onecpt, multiple_bolus_cent) {
  this -> cmt[0] = 2;
  this -> compare_overload(1.e-10);
}

TYPED_TEST_P(test_onecpt, multiple_bolus_addl) {
  this -> ii[0] = 1.3;          // ensure test II + ADDL by end of time
  this -> compare_overload(1.e-10);
}

TYPED_TEST_P(test_onecpt, multiple_bolus_with_tlag) {
  this -> tlag[0][0] = 5.7;  // FIXME: tlag-generated event coincidents another event the test would fail
  this -> tlag[0][1] = 0;  // tlag2
  this -> time[0] = 0;
  for(int i = 1; i < 10; ++i) this -> time[i] = this -> time[i - 1] + 1;
  this -> amt[0] = 1200;

  this -> compare_param_overload(1.e-10);
}

TYPED_TEST_P(test_onecpt, multiple_bolus_ss) {
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < 10; i++) this -> time[i] = this -> time[i - 1] + 5;
  this -> amt[0] = 1200;
  this -> addl[0] = 10;
  this -> ss[0] = 1;

  this -> compare_overload(1.e-10);
}

TYPED_TEST_P(test_onecpt, multiple_infusion) {
  this -> cmt[0] = 2;
  this -> rate[0] = 350;
  this -> addl[0] = 2;

  this -> biovar[0] = {0.8, 0.9};
  this -> tlag[0] = {0.4, 0.8};

  this -> compare_param_overload(1.e-10);
}

TYPED_TEST_P(test_onecpt, multiple_infusion_ss) {
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < 10; i++) this -> time[i] = this -> time[i - 1] + 5;
  this -> amt[0] = 1200;
  this -> rate[0] = 190;        // amt/rate = even time will cause test failure
  this -> addl[0] = 10;
  this -> ss[0] = 1;

  this -> compare_overload(1.e-10);
}

REGISTER_TYPED_TEST_SUITE_P(test_onecpt,
                            single_bolus_with_tlag,
                            multiple_bolus,
                            multiple_bolus_cent,
                            multiple_bolus_addl,
                            multiple_bolus_with_tlag,
                            multiple_bolus_ss,
                            multiple_infusion,
                            multiple_infusion_ss);

using onecpt_test_types = boost::mp11::mp_product<
  std::tuple,
  ::testing::Types<pmx_solve_onecpt_functor>, // solver 1
  ::testing::Types<pmx_solve_onecpt_functor>, // solver 2
  ::testing::Types<stan::math::var_value<double>>,  // TIME
  ::testing::Types<double, stan::math::var_value<double>>,  // AMT
  ::testing::Types<double, stan::math::var_value<double>> , // RATE
  ::testing::Types<stan::math::var_value<double>> , // II
  ::testing::Types<double, stan::math::var_value<double>> , // PARAM
  ::testing::Types<double, stan::math::var_value<double>> , // BIOVAR
  ::testing::Types<double> , // TLAG
  ::testing::Types<torsten::PMXOneCptODE>                   // ODE
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(PMX, test_onecpt, onecpt_test_types);
