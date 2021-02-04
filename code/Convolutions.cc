//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>
#include <apfel/splittingfunctionsunp_sl.h>
#include <apfel/gpdsplittingfunctionsunp_sl.h>

// Test distribution
double fy(double const& y)
{
  return y * ( 1 - y );
}

// Main program
int main()
{
  // Tabulation parameters
  const int nx = 100;
  const double xmin = 0.00001;
  const double xmax = 0.99;
  const double xstp = exp( log( xmax /  xmin ) / ( nx - 1 ) );
  const double xi = 0.3;

  // x-space grid
  const apfel::Grid g{{{100, 1e-7, 3}, {200, 1e-1, 3}, {100, 9e-1, 3}}};

  const apfel::Pgpd0ns pns{xi};
  const apfel::Distribution dfns = apfel::OperatorGPD{g, pns}  * ( [] (double const& y) -> double { return y; } * apfel::Distribution{g, fy} );

  const apfel::Pgpd0qq pqq{xi};
  const apfel::Distribution dfqq = apfel::OperatorGPD{g, pqq}  * ( [] (double const& y) -> double { return y; } * apfel::Distribution{g, fy} );

  const apfel::Pgpd0gq pgq{xi};
  const apfel::Distribution dfgq = apfel::OperatorGPD{g, pgq}  * ( [] (double const& y) -> double { return y; } * apfel::Distribution{g, fy} );

  const apfel::Pgpd0qg pqg{3, xi};
  const apfel::Distribution dfqg = apfel::OperatorGPD{g, pqg}  * ( [] (double const& y) -> double { return y; } * apfel::Distribution{g, fy} );

  const apfel::Pgpd0gg pgg{3, xi};
  const apfel::Distribution dfgg = apfel::OperatorGPD{g, pgg}  * ( [] (double const& y) -> double { return y; } * apfel::Distribution{g, fy} );

  // DGLAP kernels for comparison
  const apfel::P0ns pnsf{};
  const apfel::Distribution dfnsf = apfel::Operator{g, pnsf}  * ( [] (double const& y) -> double { return y; } * apfel::Distribution{g, fy} );

  const apfel::P0gq pgqf{};
  const apfel::Distribution dfgqf = apfel::Operator{g, pgqf}  * ( [] (double const& y) -> double { return y; } * apfel::Distribution{g, fy} );

  const apfel::P0qg pqgf{3};
  const apfel::Distribution dfqgf = apfel::Operator{g, pqgf}  * ( [] (double const& y) -> double { return y; } * apfel::Distribution{g, fy} );

  const apfel::P0gg pggf{3};
  const apfel::Distribution dfggf = apfel::Operator{g, pggf}  * ( [] (double const& y) -> double { return y; } * apfel::Distribution{g, fy} );

  for (double x = xmin; x <= xmax; x *= xstp)
    {
      // Explicit integration
      pns.SetExternalVariable(x);
      const apfel::Integrator Ins{[=] (double const& y) -> double { return x * ( pns.Regular(x/y) * fy(y)
										 + pns.Singular(x/y) * ( fy(y) - x / y * fy(x) * ( 1 + (x > y ? y / x - 1 : 0) ) ) ) / y; }};

      pqq.SetExternalVariable(x);
      const apfel::Integrator Iqq{[=] (double const& y) -> double { return x * ( pqq.Regular(x/y) * fy(y)
										 + pqq.Singular(x/y) * ( fy(y) - x / y * fy(x) * ( 1 + (x > y ? y / x - 1 : 0) ) ) ) / y; }};

      pgq.SetExternalVariable(x);
      const apfel::Integrator Igq{[=] (double const& y) -> double { return x * pgq.Regular(x/y) * fy(y) / y; }};

      pqg.SetExternalVariable(x);
      const apfel::Integrator Iqg{[=] (double const& y) -> double { return x * pqg.Regular(x/y) * fy(y) / y; }};

      pgg.SetExternalVariable(x);
      const apfel::Integrator Igg{[=] (double const& y) -> double { return x * ( pgg.Regular(x/y) * fy(y)
										 + pgg.Singular(x/y) * ( fy(y) - x / y * fy(x) * ( 1 + (x > y ? y / x - 1 : 0) ) ) ) / y; }};

      // Compute the convolutions
      const double interns  = dfns.Evaluate(x);
      const double internsf = dfnsf.Evaluate(x);
      const double integns  = Ins.integrate(0, 1, 1e-9) + x * fy(x) * ( pns.Local(x) + pns.LocalPV(x) );

      const double interqq  = dfqq.Evaluate(x);
      const double interqqf = dfnsf.Evaluate(x);
      const double integqq  = Iqq.integrate(0, 1, 1e-9) + x * fy(x) * ( pqq.Local(x) + pqq.LocalPV(x) );

      const double intergq  = dfgq.Evaluate(x);
      const double intergqf = dfgqf.Evaluate(x);
      const double integgq  = Igq.integrate(0, 1, 1e-9);

      const double interqg  = dfqg.Evaluate(x);
      const double interqgf = dfqgf.Evaluate(x);
      const double integqg  = Iqg.integrate(0, 1, 1e-9);

      const double intergg  = dfgg.Evaluate(x);
      const double interggf = dfggf.Evaluate(x);
      const double integgg  = Igg.integrate(0, 1, 1e-9) + x * fy(x) * ( pgg.Local(x) + pgg.LocalPV(x) );

      std::cout << std::scientific << x << "\t"
		<< interns << "\t" << integns << "\t" << interns / integns << "\t"
		<< interqq << "\t" << integqq << "\t" << interqq / integqq << "\t"
		//<< intergq << "\t" << integgq << "\t" << intergq / integgq << "\t"
		//<< interqg << "\t" << integqg << "\t" << interqg / integqg << "\t"
		<< intergg << "\t" << integgg << "\t" << intergg / integgg << "\t"
		<< std::endl;
    }
  return 0;
}
