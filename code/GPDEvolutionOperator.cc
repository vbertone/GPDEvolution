//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

int main()
{
  // x-space grid (pretty small to make the code converge faster)
  //const apfel::Grid g{{apfel::SubGrid{30, 1e-5, 3}, apfel::SubGrid{20, 1e-1, 3}, apfel::SubGrid{10, 9e-1, 3}}};
  const apfel::Grid g{{apfel::SubGrid{8, 1e-5, 3}, apfel::SubGrid{8, 1e-1, 3}}};

  // Initial scale
  const double mu0 = 1;

  // Vector of thresholds (no thresholds to simplify the picture)
  const std::vector<double> Thresholds = {0, 0, 0};

  // Perturbative order
  const int PerturbativeOrder = 0;

  // Skewness
  const double xi = 0.09;

  // Running coupling
  apfel::AlphaQCD a{0.118, 91.1876, Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // Initialize GPD evolution objects
  const auto GpdObjOp = InitializeGpdObjects(g, Thresholds, xi, true);

  // Construct the evolution operators
  const auto EvolvedOps = BuildDglap(GpdObjOp, mu0, PerturbativeOrder, as);

  // Tabulate evolution operators
  const apfel::TabulateObject<apfel::Set<apfel::Operator>> TabulatedOps{*EvolvedOps, 30, 1, 100, 3};

  // Final scale
  const double mu = 10;

  // Get evolution operators at the final scale
  apfel::Set<apfel::Operator> tops = TabulatedOps.Evaluate(mu);

  // Get evolution operators for qq and gq (actually these are
  // sigma-sigma and g-sigma)
  const apfel::matrix<double> Oqq = tops.GetObjects().at(8).GetOperator()[0];
  const apfel::matrix<double> Ogq = tops.GetObjects().at(1).GetOperator()[0];

  std::cout << std::scientific;
  std::cout.precision(3);
  for (int i = 0; i < Oqq.size(0); i++)
    for (int j = 0; j < Oqq.size(1); j++)
      std::cout << Oqq(i, j) << (std::string) "\t" + (j == Oqq.size(1) - 1 ? "\n" : "");

  std::cout << "\n";
  for (int i = 0; i < Oqq.size(0); i++)
    for (int j = 0; j < Oqq.size(1); j++)
      std::cout << Ogq(i, j) << (std::string) "\t" + (j == Oqq.size(1) - 1 ? "\n" : "");

  return 0;
}
