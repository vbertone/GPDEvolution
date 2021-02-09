// APFEL++
#include <apfel/apfelxx.h>
#include <apfel/betaqcd.h>

// Boost
#include <boost/math/special_functions/gegenbauer.hpp>

double InputUp(double const& x, int const& n)
{
  const double alpha = 1.5;
  return ( 1 - pow(x, 2) ) * boost::math::gegenbauer(n, alpha, x);
}

int main()
{
  // Final scale
  const double mu = sqrt(10);

  // Retrieve evolution parameters
  const int    pto   = 0;
  const double Qref  = 1;
  const double asref = 0.513993;
  const double mu0   = 1;

  // Vector of thresholds
  const std::vector<double> Thresholds = {0, 0, 0};

  // Running coupling
  apfel::AlphaQCD a{asref, Qref, Thresholds, pto};
  const  apfel::TabulateObject<double> Alphas{a, 100, 0.9, 100, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100, 1e-7, 3}, apfel::SubGrid{200, 0.9e-1, 3}, apfel::SubGrid{100, 9e-1, 3}}};

  // ERBL evolution kernel
  const int n = 4;
  double Vn = - 0.5 + 1. / ( n + 1 ) / ( n + 2 );
  for (int k = 2; k <= n + 1; k++)
    Vn -= 2. / k;
  Vn *= 2 * apfel::CF;
  const double b0 = apfel::beta0qcd(3);
  const double GammaERBL = pow(as(mu) / as(mu0), - Vn / b0);

  // Input distributions (in the QCD evolution basis)
  const auto InPDFs = [&] (double const& x, double const&) -> std::map<int, double>
    {
      // Call all functions only once.
      const double upv  = x * InputUp(x, n);
      const double dnv  = 0;
      const double glu  = 0;
      const double dbar = 0;
      const double ubar = 0;
      const double sbar = 0;

      // Construct QCD evolution basis conbinations.
      double const Gluon   = glu;
      double const Singlet = dnv + 2 * dbar + upv + 2 * ubar + 2 * sbar;
      double const T3      = upv + 2 * ubar - dnv - 2 * dbar;
      double const T8      = upv + 2 * ubar + dnv + 2 * dbar - 4 * sbar;
      double const Valence = upv + dnv;
      double const V3      = upv - dnv;

      // Fill in map in the QCD evolution basis.
      std::map<int, double> QCDEvMap;
      QCDEvMap[0]  = Gluon;
      QCDEvMap[1]  = Singlet;
      QCDEvMap[2]  = Valence;
      QCDEvMap[3]  = T3;
      QCDEvMap[4]  = V3;
      QCDEvMap[5]  = T8;
      QCDEvMap[6]  = Valence;
      QCDEvMap[7]  = Singlet;
      QCDEvMap[8]  = Valence;
      QCDEvMap[9]  = Singlet;
      QCDEvMap[10] = Valence;
      QCDEvMap[11] = Singlet;
      QCDEvMap[12] = Valence;

      return QCDEvMap;
    };

  // DGLAP evolution
  apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*(BuildDglap(InitializeDglapObjectsQCD(g, Thresholds), InPDFs, mu0, pto, as)), 30, 1, 100, 3};

  // Construct the GPD evolution objects and tabulate GPDs
  const std::vector<double> xiv{1};
  std::vector<apfel::TabulateObject<apfel::Set<apfel::Distribution>>> TabulatedGPDs;
  for(double const& xi : xiv)
    TabulatedGPDs.push_back(apfel::TabulateObject<apfel::Set<apfel::Distribution>> {*(BuildDglap(InitializeGpdObjects(g, Thresholds, xi), InPDFs, mu0, pto, as)), 30, 1, 100, 3});

  // Print results
  const int nx = 100;
  const double xmin = 1e-3;
  const double xmax = 0.99;
  const double xstp = ( xmax - xmin ) / ( nx - 1 );
  for (double x = xmin; x <= xmax; x += xstp)
    {
      std::cout << std::scientific << x << "\t";
      const std::map<int, double> DistMapPDFs = apfel::QCDEvToPhys(TabulatedPDFs.EvaluateMapxQ(x, mu));
      std::cout << DistMapPDFs.at(2) << "\t";
      for (auto const& tgpd : TabulatedGPDs)
	{
	  const std::map<int, double> DistMapGPDs = apfel::QCDEvToPhys(tgpd.EvaluateMapxQ(x, mu));
	  std::cout << ( DistMapGPDs.at(2) - DistMapGPDs.at(-2) ) / ( x * InputUp(x, n) ) << "\t";
	}
      std::cout << GammaERBL << "\t";
      std::cout << std::endl;
    }
  std::cout << "========\n";

  return 0;
}
