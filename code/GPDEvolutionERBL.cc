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

double InputGluon(double const& x, int const& n)
{
  const double alpha = 2.5;
  return pow(1 - pow(x, 2), 2) * boost::math::gegenbauer(n - 1, alpha, x);
}

int main()
{
  apfel::SetVerbosityLevel(0);

  // Degree of the gegenbauer polynomials
  const int n = 4;

  // Retrieve evolution parameters
  const int    pto   = 0;
  const double Qref  = 1;
  const double asref = 0.513993;
  const double mu0   = 1;

  // Vector of thresholds
  const std::vector<double> Thresholds = {0, 0, 0};

  // Number of active flavours
  const int nf = Thresholds.size();

  // Running coupling
  apfel::AlphaQCD a{asref, Qref, Thresholds, pto};
  const  apfel::TabulateObject<double> Alphas{a, 100, 0.9, 100, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100, 1e-7, 3}, apfel::SubGrid{200, 0.9e-1, 3}, apfel::SubGrid{100, 9e-1, 3}}};

  // Input distributions (in the QCD evolution basis)
  const auto InPDFs = [&] (double const& x, double const&) -> std::map<int, double>
    {
      // Call all functions only once.
      const double upv  = x * InputUp(x, n);
      const double dnv  = 0;
      const double glu  = x * InputGluon(x, n);
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
  apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedGPDs{*(BuildDglap(InitializeGpdObjects(g, Thresholds, 1), InPDFs, mu0, pto, as)), 30, 1, 100, 3};

  // Print results
  const int nx = 10000;
  const double xmin = 1e-3;
  const double xmax = 0.99;
  const double xstp = ( xmax - xmin ) / ( nx - 1 );

  // Final scales
  const std::vector<double> muv = {1.000000001, 2, 5, 10, 99.99};

  // ERBL evolution kernels
  std::vector<double> GammaNS;
  std::vector<std::vector<double>> GammaSG;
  for (double mu : muv)
    {
      const double b0 = apfel::beta0qcd(3);
      double Sn = 0;
      for (int k = 1; k <= n + 1; k++)
	Sn += 1. / k;
      const double Vns = 2 * apfel::CF * ( 1.5 + 1. / ( n + 1 ) / ( n + 2 ) - 2 * Sn );
      const double Vqq = Vns;
      const double Vqg = 4 * nf * apfel::TR * ( 4. + 3. * n + n * n ) / ( 1. + n ) / ( 2. + n ) / ( 3. + n );
      const double Vgq = 2 * apfel::CF * ( 4. + 3. * n + n * n ) / n / ( 1. + n ) / ( 2. + n );
      const double Vgg = b0 + 4 * apfel::CA * ( - 1. / ( 1. + n ) / ( 2. + n ) + 3. / n / ( 3. + n ) - Sn );
      const double La  = - log(as(mu) / as(mu0)) / b0;

      // Non singlet
      GammaNS.push_back(exp(Vns * La));

      // Singlet
      const double Delta = sqrt( pow(Vqq - Vgg, 2) + 4 * Vqg * Vgq ) * La;
      const double fact  = exp( ( Vqq + Vgg ) * La / 2 );
      GammaSG.push_back({
	  fact * ( Delta * cosh( Delta / 2 ) + ( Vqq - Vgg ) * La * sinh( Delta / 2 ) ) / Delta,
	  2 * Vqg * La * fact * sinh( Delta / 2 ) / Delta,
	  2 * Vgq * La * fact * sinh( Delta / 2 ) / Delta,
	  fact * ( Delta * cosh( Delta / 2 ) - ( Vqq - Vgg ) * La * sinh( Delta / 2 ) ) / Delta});
    }

  for (double x = xmin; x <= xmax; x += xstp)
    {
      std::cout << std::scientific << x << "\t";
      for (int imu = 0; imu < (int) muv.size(); imu++)
	{
	  const std::map<int, double> DistMapGPDs = apfel::QCDEvToPhys(TabulatedGPDs.EvaluateMapxQ(x, muv[imu]));
	  std::cout << DistMapGPDs.at(2) - DistMapGPDs.at(-2) << "\t" << GammaNS[imu] * x * InputUp(x, n) << "\t";
	  //std::cout << DistMapGPDs.at(0) << "\t" << GammaSG[imu][2] * x * InputUp(x, n) + GammaSG[imu][3] * x * InputGluon(x, n) << "\t";
	}
      std::cout << std::endl;
    }

  return 0;
}
