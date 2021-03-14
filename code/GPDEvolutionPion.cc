// APFEL++
#include <apfel/apfelxx.h>

double PionUp(double const& m_x, double const& m_xi);

int main()
{
  // Final scale
  const double mu = sqrt(100);

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
  const apfel::Grid g{{apfel::SubGrid{100, 1e-9, 3}, apfel::SubGrid{200, 0.9e-1, 3}, apfel::SubGrid{100, 9e-1, 3}}};

  // Construct the GPD evolution objects and tabulate GPDs
  //const std::vector<double> xiv{0.1};
  const std::vector<double> xiv{0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999};
  std::vector<apfel::TabulateObject<apfel::Set<apfel::Distribution>>> TabulatedGPDs;
  for(double const& xi : xiv)
    {
      // Input distributions (in the QCD evolution basis)
      const auto InPDFs = [&] (double const& x, double const&) -> std::map<int, double>
	{
	  // Call all functions only once.
	  const double upv  = x * (PionUp(x, xi) + PionUp(-x, xi));
	  const double dnv  = 0;
	  const double glu  = 0;
	  const double dbar = 0;
	  const double ubar = - x * PionUp(-x, xi);
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
      TabulatedGPDs.push_back(apfel::TabulateObject<apfel::Set<apfel::Distribution>> {*(BuildDglap(InitializeGpdObjects(g, Thresholds, xi), InPDFs, mu0, pto, as)), 30, 1, 100, 3});
    }

  // Print results
  const int nx = 100;
  const double xmin = 1e-3;
  const double xmax = 0.99;
  const double xstp = ( xmax - xmin ) / ( nx - 1 );
  for (double x = xmin; x <= xmax; x += xstp)
    {
      std::cout << std::scientific << x << "\t";
      for (auto const& tgpd : TabulatedGPDs)
	{
	  const std::map<int, double> DistMapGPDs = apfel::QCDEvToPhys(tgpd.EvaluateMapxQ(x, mu));
	  std::cout << DistMapGPDs.at(0) << "\t";
	}
      std::cout << std::endl;
    }
  std::cout << "========\n";

  // Check polynomiality
  std::vector<std::map<int, apfel::Distribution>> DistGPDs;
  for (auto const& tgpd : TabulatedGPDs)
    DistGPDs.push_back(apfel::QCDEvToPhys(tgpd.Evaluate(mu).GetObjects()));

  for (int i = 0; i < (int) xiv.size(); i++)
    std::cout << "Up-valence 1st moment (xi = " << xiv[i] << ")   : "
	      << ( [] (double const& x) -> double { return 1 / x; } * ( DistGPDs[i].at(2) - DistGPDs[i].at(-2) ) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";
  for (int i = 0; i < (int) xiv.size(); i++)
    std::cout << "Up-valence 2st moment (xi = " << xiv[i] << ")   : "
	      << ( DistGPDs[i].at(3) + DistGPDs[i].at(-3) +
		   DistGPDs[i].at(2) + DistGPDs[i].at(-2) +
		   DistGPDs[i].at(1) + DistGPDs[i].at(-1) +
		   DistGPDs[i].at(0) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";
  for (int i = 0; i < (int) xiv.size(); i++)
    std::cout << "Up-valence 3rd moment (xi = " << xiv[i] << ")   : "
	      << ( [] (double const& x) -> double { return x; } * ( DistGPDs[i].at(2) - DistGPDs[i].at(-2) ) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";


  return 0;
}

double PionUp(double const& m_x, double const& m_xi)
{
  double trans = 0.;
  if (fabs(m_xi) >= 1. || fabs(m_x) >= 1.)
    trans = 0.;
  else if(m_x >= fabs(m_xi))
    trans = 30. * pow(1 - m_x, 2) * (pow(m_x, 2) - pow(m_xi, 2)) / pow(1 - pow(m_xi, 2), 2);
  else if(m_xi >= fabs(m_x))
      trans = 15. * (m_x - 1) * (pow(m_x, 2) - pow(m_xi, 2)) * (m_x + 2 * m_x * m_xi + pow(m_xi, 2)) / pow(m_xi, 3) / pow(1 + m_xi, 2) / 2.;
  else if(m_xi <= -fabs(m_x))
      trans = 15 * (1 - m_x) * (pow(m_x, 2) - pow(m_xi, 2)) * (m_x - 2 * m_x * m_xi + pow(m_xi, 2)) / pow(m_xi, 3) / pow(1 - m_xi, 2) / 2;

  return trans;
}
