// APFEL++
#include <apfel/apfelxx.h>

// LHAPDF libs
#include "LHAPDF/LHAPDF.h"

int main()
{
  // Open LHAPDF set
  LHAPDF::PDF* dist = LHAPDF::mkPDF("MMHT2014lo68cl");

  // Final scale
  const double mu = 10;

  // Retrieve evolution parameters
  const int    pto   = dist->orderQCD();
  const double Qref  = 91.1876;
  const double asref = dist->alphasQ(Qref);
  const double mc    = dist->quarkThreshold(4);
  const double mb    = dist->quarkThreshold(5);
  const double mt    = dist->quarkThreshold(6);
  const double mu0   = dist->qMin();

  // Vector of thresholds
  const std::vector<double> Thresholds = {0, 0, 0, mc, mb, mt};

  // Running coupling
  apfel::AlphaQCD a{asref, Qref, Thresholds, pto};
  const  apfel::TabulateObject<double> Alphas{a, 100, 0.9, 100, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // Input distributions (in the QCD evolution basis)
  const auto InPDFs = [&] (double const& x, double const& Q) -> std::map<int, double> { return apfel::PhysToQCDEv(dist->xfxQ(x,Q)); };

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100, 1e-7, 3}, apfel::SubGrid{200, 1e-1, 3}, apfel::SubGrid{100, 9e-1, 3}}};

  // Construct the PDF evolution objects
  const auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCDpol(g, Thresholds), InPDFs, mu0, pto, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs,  30, 1, 100, 3};

  // Construct the GPD evolution objects and tabulate GPDs
  const std::vector<double> xiv{0, 0.02, 0.05, 0.11, 0.2, 0.35, 0.5, 0.7, 0.9, 1};
  std::vector<apfel::TabulateObject<apfel::Set<apfel::Distribution>>> TabulatedGPDs;
  for(double const& xi : xiv)
    TabulatedGPDs.push_back(apfel::TabulateObject<apfel::Set<apfel::Distribution>> {*(BuildDglap(InitializeGpdObjectsPol(g, Thresholds, xi), InPDFs, mu0, pto, as)), 30, 1, 100, 3});

  // Print results
  const int nx = 1000;
  const double xmin = 1e-3;
  const double xmax = 0.99;
  const double xstp = exp( log(xmax / xmin) / ( nx - 1 ) ); 
  for (double x = xmin; x <= xmax; x *= xstp)
    {
      const std::map<int, double> DistMapPDFs = apfel::QCDEvToPhys(TabulatedPDFs.EvaluateMapxQ(x, mu));
      //std::cout << std::scientific << x << "\t" << DistMapPDFs.at(2) - DistMapPDFs.at(-2) << "\t";
      std::cout << std::scientific << x << "\t" << DistMapPDFs.at(2) + DistMapPDFs.at(-2) << "\t";
      //std::cout << std::scientific << x << "\t" << DistMapPDFs.at(0) << "\t";
      //std::cout << std::scientific << dist->xfxQ(2, x, mu) + dist->xfxQ(-2, x, mu) << "\t";
      for (auto const& tgpd : TabulatedGPDs)
	{
	  const std::map<int, double> DistMapGPDs = apfel::QCDEvToPhys(tgpd.EvaluateMapxQ(x, mu));
	  //std::cout << DistMapGPDs.at(2) - DistMapGPDs.at(-2) << "\t";
	  std::cout << DistMapGPDs.at(2) + DistMapGPDs.at(-2) << "\t";
	  //std::cout << DistMapGPDs.at(0) << "\t";
	}
      std::cout << std::endl;
    }
  std::cout << "\n";

  // Check polynomiality
  const std::map<int, apfel::Distribution> DistPDFs = apfel::QCDEvToPhys(TabulatedPDFs.Evaluate(mu).GetObjects());
 
  std::vector<std::map<int, apfel::Distribution>> DistGPDs;
  for (auto const& tgpd : TabulatedGPDs)
    DistGPDs.push_back(apfel::QCDEvToPhys(tgpd.Evaluate(mu).GetObjects()));

  const apfel::Distribution msr = TabulatedPDFs.Evaluate(mu).at(0) + TabulatedPDFs.Evaluate(mu).at(1);
  std::cout << "Total PDF momentum fraction: " << msr.Integrate(1e-7, 1) << std::endl;
  std::cout << "Up-valence PDF fraction: " << ( [] (double const& x) -> double { return 1 / x; } * ( DistPDFs.at(2) - DistPDFs.at(-2) ) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";
  for (int i = 0; i < (int) xiv.size(); i++)
    std::cout << "Up 1st moment (xi = " << xiv[i] << ")   : "
	      << ( [] (double const& x) -> double { return 1 / x; } * ( DistGPDs[i].at(2) - DistGPDs[i].at(-2) ) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";
  for (int i = 0; i < (int) xiv.size(); i++)
    std::cout << "Up 2nd moment (xi = " << xiv[i] << ")   : "
	      << ( DistGPDs[i].at(2) + DistGPDs[i].at(-2) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";
  for (int i = 0; i < (int) xiv.size(); i++)
    std::cout << "Up 3rd moment (xi = " << xiv[i] << ")   : "
	      << ( [] (double const& x) -> double { return x; } * ( DistGPDs[i].at(2) - DistGPDs[i].at(-2) ) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";
  for (int i = 0; i < (int) xiv.size(); i++)
    std::cout << "Up 4th moment (xi = " << xiv[i] << ")   : "
	      << ( [] (double const& x) -> double { return x * x; } * ( DistGPDs[i].at(2) + DistGPDs[i].at(-2) ) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";
  for (int i = 0; i < (int) xiv.size(); i++)
    std::cout << "Up 5rd moment (xi = " << xiv[i] << ")   : "
	      << ( [] (double const& x) -> double { return pow(x, 3); } * ( DistGPDs[i].at(2) - DistGPDs[i].at(-2) ) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";
  for (int i = 0; i < (int) xiv.size(); i++)
    std::cout << "Up 6th moment (xi = " << xiv[i] << ")   : "
	      << ( [] (double const& x) -> double { return pow(x, 4); } * ( DistGPDs[i].at(2) + DistGPDs[i].at(-2) ) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";

  return 0;
}
