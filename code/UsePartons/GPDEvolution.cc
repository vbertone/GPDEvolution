// APFEL++
#include <apfel/apfelxx.h>

// LHAPDF libs
#include <LHAPDF/LHAPDF.h>

// PARTONS
#include <partons/Partons.h>
#include <ElementaryUtils/logger/LoggerManager.h>
#include <partons/ServiceObjectRegistry.h>
#include <partons/beans/KinematicUtils.h>
#include <partons/modules/gpd/GPDGK19.h>
#include <partons/ModuleObjectFactory.h>
#include <partons/services/GPDService.h>

int main(int argc, char** argv)
{
  // Final scale
  const double mu = 10;

  // Evolution parameters
  const int    pto   = 0;
  const double Qref  = 91.1876;
  const double asref = 0.135;
  const double mc    = 1.5;
  const double mb    = 4.75;
  const double mt    = 175;
  const double mu0   = sqrt(2);
  const double t     = -0.1;

  // Vector of thresholds
  const std::vector<double> Thresholds = {0, 0, 0, mc, mb, mt};

  // Running coupling
  apfel::AlphaQCD a{asref, Qref, Thresholds, pto};
  const  apfel::TabulateObject<double> Alphas{a, 100, 0.9, 100, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100, 1e-7, 3}, apfel::SubGrid{200, 1e-1, 3}, apfel::SubGrid{100, 9e-1, 3}}};

  // Initialise PARTONS
  PARTONS::Partons* pPartons = PARTONS::Partons::getInstance();
  pPartons->init(argc, argv);

  // Retrieve GPD service
  PARTONS::GPDService* pGPDService = PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

  // Create GPD module with the BaseModuleFactory
  PARTONS::GPDModule* pGPDModel = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(PARTONS::GPDGK16::classId);

  // Construct the GPD evolution objects and tabulate GPDs
  //const std::vector<double> xiv{0, 0.02, 0.05, 0.11, 0.2, 0.35, 0.5, 0.7, 0.9, 1};
  const std::vector<double> xiv{0.00001, 0.02, 0.05, 0.11, 0.2, 0.35, 0.5, 0.7, 0.9, 1};
  std::vector<apfel::TabulateObject<apfel::Set<apfel::Distribution>>> TabulatedGPDs;
  for(double const& xi : xiv)
    {
      // Input distributions (in the QCD evolution basis)
      const auto InPDFs = [=] (double const& x, double const& Q) -> std::map<int, double>
	{
	  PARTONS::GPDResult gpdResult = pGPDService->computeSingleKinematic(PARTONS::GPDKinematic{x, xi, t, Q * Q, Q * Q}, pGPDModel);
	  std::map<int, double> PhysMap{{0, 0}};
	  PhysMap[-6] = 0;
	  PhysMap[-5] = 0;
	  PhysMap[-4] = 0;
	  PhysMap[-3] = x * ( gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE).getQuarkDistributionPlus() -
			      gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE).getQuarkDistributionMinus() ) / 2;
	  PhysMap[-2] = x * ( gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistributionPlus() -
			      gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistributionMinus() ) / 2;
	  PhysMap[-1] = x * ( gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN).getQuarkDistributionPlus() -
			      gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN).getQuarkDistributionMinus() ) / 2;
	  PhysMap[0] = x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getGluonDistribution().getGluonDistribution();
	  PhysMap[1] = x * ( gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN).getQuarkDistributionPlus() +
			     gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN).getQuarkDistributionMinus() ) / 2;
	  PhysMap[2] = x * ( gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistributionPlus() +
			     gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistributionMinus() ) / 2;
	  PhysMap[3] = x * ( gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE).getQuarkDistributionPlus() +
			     gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE).getQuarkDistributionMinus() ) / 2;
	  PhysMap[4] = 0;
	  PhysMap[5] = 0;
	  PhysMap[6] = 0;
	  return apfel::PhysToQCDEv(PhysMap);
	};
      TabulatedGPDs.push_back(apfel::TabulateObject<apfel::Set<apfel::Distribution>> {*(BuildDglap(InitializeGpdObjects(g, Thresholds, xi), InPDFs, mu0, pto, as)), 30, 1, 100, 3});
    }

  // Check polynomiality
  std::vector<std::map<int, apfel::Distribution>> DistGPDs;
  for (auto const& tgpd : TabulatedGPDs)
    DistGPDs.push_back(apfel::QCDEvToPhys(tgpd.Evaluate(mu).GetObjects()));

  std::cout << std::scientific << "\n";
  for (int i = 0; i < (int) xiv.size(); i++)
    std::cout << "Up-valence 1st moment (xi = " << xiv[i] << ")   : "
	      << ( [] (double const& x) -> double { return 1 / x; } * ( DistGPDs[i].at(2) - DistGPDs[i].at(-2) ) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";
  for (int i = 0; i < (int) xiv.size(); i++)
    std::cout << "Up-valence 2st moment (xi = " << xiv[i] << ")   : "
	      << ( DistGPDs[i].at(2) + DistGPDs[i].at(-2) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";
  for (int i = 0; i < (int) xiv.size(); i++)
    std::cout << "Up-valence 3rd moment (xi = " << xiv[i] << ")   : "
	      << ( [] (double const& x) -> double { return x; } * ( DistGPDs[i].at(2) - DistGPDs[i].at(-2) ) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";
  for (int i = 0; i < (int) xiv.size(); i++)
    std::cout << "Up-valence 4th moment (xi = " << xiv[i] << ")   : "
	      << ( [] (double const& x) -> double { return x * x; } * ( DistGPDs[i].at(2) + DistGPDs[i].at(-2) ) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";

  // Remove pointer references
  // Module pointers are managed by PARTONS
  PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(pGPDModel, 0);
  pGPDModel = 0;

  // Close PARTONS application properly
  if (pPartons)
    pPartons->close();

  return 0;
}
