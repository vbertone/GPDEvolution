// APFEL++
#include <apfel/apfelxx.h>

double GhostUp(double const& m_x, double const& m_xi);
double GhostLO(double const& m_x, double const& m_xi);

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

  // Construct the GPD evolution objects and tabulate GPDs
  //const std::vector<double> xiv{0.1};
  const std::vector<double> xiv{0, 0.02, 0.05, 0.11, 0.2, 0.35, 0.5, 0.7, 0.9, 1};
  std::vector<apfel::TabulateObject<apfel::Set<apfel::Distribution>>> TabulatedGPDs;
  for(double const& xi : xiv)
    {
      // Input distributions (in the QCD evolution basis)
      const auto InPDFs = [&] (double const& x, double const&) -> std::map<int, double>
	{
	  // Call all functions only once.
	  const double upv  = 0;
	  const double dnv  = 0;
	  const double glu  = 0;
	  const double dbar = 0;
	  const double ubar = x * GhostUp(x, xi);
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
	      << ( DistGPDs[i].at(2) + DistGPDs[i].at(-2) + DistGPDs[i].at(-2) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";
  for (int i = 0; i < (int) xiv.size(); i++)
    std::cout << "Up-valence 3rd moment (xi = " << xiv[i] << ")   : "
	      << ( [] (double const& x) -> double { return x; } * ( DistGPDs[i].at(2) - DistGPDs[i].at(-2) ) ).Integrate(1e-7, 1)<< std::endl;
  std::cout << "\n";


  return 0;
}

double GhostUp(double const& m_x, double const& m_xi)
{
  double trans = 0.;
  if (fabs(m_xi) >= 1. || fabs(m_x) >= 1.)
    trans = 0.;
  else if(fabs(m_x) >= fabs(m_xi))
    trans = pow(m_xi, 2) * (pow(fabs(m_x), 2) - pow(m_xi, 2)) * pow(m_x - 1, 5) * (1.9781153654109946851*pow(10., -8)* m_x + 9.8905768270549734256*pow(10., -8)* pow(fabs(m_x), 2)  +  5.4207095851873525136*pow(10., -6)* pow(fabs(m_x), 3)  + 0.000026312301779772364694* pow(fabs(m_x), 4)  -  0.0021598622937353977786* pow(fabs(m_x), 5)  - 0.011008702138310537654* pow(fabs(m_x), 6)  -  0.56305855622681362802* pow(fabs(m_x), 7)  + 7.0072252445919459353* pow(fabs(m_x), 8)  -  15.538399226028916540* pow(fabs(m_x), 9)  - 44.371667827646332913* pow(fabs(m_x), 10)  +  181.94334698775348832* pow(fabs(m_x), 11)  - 97.011276613079764504* pow(fabs(m_x), 12)  -  181.72143366229912121* pow(fabs(m_x), 13)  + 53.025592131162241212* pow(fabs(m_x), 14)  +  158.14112838159099855* pow(fabs(m_x), 15)  + 36.688915854074584821* pow(fabs(m_x), 16)  -  56.114891535498807746* pow(fabs(m_x), 17)  - 39.657921691446077944* pow(fabs(m_x), 18)  -  7.9315843382892155887* pow(fabs(m_x), 19)  - 4.5799625973059142699*pow(10., -6)* m_x* pow(m_xi, 2)  -  0.000022899812986529571349* pow(fabs(m_x), 2) * pow(m_xi, 2)  +  0.0027154455327521476970* pow(fabs(m_x), 3) * pow(m_xi, 2)  + 0.013760426167652975056* pow(fabs(m_x), 4) * pow(m_xi, 2)  -  0.28836205489179660897* pow(fabs(m_x), 5) * pow(m_xi, 2)  + 10.591771592324523145* pow(fabs(m_x), 6) * pow(m_xi, 2)  -  98.315849249571176612* pow(fabs(m_x), 7) * pow(m_xi, 2)  + 155.79653880803078638* pow(fabs(m_x), 8) * pow(m_xi, 2)  +  1157.8242830169957934* pow(fabs(m_x), 9) * pow(m_xi, 2)  - 4741.3251010147422627* pow(fabs(m_x), 10) * pow(m_xi, 2)  +  5284.8979877142756907* pow(fabs(m_x), 11) * pow(m_xi, 2)  + 728.14359319075315680* pow(fabs(m_x), 12) * pow(m_xi, 2)  -  2983.4027503189065284* pow(fabs(m_x), 13) * pow(m_xi, 2)  - 1018.8906581020362782* pow(fabs(m_x), 14) * pow(m_xi, 2)  +  855.19510500029474521* pow(fabs(m_x), 15) * pow(m_xi, 2)  + 639.64983674072894688* pow(fabs(m_x), 16) * pow(m_xi, 2)  +  102.54889746562029949* pow(fabs(m_x), 17) * pow(m_xi, 2)  - 15.863168676578431177* pow(fabs(m_x), 18) * pow(m_xi, 2)  -  3.1726337353156862355* pow(fabs(m_x), 19) * pow(m_xi, 2)  - 0.00079628510480319772166* m_x* pow(m_xi, 4)  -  0.0039814255240159886083* pow(fabs(m_x), 2) * pow(m_xi, 4)  + 0.97810774434163170660* pow(fabs(m_x), 3) * pow(m_xi, 4)  -  18.553693906623295123* pow(fabs(m_x), 4) * pow(m_xi, 4)  + 80.019322701584110464* pow(fabs(m_x), 5) * pow(m_xi, 4)  +  205.57669342342535603* pow(fabs(m_x), 6) * pow(m_xi, 4)  - 1108.9720925829313412* pow(fabs(m_x), 7) * pow(m_xi, 4)  -  5292.4515355562200669* pow(fabs(m_x), 8) * pow(m_xi, 4)  + 32931.486866169982130* pow(fabs(m_x), 9) * pow(m_xi, 4)  -  55942.856316219247591* pow(fabs(m_x), 10) * pow(m_xi, 4)  + 25457.226013562047867* pow(fabs(m_x), 11) * pow(m_xi, 4)  +  14923.022390386173170* pow(fabs(m_x), 12) * pow(m_xi, 4)  - 5965.5382819086752892* pow(fabs(m_x), 13) * pow(m_xi, 4)  -  5267.5836235431551521* pow(fabs(m_x), 14) * pow(m_xi, 4)  - 681.29777879677809702* pow(fabs(m_x), 15) * pow(m_xi, 4)  +  231.40303918672976096* pow(fabs(m_x), 16) * pow(m_xi, 4)  + 44.870548399427869420* pow(fabs(m_x), 17) * pow(m_xi, 4)  -  0.88128714869880173208* pow(fabs(m_x), 18) * pow(m_xi, 4)  - 0.17625742973976034642* pow(fabs(m_x), 19) * pow(m_xi, 4)  -  0.24125260424052791049* m_x* pow(m_xi, 6)  + 1.9943886251995084972* pow(fabs(m_x), 2) * pow(m_xi, 6)  +  44.423991242391992362* pow(fabs(m_x), 3) * pow(m_xi, 6)  - 388.41593759720796164* pow(fabs(m_x), 4) * pow(m_xi, 6)  -  524.67033235850301142* pow(fabs(m_x), 5) * pow(m_xi, 6)  + 10202.070562178728021* pow(fabs(m_x), 6) * pow(m_xi, 6)  -  15249.732327475805994* pow(fabs(m_x), 7) * pow(m_xi, 6)  - 57095.368928101023144* pow(fabs(m_x), 8) * pow(m_xi, 6)  +  192110.31232352455849* pow(fabs(m_x), 9) * pow(m_xi, 6)  - 179121.24408611189640* pow(fabs(m_x), 10) * pow(m_xi, 6)  +  18406.489385218757264* pow(fabs(m_x), 11) * pow(m_xi, 6)  + 31455.322140932687386* pow(fabs(m_x), 12) * pow(m_xi, 6)  +  3515.3877982545091957* pow(fabs(m_x), 13) * pow(m_xi, 6)  - 1712.5462403270528480* pow(fabs(m_x), 14) * pow(m_xi, 6)  -  317.07878705916505185* pow(fabs(m_x), 15) * pow(m_xi, 6)  + 15.894038128903448588* pow(fabs(m_x), 16) * pow(m_xi, 6)  +  3.1788076257806897176* pow(fabs(m_x), 17) * pow(m_xi, 6)  + 0.97024584933205212617* pow(m_xi, 8)  -  22.821221936891537592* m_x* pow(m_xi, 8)  + 80.546663366834170180* pow(fabs(m_x), 2) * pow(m_xi, 8)  +  461.56501780904158235* pow(fabs(m_x), 3) * pow(m_xi, 8)  + 587.65889317797101130* pow(fabs(m_x), 4) * pow(m_xi, 8)  -  30942.846621894570688* pow(fabs(m_x), 5) * pow(m_xi, 8)  + 116058.28987373951330* pow(fabs(m_x), 6) * pow(m_xi, 8)  -  106136.24140963659717* pow(fabs(m_x), 7) * pow(m_xi, 8)  - 172560.74571146452891* pow(fabs(m_x), 8) * pow(m_xi, 8)  +  374932.30846467831723* pow(fabs(m_x), 9) * pow(m_xi, 8)  - 175075.44024096178211* pow(fabs(m_x), 10) * pow(m_xi, 8)  -  20654.346082834079890* pow(fabs(m_x), 11) * pow(m_xi, 8)  + 8810.6358164337153885* pow(fabs(m_x), 12) * pow(m_xi, 8)  +  1573.7524068122202844* pow(fabs(m_x), 13) * pow(m_xi, 8)  - 117.73422279657674585* pow(fabs(m_x), 14) * pow(m_xi, 8)  -  23.546844559315349169* pow(fabs(m_x), 15) * pow(m_xi, 8)  + 14.952833221954524282* pow(m_xi, 10)  -  56.315554285042311432* m_x* pow(m_xi, 10)  - 1594.9038256576631971* pow(fabs(m_x), 2) * pow(m_xi, 10)  +  9533.2068012653309691* pow(fabs(m_x), 3) * pow(m_xi, 10)  + 8797.1674504581463096* pow(fabs(m_x), 4) * pow(m_xi, 10)  -  164654.50613438123593* pow(fabs(m_x), 5) * pow(m_xi, 10)  + 378000.78902134162791* pow(fabs(m_x), 6) * pow(m_xi, 10)  -  217468.08809982818278* pow(fabs(m_x), 7) * pow(m_xi, 10)  - 217336.82472243109333* pow(fabs(m_x), 8) * pow(m_xi, 10)  +  254648.05826653967989* pow(fabs(m_x), 9) * pow(m_xi, 10)  - 38774.842040282214100* pow(fabs(m_x), 10) * pow(m_xi, 10)  -  6902.9304818604595745* pow(fabs(m_x), 11) * pow(m_xi, 10)  + 532.52370387248952846* pow(fabs(m_x), 12) * pow(m_xi, 10)  +  106.50474077449790569* pow(fabs(m_x), 13) * pow(m_xi, 10)  - 53.007998147416263552* pow(m_xi, 12)  +  2429.8528315988857087* m_x* pow(m_xi, 12)  - 21351.325462984932197* pow(fabs(m_x), 2) * pow(m_xi, 12)  +  60407.737431538303098* pow(fabs(m_x), 3) * pow(m_xi, 12)  + 3402.5646475321875629* pow(fabs(m_x), 4) * pow(m_xi, 12)  -  286936.74861146392456* pow(fabs(m_x), 5) * pow(m_xi, 12)  + 428790.42237029855502* pow(fabs(m_x), 6) * pow(m_xi, 12)  -  125978.51411127411373* pow(fabs(m_x), 7) * pow(m_xi, 12)  - 108872.51935779705046* pow(fabs(m_x), 8) * pow(m_xi, 12)  +  44740.517983802130789* pow(fabs(m_x), 9) * pow(m_xi, 12)  - 1859.3641735903611877* pow(fabs(m_x), 10) * pow(m_xi, 12)  -  371.87283471807223753* pow(fabs(m_x), 11) * pow(m_xi, 12)  - 853.29844833922602184* pow(m_xi, 14)  +  13236.674848582701410* m_x* pow(m_xi, 14)  - 65719.669334875971679* pow(fabs(m_x), 2) * pow(m_xi, 14)  +  119859.88542902686035* pow(fabs(m_x), 3) * pow(m_xi, 14)  - 10652.576069538706354* pow(fabs(m_x), 4) * pow(m_xi, 14)  -  194199.58887279685824* pow(fabs(m_x), 5) * pow(m_xi, 14)  + 164028.73526016797471* pow(fabs(m_x), 6) * pow(m_xi, 14)  -  10607.750816680038307* pow(fabs(m_x), 7) * pow(m_xi, 14)  - 11848.731198908898554* pow(fabs(m_x), 8) * pow(m_xi, 14)  +  1601.1111621099985052* pow(fabs(m_x), 9) * pow(m_xi, 14)  - 2298.8021027316481529* pow(m_xi, 16)  +  22302.756908501618752* m_x* pow(m_xi, 16)  - 71566.303830252158958* pow(fabs(m_x), 2) * pow(m_xi, 16)  +  82286.518749752046173* pow(fabs(m_x), 3) * pow(m_xi, 16)  - 442.91596095874039645* pow(fabs(m_x), 4) * pow(m_xi, 16)  -  47196.680019254635270* pow(fabs(m_x), 5) * pow(m_xi, 16)  + 13994.198385616129654* pow(fabs(m_x), 6) * pow(m_xi, 16)  +  39.220299591306186264* pow(fabs(m_x), 7) * pow(m_xi, 16)  - 145.98740448131537559* pow(fabs(m_x), 8) * pow(m_xi, 16)  -  2275.5125529794829239* pow(m_xi, 18)  + 14444.913950371349354* m_x* pow(m_xi, 18)  -  28770.250505356937508* pow(fabs(m_x), 2) * pow(m_xi, 18)  + 17387.476362618452112* pow(fabs(m_x), 3) * pow(m_xi, 18)  +  2446.0343387147832737* pow(fabs(m_x), 4) * pow(m_xi, 18)  - 2021.5120025470624955* pow(fabs(m_x), 5) * pow(m_xi, 18)  +  134.62595327551924933* pow(fabs(m_x), 6) * pow(m_xi, 18)  - 913.10064154886320225* pow(m_xi, 20)  +  3527.0113185998332592* m_x* pow(m_xi, 20)  - 3564.3766684119398832* pow(fabs(m_x), 2) * pow(m_xi, 20)  +  495.42631216184088147* pow(fabs(m_x), 3) * pow(m_xi, 20)  + 51.307015970142626343* pow(fabs(m_x), 4) * pow(m_xi, 20)  -  139.70810880715583728* pow(m_xi, 22)  + 252.57036625594730332* m_x* pow(m_xi, 22)  -  39.456318679884862690* pow(fabs(m_x), 2) * pow(m_xi, 22)  - 6.1171615640755502789* pow(m_xi, 24)) / pow(1 - pow(m_xi, 2), 17);
  else
    trans = m_x * m_xi * (pow(m_x, 2) - pow(m_xi, 2)) * (-1.7773666992943329387 * pow(m_x, 8) - 10.370510860523957791 * pow(m_x, 10) +  145.98740448131537559 * pow(m_x, 12) - 1.9781153654109946851*pow(10., -8) * m_xi -  5.1239922803757033108*pow(10., -6) * pow(m_x, 2) * m_xi + 0.0022381068586968210246 * pow(m_x, 4) * m_xi +  0.52984978684825497497 * pow(m_x, 6) * m_xi + 25.890615725640473924 * pow(m_x, 8) * m_xi -  351.83984196968805298 * pow(m_x, 10) * m_xi - 3.3627961211986909647*pow(10., -7) * pow(m_xi, 2) -  0.000087107868766386956284 * pow(m_x, 2) * pow(m_xi, 2) +  0.038047816597845957418 * pow(m_x, 4) * pow(m_xi, 2) - 10.335622269227429726 * pow(m_x, 6) * pow(m_xi, 2) +  261.21121208806439382 * pow(m_x, 8) * pow(m_xi, 2) - 134.62595327551924933 * pow(m_x, 10) * pow(m_xi, 2) +  1.5534460882270924017*pow(10., -6) * pow(m_xi, 3) - 0.0035681157906092190176 * pow(m_x, 2) * pow(m_xi, 3) +  0.67223398146627485790 * pow(m_x, 4) * pow(m_xi, 3) - 93.482817693917441024 * pow(m_x, 6) * pow(m_xi, 3) +  422.54770143258561504 * pow(m_x, 8) * pow(m_xi, 3) + 0.000058691426263368004089 * pow(m_xi, 4) -  0.052295613038783575496 * pow(m_x, 2) * pow(m_xi, 4) + 19.951779424585110861 * pow(m_x, 4) * pow(m_xi, 4) -  329.38563291057568950 * pow(m_x, 6) * pow(m_xi, 4) - 51.307015970142626343 * pow(m_x, 8) * pow(m_xi, 4) +  0.0014011796927368399125 * pow(m_xi, 5) - 1.4408519441839956259 * pow(m_x, 2) * pow(m_xi, 5) +  85.697547193537273931 * pow(m_x, 4) * pow(m_xi, 5) - 28.927315409761905341 * pow(m_x, 6) * pow(m_xi, 5) +  0.017572303842736308887 * pow(m_xi, 6) - 9.1775269071209085712 * pow(m_x, 2) * pow(m_xi, 6) +  74.584872853564762301 * pow(m_x, 4) * pow(m_xi, 6) + 39.456318679884862690 * pow(m_x, 6) * pow(m_xi, 6) +  0.38379821284177021108 * pow(m_xi, 7) - 19.646125334208574221 * pow(m_x, 2) * pow(m_xi, 7) -  52.388487136347140032 * pow(m_x, 4) * pow(m_xi, 7) + 0.87693289090585982442 * pow(m_xi, 8) +  14.045782444018862029 * pow(m_x, 2) * pow(m_xi, 8) + 6.1171615640755502789 * pow(m_x, 4) * pow(m_xi, 8) -  0.66398966230367550442 * pow(m_xi, 9) - 2.7721078351021553388 * pow(m_x, 2) * pow(m_xi, 9) +  0.19188594448238966661 * pow(m_xi, 10) + 0.011287408498964098036 * pow(m_xi, 11)) / pow(1 + m_xi, 17);

  if (m_x < -fabs(m_xi))
    trans = -trans;

  return trans;
}

double GhostLO(double const& m_x, double const& m_xi)
{
  double trans = 0.;
  if (fabs(m_xi) >= 1. || fabs(m_x) >= 1.)
    trans = 0.;
  else if(fabs(m_x) >= fabs(m_xi))
    trans = - 1024. * pow(1 - m_x, 3) * (pow(m_x, 2) - pow(m_xi, 2)) * pow(m_xi, 2) * (3. * m_x - 12. * pow(m_x, 2) + 4. * pow(m_x, 3) + 9. * pow(m_x, 4) + 3. * pow(m_x, 5) - 7. * pow(m_xi, 2) + 42. * m_x * pow(m_xi, 2) - 42. * pow(m_x, 2) * pow(m_xi, 2) - 14. * pow(m_x, 3) * pow(m_xi, 2) - 21. * pow(m_xi, 4) + 42. * m_x * pow(m_xi, 4) - 7. * pow(m_xi, 6)) / 7. / pow(1 - pow(m_xi, 2), 7);
  else
    trans = 1024. * m_x * m_xi * (pow(m_x, 2) - pow(m_xi, 2)) * (7. * pow(m_x, 2) - 3. * m_xi) / 7. / pow(1 + m_xi, 7);

  if (m_x < -fabs(m_xi))
    trans = -trans;

  return trans;
}
