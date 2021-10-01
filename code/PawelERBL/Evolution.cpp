/*
 * Evolution.cpp
 *
 *  Created on: May 8, 2020
 *      Author: partons
 */

#include "Evolution.h"

#include <ElementaryUtils/logger/CustomException.h>
#include <ElementaryUtils/string_utils/Formatter.h>
#include <cmath>
#include <cstdlib>
#include <iostream> // Delete this later
#include <iterator>
#include <limits>
#include <utility>

namespace PARTONS {

  Evolution::Evolution() :
    BaseObject("Evolution")
  {
    //thresholds
    m_thresholds.push_back(0.);                 //0
    m_thresholds.push_back(pow(2.2 / 1.E3, 2)); //u
    m_thresholds.push_back(pow(4.7 / 1.E3, 2)); //d
    m_thresholds.push_back(pow(96. / 1.E3, 2)); //s
    m_thresholds.push_back(pow(1.4, 2));        //c
    m_thresholds.push_back(pow(4.75, 2));       //b
    m_thresholds.push_back(pow(175, 2));        //t
    m_thresholds.push_back(std::numeric_limits<double>::max()); //inf
  }

  Evolution::Evolution(const Evolution& other) :
    BaseObject(other)
  {
    //thresholds
    m_thresholds = other.m_thresholds;
  }

  Evolution::~Evolution()
  {
  }

  double Evolution::gegenbauerCoefficient32(size_t n, size_t k) const
  {
    double result = 1.;

    //c_n_0
    if (k == 0) {

      for (int i = 1; i <= (n - (n % 2)) / 2; i++) 
	{
	  int a = n - ((n + 1) % 2) + 2 * i;
	  int b = 2 * i;
	  b *= 2;
	  result *= double(a) / double(b);
	}
      result *= (2 * n + 1) * pow(2, (n - (n % 2)) / 2);
    }

    //c_n_k
    else
      {
	result = gegenbauerCoefficient32(n, k - 1);
	result *= (-2. + 2 * k - n) * (-1. + 2 * k - n) / (2 * k * (-3. + 2 * k - 2 * n));
      }
    return result;
  }

  double Evolution::gegenbauerCoefficient52(size_t n, size_t k) const
  {
    return (3. - 2 * k + 2 * n) / 3. * gegenbauerCoefficient32(n, k);
  }

  void Evolution::evolve2(std::vector<std::vector<double>>& gMell,
			  std::map<PARTONS::QuarkFlavor::Type,
			  std::vector<std::vector<double>>>& qMell,
			  double muF20,
			  double muF2) const
  {
    //number of flavors
    size_t nf = qMell.size();
    // nf = 3;

    //iterator
    std::map<PARTONS::QuarkFlavor::Type, std::vector<std::vector<double>>>::iterator itQ;

    //check input vectors
    checkInputVectors(gMell);

    size_t nMoments = gMell.size();

    for (itQ = qMell.begin(); itQ != qMell.end(); itQ++)
      {
	checkInputVectors(itQ->second);

	if (itQ->second.size() != nMoments)
	  throw ElemUtils::CustomException(getClassName(), __func__, "Inconsistent input");
      }

    //evaluate conformal moments
    std::vector<std::vector<double> > gConf = mellinToConformal(gMell);

    std::map<PARTONS::QuarkFlavor::Type, std::vector<std::vector<double>>> qConf;

    for (itQ = qMell.begin(); itQ != qMell.end(); itQ++) {
      qConf.insert(
		   std::make_pair(itQ->first, mellinToConformal(itQ->second)));
    }

    //evolve
    for (size_t iMoment = 0; iMoment < nMoments; iMoment++) {

      //odd moments
      if (iMoment % 2) {

	std::map<PARTONS::QuarkFlavor::Type, std::vector<double> > q;

	for (itQ = qConf.begin(); itQ != qConf.end(); itQ++) {
	  q.insert(std::make_pair(itQ->first, itQ->second.at(iMoment)));
	}

	evolveOdd(iMoment, gConf.at(iMoment), q, muF20, muF2, nf);

	for (itQ = qConf.begin(); itQ != qConf.end(); itQ++) {
	  itQ->second.at(iMoment) = q.find(itQ->first)->second;
	}
      }

      //even moments
      else {

	//quarks
	for (itQ = qConf.begin(); itQ != qConf.end(); itQ++) {
	  evolveEven(iMoment, itQ->second.at(iMoment), muF20, muF2, nf);
	}
      }
    }

    //return
    gMell = conformalToMellin(gConf);

    qMell.clear();

    for (itQ = qConf.begin(); itQ != qConf.end(); itQ++) {
      qMell.insert(
		   std::make_pair(itQ->first, conformalToMellin(itQ->second)));
    }
  }

  void Evolution::evolveEven(size_t n, std::vector<double>& mom, double muF20,
			     double muF2, size_t nf) const {

    //CF
    const double CF = 4. / 3.;

    //gamma
    double gamma;

    gamma = 0.5 - 1. / (n + 1.) / (n + 2.);
    for (size_t j = 2; j <= n + 1; j++) {
      gamma += 2. / double(j);
    }
    gamma *= CF;

    //loop
    std::vector<double>::iterator it;

    for (it = mom.begin(); it != mom.end(); it++) {
      (*it) *= pow(alphaS(muF2, nf) / alphaS(muF20, nf),
		   2 * gamma / beta0(nf));
    }
  }

  void Evolution::evolveOdd(size_t n, std::vector<double>& gMom,
			    std::map<PARTONS::QuarkFlavor::Type, std::vector<double> >& qMom, double muF20,
			    double muF2, size_t nf) const {

    //check consistency
    if (qMom.size() != nf) {
      throw ElemUtils::CustomException(getClassName(), __func__,
				       ElemUtils::Formatter() << "Inconsistent input, map size: "
				       << qMom.size() << ", nf: " << nf);
    }

    //iterator
    std::map<PARTONS::QuarkFlavor::Type, std::vector<double> >::iterator itQ;

    //CF
    const double CF = 4. / 3.;

    //TF
    const double TF = 1. / 2.;

    //CA
    const double CA = 3.;

    //gamma
    double sum = 0.;

    for (size_t k = 2; k <= (n + 1); k++) {
      sum += 1. / k;
    }

    double gammaQQ = CF * (0.5 - (1. / ((n + 1.) * (n + 2.))) + 2 * sum);
    double gammaQG = -1. * nf * TF * (n * n + 3 * n + 4.)
      / (n * (n + 1.) * (n + 2.));
    double gammaGQ = -2. * CF * (n * n + 3 * n + 4.)
      / ((n + 1.) * (n + 2.) * (n + 3.));
    double gammaGG =
      CA
      * (1. / 6. - 2. / (n * (n + 1)) - 2. / ((n + 2) * (n + 3))
	 + 2 * sum) + (2. / 3.) * nf * TF;

    double gammaPlus = 0.5
      * (gammaQQ + gammaGG
	 + sqrt(pow(gammaQQ - gammaGG, 2) + 4 * gammaQG * gammaGQ));

    double gammaMinus = 0.5
      * (gammaQQ + gammaGG
	 - sqrt(pow(gammaQQ - gammaGG, 2) + 4 * gammaQG * gammaGQ));

    //mixing coefficients
    double aPlus = 2 * nf / double(n) * (gammaPlus - gammaQQ) / gammaQG;
    double aMinus = 2 * nf / double(n) * (gammaMinus - gammaQQ) / gammaQG;

    //== singlet ================================================================
    std::vector<std::vector<double> > qMomSinglet;

    for (itQ = qMom.begin(); itQ != qMom.end(); itQ++) {

      if (itQ == qMom.begin())
	continue;

      std::vector<double> qMomSingletThis;

      for (size_t i = 0; i < itQ->second.size(); i++) {
	qMomSingletThis.push_back(
				  qMom.begin()->second.at(i) - itQ->second.at(i));
      }

      evolveEven(n, qMomSingletThis, muF20, muF2, nf);

      qMomSinglet.push_back(qMomSingletThis);
    }

    //== non-singlet ============================================================

    //sum over quark flavors
    std::vector<double> qMomNonSinglet(qMom.begin()->second.size(), 0.);

    for (itQ = qMom.begin(); itQ != qMom.end(); itQ++) {
      for (size_t i = 0; i < itQ->second.size(); i++) {
	qMomNonSinglet.at(i) += itQ->second.at(i);
      }
    }

    for (size_t i = 0; i < qMomNonSinglet.size(); i++) {

      //quark-gluon combinations
      double qgPlus =
	(gMom.at(i) - aMinus * qMomNonSinglet.at(i) / double(nf))
	/ (aPlus - aMinus);

      double qgMinus = -1
	* (gMom.at(i) - aPlus * qMomNonSinglet.at(i) / double(nf))
	/ (aPlus - aMinus);

      //evolve
      qgPlus *= pow(alphaS(muF2, nf) / alphaS(muF20, nf),
		    2 * gammaPlus / beta0(nf));

      qgMinus *= pow(alphaS(muF2, nf) / alphaS(muF20, nf),
		     2 * gammaMinus / beta0(nf));

      //back
      qMomNonSinglet.at(i) = nf * (qgPlus + qgMinus);
      gMom.at(i) = aPlus * qgPlus + aMinus * qgMinus;
    }

    //== back ================================================================

    size_t a = 0;

    for (itQ = qMom.begin(); itQ != qMom.end(); itQ++) {
      for (size_t i = 0; i < itQ->second.size(); i++) {

	if (itQ == qMom.begin()) {

	  itQ->second.at(i) = qMomNonSinglet.at(i);

	  for (size_t j = 0; j < qMomSinglet.size(); j++) {
	    itQ->second.at(i) += qMomSinglet.at(j).at(i);
	  }

	  itQ->second.at(i) /= double(nf);
	} else {
	  itQ->second.at(i) = qMom.begin()->second.at(i)
	    - qMomSinglet.at(a-1).at(i);

	}
      }
      a++;
    }
  }

  double Evolution::lambdaQCD(size_t nf) const {
    return 0.22;
  }

  double Evolution::beta0(size_t nf) const {
    return 11. - 2. * nf / 3.;
  }

  double Evolution::alphaS(double muF2, size_t nf) const {
    return 4 * M_PI / (beta0(nf) * log(muF2 / (pow(lambdaQCD(nf), 2))));
  }

  size_t Evolution::getNActiveFlavors(double muF2) const {

    //iterator
    std::vector<double>::const_iterator it;

    for (it = m_thresholds.begin() + 1; it != m_thresholds.end() - 1; it++) {
      if (muF2 >= (*it) && muF2 < (*(it + 1)))
	return it - m_thresholds.begin();
    }

    throw ElemUtils::CustomException(getClassName(), __func__,
				     ElemUtils::Formatter()
				     << "Unable to identify number of active flavors for muF2 = "
				     << muF2);

    return 0;
  }

  void Evolution::checkInputVectors(const std::vector<std::vector<double> >& v) const {

    //iterator
    std::vector<std::vector<double> >::const_iterator it;

    //loop
    for (it = v.begin(); it != v.end(); it++) {

      //which Mellin moment?
      size_t n = (it - v.begin());

      //number of expected Mellin coefficients.
      //n:  0 1 2 3 4 5 6
      //nI: 1 2 2 3 3 4 4
      size_t nI = 1 + (n + (n % 2)) / 2;

      //check number of elements
      if (it->size() != nI) {
	throw ElemUtils::CustomException(getClassName(), __func__,
					 ElemUtils::Formatter()
					 << "Wrong number of Mellin coefficients for moment: "
					 << n << ", is: " << it->size() << ", expected: "
					 << nI);
      }
    }
  }

  std::vector<std::vector<double>> Evolution::mellinToConformal(const std::vector<std::vector<double>>& mell) const
  {
    //the same size
    std::vector<std::vector<double> > conf = mell;

    //zero
    for (size_t n = 0; n < conf.size(); n++) {
      for (size_t k = 0; k < conf.at(n).size(); k++) {
	conf[n][k] = 0.;
      }
    }

    //set
    for (size_t n = 0; n < conf.size(); n++) {
      for (size_t k = 0; k <= (n - (n % 2)) / 2; k++) {

	double c_nk = gegenbauerCoefficient32(n, k);

	for (size_t j = 0; j <= (n + (n % 2)) / 2 - k; j++) {
	  conf.at(n).at(k + j) += c_nk * mell.at(n - 2 * k).at(j);
	}
      }
    }

    //return
    return conf;
  }

  std::vector<std::vector<double> > Evolution::conformalToMellin(const std::vector<std::vector<double> >& conf) const {

    //the same
    std::vector<std::vector<double> > mell = conf;

    //set
    for (size_t n = 0; n < conf.size(); n++) {
      for (size_t k = 1; k <= (n - (n % 2)) / 2; k++) {

	double c_nk = gegenbauerCoefficient32(n, k);

	for (size_t j = 0; j <= (n + (n % 2)) / 2 - k; j++) {
	  mell.at(n).at(k + j) -= c_nk * mell.at(n - 2 * k).at(j);
	}
      }

      double c_n0 = gegenbauerCoefficient32(n, 0);

      for (size_t k = 0; k <= (n + (n % 2)) / 2; k++) {
	mell.at(n).at(k) /= c_n0;
      }
    }

    //return
    return mell;
  }

  void Evolution::printMoments(const std::vector<std::vector<double> >& moms) const
  {
    for (size_t n = 0; n < moms.size(); n++)
      for (size_t k = 0; k < moms.at(n).size(); k++)
	std::cout << __func__ << ": " << n << '\t' << k << '\t' << moms[n][k] << std::endl;
  }
}

std::vector<std::vector<double>> setMellinMoments(const std::vector<double>& v) {

  //iterator
  std::vector<double>::const_iterator it;

  //Mellin moment
  size_t n = 0;

  //Mellin coefficient
  size_t i = 0;

  //number of Mellin coefficients for given moment
  size_t nI;

  //vector
  std::vector<double> nVec;

  //result
  std::vector<std::vector<double> > result;

  //loop
  for (it = v.begin(); it != v.end(); it++) {

    if (i == 0) {

      nI = 1 + (n + (n % 2)) / 2;
      nVec = std::vector<double>(nI);
    }

    nVec.at(i) = (*it);

    i++;

    if (i == nI) {

      result.push_back(nVec);
      i = 0;
      n++;
    }
  }

  //check number of elements
  if (i != 0)
    throw ElemUtils::CustomException("main", __func__, "Wrong number of Mellin coefficients");

  //return
  return result;
}

int main()
{
  //quarks (u, d, s)
  //for u we set 3 first moments, where:
  //A0,0 = 1			        int dx x^0 H(x, xi) = A0,0
  //A1,0 = 2	A1,1 = 2.1		int dx x^1 H(x, xi) = A1,0 + A1,1 xi^2
  //A2,0 = 3	A2,1 = 3.1		int dx x^2 H(x, xi) = A2,0 + A2,1 xi^2
  std::map<PARTONS::QuarkFlavor::Type, std::vector<std::vector<double>>> q;
  q.insert(std::make_pair(PARTONS::QuarkFlavor::UP,      setMellinMoments( { 1., 2., 2.1, 3., 3.1 })));
  q.insert(std::make_pair(PARTONS::QuarkFlavor::DOWN,    setMellinMoments( { 0., 0., 0.,  0., 0. })));
  q.insert(std::make_pair(PARTONS::QuarkFlavor::STRANGE, setMellinMoments( { 0., 0., 0.,  0., 0. })));

  //gluons
  std::vector<std::vector<double> > g = setMellinMoments( { 1., 2., 2.2, 3., 3.2 });

  //evolve from 1 GeV2 to 2 GeV2
  PARTONS::Evolution evolution;
  evolution.evolve2(g, q, 1., 2.);

  //print
  std::vector<std::vector<double> >::const_iterator it0;
  std::vector<double>::const_iterator it1;

  std::cout << "up: " << std::endl;

  for (it0 = q.find(PARTONS::QuarkFlavor::UP)->second.begin();
       it0 != q.find(PARTONS::QuarkFlavor::UP)->second.end(); it0++) {
    for (it1 = it0->begin(); it1 != it0->end(); it1++) {
      std::cout << (*it1) << std::endl;
    }
  }

  std::cout << "down (and strange): " << std::endl;

  for (it0 = q.find(PARTONS::QuarkFlavor::DOWN)->second.begin();
       it0 != q.find(PARTONS::QuarkFlavor::DOWN)->second.end(); it0++) {
    for (it1 = it0->begin(); it1 != it0->end(); it1++) {
      std::cout << (*it1) << std::endl;
    }
  }

  std::cout << "gluon: " << std::endl;

  for (it0 = g.begin(); it0 != g.end(); it0++) {
    for (it1 = it0->begin(); it1 != it0->end(); it1++) {
      std::cout << (*it1) << std::endl;
    }
  }
  return 0;
}
