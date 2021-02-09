/*
 * Evolution.h
 *
 *  Created on: May 8, 2020
 *      Author: partons
 */

#ifndef EVOLUTION_EVOLUTION_H_
#define EVOLUTION_EVOLUTION_H_

#include <partons/beans/QuarkFlavor.h>
#include <partons/BaseObject.h>
#include <stddef.h>
#include <map>
#include <vector>

namespace PARTONS {

class Evolution: public BaseObject {

public:

    /**
     * Constructor.
     */
    Evolution();

    /**
     * Copy constructor.
     */
    Evolution(const Evolution& other);

    /**
     * Destructor.
     */
    virtual ~Evolution();

    /**
     * Perform the evolution.
     */
    void evolve2(std::vector<std::vector<double> >& g,
            std::map<QuarkFlavor::Type, std::vector<std::vector<double> > >& q,
            double muF20, double muF2) const;

private:

    /**
     * Evaluate lambdaQCD.
     */
    double lambdaQCD(size_t nf) const;

    /**
     * Evaluate beta0.
     */
    double beta0(size_t nf) const;

    /**
     * Evaluate alpha_s for given muF2.
     */
    double alphaS(double muF2, size_t nf) const;

    /**
     * Get number of active flavors for given muF2.
     */
    size_t getNActiveFlavors(double muF2) const;

    /**
     * Check input vectors.
     */
    void checkInputVectors(const std::vector<std::vector<double> >& v) const;

    std::vector<double> m_thresholds; ///< Thresholds.

    /**
     * Gegenbauer coefficient for C^(3/2)_{n}(x), k = 0 corresponds to the highest possible power of x.
     */
    double gegenbauerCoefficient32(size_t n, size_t k) const;

    /**
     * Gegenbauer coefficient for C^(5/2)_{n}(x), k = 0 corresponds to the highest possible power of x.
     */
    double gegenbauerCoefficient52(size_t n, size_t k) const;

    /**
     * Evaluate conformal moments from Mellin ones.
     */
    std::vector<std::vector<double> > mellinToConformal(
            const std::vector<std::vector<double> >& mell) const;

    /**
     * Evaluate Mellin moments from conformal ones.
     */
    std::vector<std::vector<double> > conformalToMellin(
            const std::vector<std::vector<double> >& conf) const;

    /**
     * Evaluate even conformal moments.
     */
    void evolveEven(size_t iMoment, std::vector<double>& mom, double muF20,
            double muF2, size_t nf) const;

    /**
     * Evaluate odd conformal moments.
     */
    void evolveOdd(size_t iMoment, std::vector<double>& gMom,
            std::map<QuarkFlavor::Type, std::vector<double> >& qMom,
            double muF20, double muF2, size_t nf) const;

    /**
     * STDOUT Print given moments.
     */
    void printMoments(const std::vector<std::vector<double> >& moms) const;
};

} /* namespace PARTONS */

#endif /* EVOLUTION_EVOLUTION_H_ */
