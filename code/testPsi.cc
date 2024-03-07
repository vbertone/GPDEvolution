// APFEL++
#include <apfel/apfelxx.h>

// Boost
#include <boost/math/special_functions/digamma.hpp>

int main()
{
  const int nx = 100;
  const double xmin = 1e-5;
  const double xmax = 10;
  const double xstp = exp(log(xmax / xmin) / ( nx - 1 ));
  for (double x = xmin; x <= 1.0000001 * xmax; x *= xstp)
    std::cout << std::scientific << x << "  " << apfel::digamma(x) << "  " << boost::math::digamma(x) << "  " << apfel::digamma(x) / boost::math::digamma(x) << std::endl;

  return 0;
}
