#include "argsolver.h"
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/pow.hpp>


struct equation_to_solve
{
    // Functor returning both 1st and 2nd derivatives.
    equation_to_solve(double a) : m_arg(a) {}

    std::tuple<double, double, double> operator()(double x)
    {
        using boost::math::pow;
        // Return both f(x) and f'(x) and f''(x).
        double fx = pow<2>(x) * (1 - x) - pow<2>(1 + 2 * x) * (1 + x / 2) * pow<2>(1 + x) * m_arg;
        double dx = (2 - 3 * x) * x - m_arg * (6.5 + x * (34 + x * (55.5 + x * (40 + x * 10))));
        double d2x = 2 - 6 * x - m_arg * (34 + x * (111 + x * 120 + x * 40));
        return std::make_tuple(fx, dx, d2x); // 'return' fx, dx and d2x.
    }

private:
    double m_arg; // to be 'fifth_rooted'.
};


// return -1 if no solution exists
// For one-level, a=2NB/M^2 *(kappa+1+2*log(N/M))
// For two-level, a=2B*sqrt(N/M)/M *(kappa+2+1.5*log(N/M))
double argsolver(double a)
{
    // find beta
    using namespace std;                // Help ADL of std functions.
    using namespace boost::math::tools; // for halley_iterate.
    double guess = 0.2;                 // Rough guess is to divide the exponent by five.
    double min = 0.0;                   // Minimum possible value is half our guess.
    double max = 1.0;                   // Maximum possible value is twice our guess.
    const int digits = 5;
    const boost::uintmax_t maxit = 10;
    boost::uintmax_t it = maxit;
    double result = halley_iterate(equation_to_solve(a), guess, min, max, digits, it);
    if (result < 1e-5 || result > 1 - 1e-5) // fail to solve
        result = -1;
    return result;
}