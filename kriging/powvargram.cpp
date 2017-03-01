#include "powvargram.h"
#include <math.h>

#include <QDebug>

inline double sqr(double x) { return x*x; }

/*
 * Based on Numerical Recipes, 3rd edition, p.147
 *
 * The model for the variogram is v(r) = apha*r^beta, where "a" is fitted
 * by unweighted least square over all data pairs.
 */

PowVargram::PowVargram(MatrixXd &x, VectorXd &y, const double beta,
                       const double nug) :
    bet(beta), nugsq(nug * nug)
{
    int npt = x.rows();
    int ndim = x.cols();
    double rb;
    double num = 0;
    double denom = 0;

    for (int i = 0; i < npt; i++) {
        for (int j = i + 1; j < npt; j++) {
            rb = 0;
            for (int k = 0; k < ndim; k++) {
                rb += sqr(x(i, k) - x(j, k));
            }
            rb = pow(rb, 0.5*beta);
            num += rb * (0.5 * sqr(y(i) - y(j)) - nugsq);
            denom += sqr(rb);
        }
    }
    alph = num/denom;
    qDebug() << "alpha" << alph;

}

double PowVargram::operator() (const double r) const
{
    return nugsq + alph * std::pow(r, bet);
}
