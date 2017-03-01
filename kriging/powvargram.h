#ifndef POWVARGRAM_H
#define POWVARGRAM_H

#include <eigen3/Eigen/Dense>

using namespace Eigen;

class PowVargram
{
public:
    PowVargram(MatrixXd &x, VectorXd &y, const double beta = 1.5,
               const double nug = 0.);
    double operator() (const double r) const;
private:
    double alph;
    const double bet;
    const double nugsq;
};

#endif // POWVARGRAM_H
