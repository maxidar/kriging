#ifndef KRIG_H
#define KRIG_H

#include <eigen3/Eigen/Dense>

using namespace Eigen;

typedef std::function<double (double)> Vargram;

class Krig
{
public:
    Krig(MatrixXd &xx, VectorXd &yy, Vargram vargram, const double *err = nullptr);
private:
    MatrixXd &x;
    Vargram vgram;
    int ndim;
    int npt;

    VectorXd dstar;
    VectorXd vstar;
    VectorXd yvi;
    VectorXd y;
    MatrixXd v;

    double lastval;
    double lasterr;
};

#endif // KRIG_H
