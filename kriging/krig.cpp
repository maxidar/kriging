#include "krig.h"

Krig::Krig(MatrixXd &xx, VectorXd &yy, Vargram vargram, const double *err) :
    x(xx), vgram(vargram), ndim(xx.rows()), npt(xx.cols()), dstar(npt+1),
    vstar(npt+1), yvi(npt+1), y(npt+1), v(npt+1,npt+1)
{
    (void) err;
    (void) yy;

    int i, j;
    for (int i = 0; i < npt; i++) {
        y[i] = yy[i];
    }
}

