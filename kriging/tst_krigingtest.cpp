#include <QString>
#include <QtTest>

#include <iostream>

#include "krig.h"
#include "powvargram.h"
#include <eigen3/Eigen/Dense>

/*
 * FIXME: is it a bad thing to call "using namespace" in a header?
 * It may have side effects on code including this header.
 */
using namespace Eigen;
using namespace std;

class KrigingTest : public QObject
{
    Q_OBJECT

public:
    KrigingTest();

private Q_SLOTS:
    void testCase1();
};

KrigingTest::KrigingTest()
{
}

void KrigingTest::testCase1()
{
    int npt = 100;
    MatrixXd x(npt, 1);
    VectorXd y(npt);

    x.col(0) = ArrayXd::LinSpaced(npt, 0, M_PI * 10);
    for (int i = 0; i < x.rows(); i++) {
        y(i) = std::cos(x(i, 0));
    }
    /*
     * FIXME: invalid use of incomplete type const class
     * Eigen::MatrixFunctionReturnValue...:
     * y = x.col(0).cos();
     *
     * This one works, but is redundant:
     * y = ArrayXd::LinSpaced(npt, 0, M_PI_2).cos();
     *
     * Therefore, we use explicit loop
     */

    cout << x << endl;
    cout << y << endl;

    /*
     * FIXME: how to validate the PowVargram implementation?
     */
    PowVargram vgram(x, y);
    for (int i = 0; i < x.rows(); i++) {
        cout << x(i, 0) << " " << vgram(x(i, 0)) << endl;
    }

    //Krig k(x, y, [&](double r) -> double { return vgram(r); }, nullptr);
}

QTEST_APPLESS_MAIN(KrigingTest)

#include "tst_krigingtest.moc"
