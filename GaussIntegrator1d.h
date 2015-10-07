#ifndef GAUSS_INTEGRATOR_1D_H
#define GAUSS_INTEGRATOR_1D_H

#include <iostream>
#include <cmath>

class GaussIntegrator1d {
public:
    GaussIntegrator1d(double a, double b);

    double getFirstPoint();

    double getSecondPoint();

    double calculate(double point1, double point2);

private:
    double* _points;
    double _h;
    double sqrt3 = 1/ std::sqrt(3);
};
#endif
