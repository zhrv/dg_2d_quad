#ifndef GAUSSINTEGRATOR2D_H
#define GAUSSINTEGRATOR2D_H

#include <cmath>

class GaussIntegrator2d {
public:
    GaussIntegrator2d(double a, double b, double c, double d);

    GaussIntegrator2d(double *bounds);

    ~GaussIntegrator2d();

    double** getPoints();

    double calculate(double* values);

    const int N = 4;
private:
    double sqrt3 = 1/sqrt(3);
    double g_points[4][2] = { { -sqrt3, -sqrt3 },
                              { -sqrt3,  sqrt3 },
                              {  sqrt3, -sqrt3 },
                              {  sqrt3,  sqrt3 } };
    double** _points;
    double _h[2];

    void init(double* bounds);
};


#endif
