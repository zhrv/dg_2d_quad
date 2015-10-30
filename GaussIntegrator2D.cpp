#include "GaussIntegrator2D.h"

GaussIntegrator2d::GaussIntegrator2d(double* bounds) {
    init(bounds);
}

double **GaussIntegrator2d::getPoints() {
    return _points;
}

GaussIntegrator2d::GaussIntegrator2d(double a, double b, double c, double d) {
    double array[4] = {a, b, c, d};
    init(array);
}

void GaussIntegrator2d::init(double *bounds) {
    double c[2];

    for (int i = 0; i < 2; ++i) {
        double a = bounds[i*2];
        double b = bounds[i*2+1];

        c[i] = (a+b)/2;
        _h[i] = (b-a)/2;
    }

    _points = new double*[N];
    for (int i = 0; i < N; ++i) {
        _points[i] = new double[2];
        for (int j = 0; j < 2; ++j) {
            _points[i][j] = c[j] + _h[j]*g_points[i][j];
        }
    }
}

GaussIntegrator2d::~GaussIntegrator2d() {
    for (int i = 0; i < N; i++)
        delete[] _points[i];
    delete[] _points;
}

double GaussIntegrator2d::calculate(double *values) {
    double result = 0;

    for (int i = 0; i < N; ++i) {
        result += values[i];
    }

    result *= _h[0]*_h[1];

    return result;
}
