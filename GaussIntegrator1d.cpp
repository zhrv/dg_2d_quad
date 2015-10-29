#include "GaussIntegrator1d.h"

GaussIntegrator1d::GaussIntegrator1d(double a, double b) {
    _h = (b - a)/2;

    _points = new double[2];
    double c = (a + b)/2;
    _points[0] = c - _h*sqrt3;
    _points[1] = c + _h*sqrt3;
}

double GaussIntegrator1d::getFirstPoint() {
    return _points[0];
}

double GaussIntegrator1d::getSecondPoint() {
    return _points[1];
}

double GaussIntegrator1d::calculate(double point1, double point2) {
    return _h * (point1 + point2);
}