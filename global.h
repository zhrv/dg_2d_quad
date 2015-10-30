#ifndef _DG_2D_QUAD_GLOBAL_H_
#define _DG_2D_QUAD_GLOBAL_H_

#include <cstring>

struct Point {
    double x;
    double y;

    Point() : x(0.0), y(0.0) { }

    Point(double ax, double ay) : x(ax), y(ay) { }

    inline void operator*=(double q) {
        x *= q;
        y *= q;
    }

    inline void operator/=(double q) {
        x /= q;
        y /= q;
    }

    inline void operator=(double q) {
        x = q;
        y = q;
    }

    inline void operator=(Point p) {
        x = p.x;
        y = p.y;
    }

    inline void operator+=(Point p) {
        x += p.x;
        y += p.y;
    }

    inline void operator-=(Point p) {
        x -= p.x;
        y -= p.y;
    }

    inline void operator+=(double q) {
        x += q;
        y += q;
    }

    inline void operator-=(double q) {
        x -= q;
        y -= q;
    }

    inline double length() { return sqrt(x * x + y * y); }

    friend inline double operator*(Point p1, Point p2) { return p1.x * p2.x + p1.y * p2.y; }

    friend inline Point operator+(Point p1, Point p2) { return Point(p1.x + p2.x, p1.y + p2.y); }

    friend inline Point operator-(Point p1, Point p2) { return Point(p1.x - p2.x, p1.y - p2.y); }
};


typedef Point Vector;

struct VECTOR {
    VECTOR(int an = 0) : n(an) {
        if (n) {
            elem = new double[n];
            memset(elem, 0, sizeof(double) * n);
        } else elem = NULL;
    }

    VECTOR(const VECTOR &v) : n(v.n) {
        elem = new double[n];
        memcpy(elem, v.elem, n * sizeof(double));
    }

    ~VECTOR() {
        if (elem) delete[] elem;
        n = 0;
    }

    void init(int an) {
        n = an;
        if (elem) delete[] elem;
        elem = new double[n];
    }

    double &operator[](int i) { return elem[i]; }

    double &operator()(int i) { return elem[i]; }

    void operator=(const VECTOR &v) {
        if (elem) delete[] elem;
        n = v.n;
        elem = new double[n];
        memcpy(elem, v.elem, n * sizeof(double));
    }

    void operator=(const double &x) { for (int i = 0; i < n; i++) elem[i] = x; }

    void operator+=(const VECTOR &v) { for (int i = 0; i < n; i++) elem[i] += v.elem[i]; }

    void operator+=(const double &x) { for (int i = 0; i < n; i++) elem[i] += x; }

    void operator-=(const VECTOR &v) { for (int i = 0; i < n; i++) elem[i] -= v.elem[i]; }

    void operator-=(const double &x) { for (int i = 0; i < n; i++) elem[i] -= x; }

    void operator*=(const double &x) { for (int i = 0; i < n; i++) elem[i] *= x; }

    void operator*=(double **matr) {
        double *tmp = new double[n];
        for (int i = 0; i < n; i++) {
            tmp[i] = 0.0;
            for (int j = 0; j < n; j++) {
                tmp[i] += matr[i][j] * elem[j];
            }
        }
        delete[] elem;
        elem = tmp;
    }

    void zero() {
        if (elem != NULL && n != 0) {
            memset(elem, 0, n * sizeof(double));
        }
    }

    void abs() {
        for (int i = 0; i < n; i++) {
            elem[i] = fabs(elem[i]);
        }
    }

    double norm() {
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            s += elem[i] * elem[i];
        }
        return sqrt(s);
    }

    static double SCALAR_PROD(const VECTOR &v1, const VECTOR &v2) {
        double s = 0;
        for (int i = 0; i < v1.n; i++) {
            s += v1.elem[i] * v2.elem[i];
        }
        return s;
    }

    double *elem;
    int n;
};

const VECTOR operator+(const VECTOR &left, const VECTOR &right) {
    VECTOR result(left);
    result += right;
    return result;
}

const VECTOR operator-(const VECTOR &left, const VECTOR &right) {
    VECTOR result(left);
    result -= right;
    return result;
}

const VECTOR operator*(const double &left, const VECTOR &right) {
    VECTOR result(right);
    result *= left;
    return result;
}

const VECTOR operator*(const VECTOR &left, const double &right) {
    return right * left;
}

#endif