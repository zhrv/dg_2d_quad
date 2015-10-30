#include <iostream>
#include "cmath"
#include "global.h"
#include "GaussIntegrator1d.h"
#include "GaussIntegrator2D.h"

/**
 * du/dt + Vx*du/dx + Vy*du/dy=0
 */

using namespace std;

const int N = 200;

const double H = 1.0/N;

const double TAU = 0.1*H;
const double TMAX = 200*TAU;

const int FILE_SAVE_STEP = 1;
const int PRINT_STEP = 1;

VECTOR u[N+2][N+2];
Point  c[N+2][N+2]; // центры ячеек

VECTOR u1[N+2][N+2]; // сюда суммируем интегралы

// матрица масс
double** massMatrix;

int step = 0;
double t = 0.0;

VECTOR getF(int i, int j, Point pt);

inline double getU(int i, int j, Point pt) { return VECTOR::SCALAR_PROD(u[i][j], getF(i, j, pt)); }

void calculateMassMatrix();
void calculateDoubleIntegral();
void calculateLineIntegral();
void calculateSolution();
void incrementTime();
void output();
void saveResult(char *fName);
void fillWithInitialData();
void calculateCellCenters();
void set_bounds();

int main() {
     // заполняем вспомогательные массивы
    calculateCellCenters();
    calculateMassMatrix();
    fillWithInitialData();
    output();

    while (t < TMAX) {
        incrementTime();

        // Граничные условия
        set_bounds();

        // Правая часть уравнения
        calculateDoubleIntegral(); // вычисляем двойной интеграл
        calculateLineIntegral(); // вычисляем криволинейный интеграл по границе квадрата

        // Вычисляем значения решения в ячейках использую предыдущий слой и правую часть уравнения
        calculateSolution();

        // вывод через определенное количество шагов
        output();
    }
    // завершение приложения, очистка памяти..
    return 0;
}

void calculateMassMatrix() {
    double diag[3] = { 1.0/(H*H), 12.0/(H*H), 12/(H*H) };
    massMatrix = new double*[3];
    for (int i = 0; i < 3; ++i) {
        massMatrix[i] = new double[3];
        for (int j = 0; j < 3; ++j) {
            massMatrix[i][j] = i == j ? diag[i] : 0.0;
        }
    }
}

void calculateDoubleIntegral() {
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {



            GaussIntegrator2d gi( c[i][j].x - H/2.0, c[i][j].x + H/2.0,
                                  c[i][j].y - H/2.0, c[i][j].y + H/2.0 );
            int x = 0;
            int y = 1;

            double** points = gi.getPoints();
            double values[gi.N];

            for (int k = 0; k < gi.N; ++k) {
                Point point(points[k][x], points[k][y]);

                values[k] = getU(i, j, point);
            }

            double s = gi.calculate(values);
            //VECTOR v(3);
            //v = getDFDX()
            u1[i][j][0] = 0.0;
            u1[i][j][1] = s/H;
            u1[i][j][2] = s/H;
        }
    }
}

#define XDIR 0
#define YDIR 1

#define LSIDE -1
#define RSIDE 1

#define TSIDE 1
#define BSIDE -1

VECTOR calculateFlux(int i, int j,int direction, int side){

    double a;
    double b;

    if(direction == XDIR) {
        a = c[i][j].y-H/2.0;
        b = c[i][j].y+H/2.0;
    } else {
        a = c[i][j].x-H/2.0;
        b = c[i][j].x+H/2.0;
    }
    const double SQRT3 = 1.0/sqrt(3.0);
    double gp[2] = {0.5*(a+b)-0.5*(b-a)*SQRT3, 0.5*(a+b)+0.5*(b-a)*SQRT3};
    double gJ = 0.5*(b-a);

    VECTOR values[2];

    Point points[2];

    if(direction == XDIR) {
        points[0] = Point(c[i][j].x+side*H/2.0, gp[0]);
        points[1] = Point(c[i][j].x+side*H/2.0, gp[1]);
    } else {
        points[0] = Point(gp[0], c[i][j].y+side*H/2.0);
        points[1] = Point(gp[1], c[i][j].y+side*H/2.0);
    }

    int i1, j1;

    if(direction == XDIR) {
        i1 = (side == LSIDE) ? i-1 : i;
        j1 = j;
    } else {
        i1 = i;
        j1 = (side == BSIDE) ? j-1 : j;
    }

    for (int k = 0; k < 2; ++k) {
        values[k]  = getF(i,j,points[k]);
        values[k] *= getU(i1, j1, points[k]);
    }

    VECTOR flux(3);

    for (int l = 0; l < 3; ++l) {
        flux[l] = gJ*(values[0][l]+values[1][l]);
    }

    return flux;
}

void calculateLineIntegral() {
    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            u1[i][j] += calculateFlux(i, j, XDIR, LSIDE); // |
            u1[i][j] -= calculateFlux(i, j, XDIR, RSIDE); //   |

            u1[i][j] -= calculateFlux(i, j, YDIR, TSIDE); //  -
            u1[i][j] += calculateFlux(i, j, YDIR, BSIDE); //  _
        }
    }
}

void calculateSolution() {
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {

            u1[i][j] *= massMatrix;
            u1[i][j] *= TAU;

            u[i][j] += u1[i][j];
        }
    }
}

VECTOR getF(int i, int j, Point pt) {
    VECTOR v(3);
    v[0] = 1.0;
    v[1] = (pt.x - c[i][j].x)/H;
    v[2] = (pt.y - c[i][j].y)/H;
    return v;
}

VECTOR getDFDX(int i, int j, Point pt) {
    VECTOR v(3);
    v[0] = 0.0;
    v[1] = 1.0/H;
    v[2] = 0.0;
    return v;
}

VECTOR getDFDY(int i, int j, Point pt) {
    VECTOR v(3);
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 1.0/H;
    return v;
}

void output() {
    if (step % FILE_SAVE_STEP  == 0)
    {
        char fName[50];
        sprintf(fName, "res_%010d.csv", step);

        saveResult(fName);
    }

    if (step % PRINT_STEP == 0)
    {
        printf("step: %d\t\ttime step: %.16f\n", step, t);
    }
}

void saveResult(char *fName) {
    FILE * fp = fopen(fName, "w");

    fprintf(fp, "x,y,u\n");
    printf("File '%s' saved...\n", fName);

    for(int i = 1; i <= N; i++) {
        for(int j = 1; j <= N; j++) {
            double val = getU(i, j, c[i][j]);
            fprintf(fp, "%f,%f,%f\n", c[i][j].x, c[i][j].y, val);
        }
    }

    fclose(fp);
}

void incrementTime() {
    step++;
    t += TAU;
}

void fillWithInitialData() {
    const double r = 0.2;
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            u[i][j] = VECTOR(3);
            u1[i][j] = VECTOR(3);

            u[i][j] = 0;
            u1[i][j] = 0;

            Point lc = c[i][j];
            if (lc.x < 0.5 + r && lc.x > 0.5 - r && lc.y < 0.5 + r && lc.y > 0.5 - r)
            {
                u[i][j][0] = 10;
                u1[i][j][0] = 10;
            }
        }
    }
}

void calculateCellCenters() {
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            c[i][j] = Point(i*H + H/2, j*H + H/2);
        }
    }
}

void set_bounds() {
    for (int i = 1; i <= N; i++) {
        u[0][i] = u[N][i];
        u[N+1][i] = u[1][i];
        u[i][0] = u[i][N];
        u[i][N+1] = u[i][0];
    }
}


