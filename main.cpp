#include <iostream>
#include "cmath"
#include "global.h"
/**
 * du/dt + Vx*du/dx + Vy*du/dy=0
 */

using namespace std;

const int N = 100;

const double H = 1.0/N;

const double TMAX = 0.3;
const double TAU = 1.e-2;

const int FILE_SAVE_STEP = 10;
const int PRINT_STEP = 1;

VECTOR u[N+2][N+2];
Point  c[N+2][N+2]; // центры ячеек

VECTOR u1[N+2][N+2]; // сюда суммируем интегралы

Point cellGP[4];
double cellGW[4] = {1, 1, 1, 1};
const double J = 0.25*H*H; //якобиан

// матрица масс
double			****matrA;
double			****matrInvA;

int step =0;
double t = 0.0;

const int FUNC_COUNT = 3;

double sqrt3 = 1.0/3.0;

VECTOR getF(int i, int j, Point pt);
VECTOR getDFDX(int iCell, Point pt);
VECTOR getDFDY(int iCell, Point pt);

inline double getU(int i, int j, Point pt) { return VECTOR::SCALAR_PROD(u[i][j], getF(i, j, pt)); }

void fillGaussPoints();
void calculateMassMatrix();
void calculateDoubleIntegral();
void calculateLineIntegral();
void calculateSolution();
void incrementTime();
void output();
void saveResult(char *fName);
void fillWithInitialData();
void calculateCellCenters();

void transposedMatr( double **a, int n) {
	for(int row = 0; row < n - 1; row++)
        for(int col = row + 1; col < n; col++)
            std::swap(a[col][row], a[row][col]);
    }

double determinantMatr(double **a, int N) { 
	int i,j; 
	double **matr1; 
	double determ=0; 

		if (N==2) 
			{ 
				determ=a[0][0]*a[1][1]-a[0][1]*a[1][0]; 
			} 
		else 
			{ 
				matr1=new double*[N-1]; 

				for(i=0;i<N;i++) { 
					for(j=0;j<N-1;j++) { 
						if(j<i) { matr1[j]=a[j]; } 
						else { matr1[j]=a[j+1]; } 
					} 
					determ+=pow(-1,(i+j))*determinantMatr(matr1,N-1)*a[i][N-1]; 
			} 
			delete matr1; 
			} 

		return determ; 
} 

double** multMatrix(double **a, double **b,int n){
	double **result = new double*[n];
	for (int i=0;i<n;i++)
	{
		result[i] = new double[n];
	}

	for (int i=0;i<n;i++)
	{
		for (int j=0;j<n;j++)
		{
			result[i][j] = 0;
			for (int k=0;k<n;k++)
			{
				result[i][j] += a[i][k]*b[k][j];
				cout<<result[i][j]<<"  ";
			}
			cout<<endl;
		}
	}

	return result;
}

void multMatrixNum(double** a, double num, int n){
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<n;j++)
		{
			a[i][j] *= num;
		}
	}

}
double* multMatrixVect(double **a, double *b, int n){
	double *result = new double[n];

	for (int i=0;i<n;i++)
	{
		result[i] = 0;
		for (int j=0;j<n;j++)
		{
			result[i] += a[i][j]*b[j];
		}
	}

	return result;
}

int main() {
    // заполняем вспомогательные массивы
    calculateCellCenters();
    fillGaussPoints();
    calculateMassMatrix();
    fillWithInitialData();

    while (t < TMAX) {
        incrementTime();

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

void fillGaussPoints()
{
    cellGP[0].x = -sqrt3;
    cellGP[0].y = -sqrt3;

    cellGP[1].x =  sqrt3;
    cellGP[1].y = -sqrt3;

    cellGP[2].x = sqrt3;
    cellGP[2].y = sqrt3;

    cellGP[3].x = -sqrt3;
    cellGP[3].y =  sqrt3;
    // todo: cellGW ??????
}

void calculateMassMatrix()
{
    matrA = new double***[N+2];
    matrInvA = new double***[N+2];
    for (int i = 0; i < N+2; i++) {
        matrA[i] = new double**[N+2];
        matrInvA[i] = new double**[N+2];
        for (int j = 0; j < N+2; j++) {
            matrA[i][j] = new double*[FUNC_COUNT];
            matrInvA[i][j] = new double*[FUNC_COUNT];
            for (int k = 0; k < FUNC_COUNT; k++) {
                matrA[i][j][k] = new double[FUNC_COUNT];
                matrInvA[i][j][k] = new double[FUNC_COUNT];
            }
        }
    }

    // вычисляем матрицу масс
    for (int i = 1; i <= N; i++)
	{
				for (int j = 1; j <= N; j++) {
					double **a = matrA[i][j];

				for (int m = 0; m < FUNC_COUNT; m++) 
					for (int n = 0; n < FUNC_COUNT; n++) 
						a[m][n] = 0.0;
						
						for (int igp = 0; igp < 4; igp++)
						{
								Point gp = cellGP[i];
								gp *= H; // todo: неоходимо изменить если не квадратные ячейки
								gp *= 0.5;
								gp += c[i][j];
								VECTOR phi1 = getF(i,j, gp);
								VECTOR phi2 = getF(i,j, gp);
								for (int m = 0; m < FUNC_COUNT; m++) {
									for (int n = 0; n < FUNC_COUNT; n++) {
										a[m][n] += H*H* phi1[m] * phi2[n] * cellGW[igp] * 0.25;
									}
								}
						}
						//обратная матрица
						double detA=determinantMatr(a,FUNC_COUNT);
						transposedMatr(a,FUNC_COUNT);
						multMatrixNum(a,detA,FUNC_COUNT);

						for (int m = 0; m < FUNC_COUNT; m++) {
									for (int n = 0; n < FUNC_COUNT; n++) {
										matrInvA[i][j][m][n]=a[m][n];
									}
								}


				}
    }
}

void calculateDoubleIntegral()
{
    for (int i = 1; i < N+1; i++) {
        for (int j = 0; j < N+1; j++) {
            VECTOR s(FUNC_COUNT);
            s = 0.0;
            for (int igp = 0; igp < 4; igp++) {
                Point gp = cellGP[i];
                gp *= H; // todo: неоходимо изменить если не квадратные ячейки
                gp *= 0.5;
				gp += c[i][j];
                s += (getU(i,j,gp)*(getDFDX(i,gp)+getDFDY(i,gp)))*H*H*0.25;
                s *= cellGW[igp];
            }

            u1[i][j] = s;
        }
    }
}

void lineIntegralXDirection()
{
    for (int i = 0; i < N+1; i++) {
        for (int j = 1; j < N+1; j++) {
            Point gp1, gp2;
            gp1 = c[i][j];
            gp1.x += 0.5*H;
            gp2 = gp1;
            gp1.y -= sqrt3*H*0.5;
            gp2.y += sqrt3*H*0.5;
            double flux1 = getU(i, j, gp1);
            double flux2 = getU(i, j, gp2);

            VECTOR sL(FUNC_COUNT), sR(FUNC_COUNT), tmp(FUNC_COUNT);
            tmp = getF(i,j,gp1);
            tmp *= flux1;
            sL += tmp;
            tmp = getF(i,j,gp1);
            tmp *= flux2;
            sL += tmp;
            sL *= H*0.5;

            tmp = getF(i,j,gp2);
            tmp *= flux1;
            sR += tmp;
            tmp = getF(i,j,gp2);
            tmp *= flux2;
            sR += tmp;
            sR *= H*0.5;

            u1[i][j]   -= sL;
            u1[i+1][j] += sR;
        }
    }
}

void lineIntegralYDirection()
{
    for (int i = 1; i < N+1; i++) {
        for (int j = 0; j < N+1; j++) {
            Point gp1, gp2;
            gp1 = c[i][j];
            gp1.y += 0.5*H;
            gp2 = gp1;
            gp1.x -= sqrt3*H*0.5;
            gp2.x += sqrt3*H*0.5;
            double flux1 = getU(i, j, gp1);
            double flux2 = getU(i, j, gp2);

            VECTOR sL(FUNC_COUNT), sR(FUNC_COUNT), tmp(FUNC_COUNT);
            tmp = getF(i,j,gp1);
            tmp *= flux1;
            sL += tmp;
            tmp = getF(i,j,gp1);
            tmp *= flux2;
            sL += tmp;
            sL *= H*0.5;

            tmp = getF(i,j,gp2);
            tmp *= flux1;
            sR += tmp;
            tmp = getF(i,j,gp2);
            tmp *= flux2;
            sR += tmp;
            sR *= H*0.5;

            u1[i][j]   -= sL;
            u1[i][j+1] += sR;
        }
    }
}

void calculateLineIntegral()
{
    lineIntegralXDirection(); // x-direction

    lineIntegralYDirection(); // y-direction
}

void calculateSolution()
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            u1[i][j] *= matrInvA[i][j];
            u1[i][j] *= TAU;

            u[i][j] += u1[i][j];
        }
    }
}

VECTOR getF(int i, int j, Point pt)
{
    if (FUNC_COUNT == 3 )
    {
        VECTOR v(3);
        v[0] = 1.0;
        v[1] = (pt.x - c[i][j].x)/H;
        v[2] = (pt.y - c[i][j].y)/H;
        return v;
    } else {
        std::cout << "ERROR: wrong basic functions count!\n";
        return 1;
    }
}

VECTOR getDFDX(int iCell, Point pt)
{
    if (FUNC_COUNT == 3 )
    {
        VECTOR v(3);
        v[0] = 0.0;
        v[1] = 1.0/H;
        v[2] = 0.0;
        return v;
    } else {
        std::cout << "ERROR: wrong basic functions count!\n";
        return 1;
    }
}

VECTOR getDFDY(int iCell, Point pt)
{
    if (FUNC_COUNT == 3 )
    {
        VECTOR v(3);
        v[0] = 0.0;
        v[1] = 0.0;
        v[2] = 1.0/H;
        return v;
    } else {
        std::cout << "ERROR: wrong basic functions count!\n";
        return 1;
    }
}

void output()
{
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

void saveResult(char *fName)
{
    FILE * fp = fopen(fName, "w");

    fprintf(fp, "x,y,u\n");
    printf("File '%s' saved...\n", fName);

   for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
           // double val = getU(i, j, с[i][j]);
            fprintf(fp, "%f,%f,%f\n", i*H, j*H, u[i][j]);
        }
    }
    fclose(fp);
}

void incrementTime()
{
    step++;
    t += TAU;
}

void fillWithInitialData()
{
    const double r = 0.2;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            //u[i][j] = VECTOR(3);
            Point lc = c[i][j];
            if ((lc.x-0.5)*(lc.x-0.5) + (lc.y-0.5)*(lc.y-0.5) < r*r)
            {
                u[i][j]= 10;
            }
        }
    }
}

void calculateCellCenters() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            c[i][j] = Point(i*H + H/2, j*H + H/2);
        }
    }
}