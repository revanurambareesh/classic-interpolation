#include <stdio.h>
#include "..\Resources\Graphics\graphics.h"

#define MAXSET 10
#define AXISLENGTH 300

//Initialization
void initInterpolation();

//CUBIC-SPLINE INTERPOLATION
double implement_CubicSplineandBackFit(double x[], double y[], int num, double pointY, double accuracy);
void cubicSplineCoeff(double x[], double y[], int num);
double poly_newtonRaphson(int power, double coeff[], float accuracy, double initApprox);
double poly_fvalue(int power, double coeff[], double x);
double poly_f1dash(int power, double coeff[], double x);
double curvefitCSI(int num_of_points, double pointX);
double backfitCSI(int num_of_points, double pointY);
void calculateCSIinterpolation(int num_of_points);

//NEVILLE'S POLYNOMIAL INTERPOLATION
double neville_poly(int num_of_points, double y[], double x[], double pointXY);
double curvefitNPoly(int num_of_points, double pointX);
double backfitNPoly(int num_of_points, double pointY);
void calculateNpolyInterpolation(int num_of_points);

//LINEAR FIT
double curvefitLinear(int n, double pointX);
double backfitLinear(int n, double pointY);
int guassJordan(double equation[][6], double solution[], int m, int n);
void calculateBestFitPow1(int num_of_points);
//double poly_newtonRaphson(int power, double coeff[], float accuracy, double initApprox);
//double poly_fvalue(int power, double coeff[], double x);
int initCoeffLinear(double x[], double y[], double LinearABsolution[], int n, double PolyCoeff1[][6]);
//double poly_f1dash(int power, double coeff[], double x);

//POLY2 FIT
void calculateBestFitPow2(int num_of_points);
double curvefitPoly2(int n, double pointX);
double backfitPoly2(int n, double pointY);
//int guassJordan(double equation[][5], double solution[], int m, int n);
//double poly_newtonRaphson(int power, double coeff[], float accuracy, double initApprox);
//double poly_fvalue(int power, double coeff[], double x);
int initCoeffPoly2(double x[], double y[], double PolyABCsolution[], int n, double PolyCoeff2[][6]);
//double poly_f1dash(int power, double coeff[], double x);

//POLY3 FIT
void calculateBestFitPow3(int num_of_points);
double curvefitPoly3(int n, double pointX);
double backfitPoly3(int n, double pointY);
//int guassJordan(double equation[][5], double solution[], int m, int n);
//double poly_newtonRaphson(int power, double coeff[], float accuracy, double initApprox);
//double poly_fvalue(int power, double coeff[], double x);
int initCoeffPoly3(double x[], double y[], double PolyABCDsolution[], int n, double PolyCoeff3[][6]);
//double poly_f1dash(int power, double coeff[], double x);

//POLY4 FIT
void calculateBestFitPow4(int num_of_points);
double curvefitPoly4(int n, double pointX);
double backfitPoly4(int n, double pointY);
//int guassJordan(double equation[][5], double solution[], int m, int n);
//double poly_newtonRaphson(int power, double coeff[], float accuracy, double initApprox);
//double poly_fvalue(int power, double coeff[], double x);
int initCoeffPoly4(double x[], double y[], double PolyABCDEsolution[], int n, double PolyCoeff4[][6]);
//double poly_f1dash(int power, double coeff[], double x);


void plotCurve(int n);

//---------------------------------------//
//VARIABLES
//INPUTS (undisturbed)
double oriX[MAXSET] = { 0, 1, 4, 9, 16},oriY[MAXSET] = { 0, 1, 2, 3, 4 };
int num_of_points = 5;
double pointY = 4, pointX = 1;
//COEFFICIENTS
double CSIcoeff[MAXSET * 4 - 4]; //cubicSpline
double PolyCoeff1[2][6];//Best fit -Linear
double PolyCoeff2[3][6];//Best fit polynomial2
double PolyCoeff3[4][6];//Best fit polynomial3
double PolyCoeff4[5][6];//Best fit polynomial4

//Y-values (interpolated)
double CSIinterpolation[MAXSET];
double NpolyInterpolation[MAXSET];
double linearFit[MAXSET];
double Poly2ndOrder[MAXSET];
double Poly3rdOrder[MAXSET];
double Poly4thOrder[MAXSET];

//Others essential variables
double Bo; //Max value on Y-axis
double x[MAXSET], y[MAXSET]; //Used variables
double LinearABsolution[2], PolyABCsolution[3], PolyABCDsolution[4], PolyABCDEsolution[5]; //These are the polynomial coefficients

int main()
{
	int gd = DETECT, gm, i;
	double backfitX, curvefitY;

	//STEP 1: initialize the values
	initInterpolation();

	//STEP 2: calculate coeff
	cubicSplineCoeff(x, y, num_of_points);
	if(initCoeffLinear(x, y, LinearABsolution, num_of_points, PolyCoeff1) == -9994)printf("Cannot Implement Linear\n");
	if(initCoeffPoly2(x, y, PolyABCsolution, num_of_points, PolyCoeff2) == -9994)printf("Cannot Implement Poly2\n");
	if(initCoeffPoly3(x, y, PolyABCDsolution, num_of_points, PolyCoeff3) == -9994)printf("Cannot Implement Poly3\n");
	if(initCoeffPoly4(x, y, PolyABCDEsolution, num_of_points, PolyCoeff4) == -9994)printf("Cannot Implement Poly4\n");

	//Step 3: Calculate interpolated values.
	calculateCSIinterpolation(num_of_points);
	calculateNpolyInterpolation(num_of_points);
	calculateBestFitPow1(num_of_points);
	calculateBestFitPow2(num_of_points);
	calculateBestFitPow3(num_of_points);
	calculateBestFitPow4(num_of_points);


	//Step 4: Backfitting
	//Cubic Spline Interpolation
	backfitX = backfitCSI(num_of_points, pointY);
	//Curve fitting
	curvefitY = curvefitCSI(num_of_points, pointX);

	printf("Backfit (Cubic Spline Interpolation): %.3lf\n", backfitX);
	printf("Curvefit (Cubic Spline Interpolation): %.3lf\n", curvefitY);

	//Neville's polynomial
	//backfitting
	backfitX = backfitNPoly(num_of_points, pointY);
	//Curve fitting
	curvefitY = curvefitNPoly(num_of_points, pointX);

	printf("Backfit (Neville's Polynomial Interpolation): %.3lf\n", backfitX);
	printf("Curvefit (Neville's Polynomial Interpolation): %.3lf\n", curvefitY);

	//LinearCurveFit
	backfitX = backfitLinear(num_of_points, pointY);
	//Curve fitting
	curvefitY = curvefitLinear(num_of_points, pointX);

	printf("Backfit (Linear Best Fit): %.3lf\n", backfitX);
	printf("Curvefit (Linear Best Fit): %.3lf\n", curvefitY);

	//Polynomial 2nd order
	backfitX = backfitPoly2(num_of_points, pointY);
	//Curve fitting
	curvefitY = curvefitPoly2(num_of_points, pointX);

	printf("Backfit (Polynomial 2nd order): %.3lf\n", backfitX);
	printf("Curvefit (Polynomial 2nd order): %.3lf\n", curvefitY);

	//Polynomial 3rd order
	backfitX = backfitPoly3(num_of_points, pointY);
	//Curve fitting
	curvefitY = curvefitPoly3(num_of_points, pointX);

	printf("Backfit (Polynomial 3rd order): %.3lf\n", backfitX);
	printf("Curvefit (Polynomial 3rd order): %.3lf\n", curvefitY);

	//Polynomial 4th order
	backfitX = backfitPoly4(num_of_points, pointY);
	//Curve fitting
	curvefitY = curvefitPoly4(num_of_points, pointX);

	printf("Backfit (Polynomial 4th order): %.3lf\n", backfitX);
	printf("Curvefit (Polynomial 4th order): %.3lf\n", curvefitY);


	//Error observation table
	printf("\n\nError Observation Table\n");
	printf("X\tY\tCubicY\tNevileY\tLinearY\tPoly2Y\tPoly3Y\tPoly4Y\n");

	for (i = 0; i < num_of_points; i++)
	{
		printf("%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.1lf\n", oriX[i], oriY[i],
         CSIinterpolation[i],
         NpolyInterpolation[i],
         linearFit[i],
         Poly2ndOrder[i],
         Poly3rdOrder[i],
         Poly4thOrder[i]);
	}

    //Backfitting accuracy
    printf("\n\nBackfitting accuracy Table\n");
    printf("X\tY\tCubicX\tNevileX\tLinearX\tPoly2X\tPoly3X\tPoly4X\n");

    for (i = 0; i < num_of_points; i++)
	{
		printf("%.3lf\t%.1lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n", oriX[i], oriY[i],
         backfitCSI(num_of_points, oriY[i]),
         backfitNPoly(num_of_points, oriY[i]),
         backfitLinear(num_of_points, oriY[i]),
         backfitPoly2(num_of_points, oriY[i]),
         backfitPoly3(num_of_points, oriY[i]),
         backfitPoly4(num_of_points, oriY[i]));
	}

	//CURVE
	initgraph(&gd, &gm, "C:\\TC\\BGI");
	plotCurve(num_of_points);
	closegraph();
	return 0;
}

void plotCurve(int n)
{

	double maxX = x[0], maxY = y[0];
	int i;
	double j;
	double tempx[MAXSET];
	double tempy[MAXSET];
	int offset = 70, netXaxis = AXISLENGTH + offset, netYaxis = AXISLENGTH + offset;
	int beta = AXISLENGTH / 100;
	double px, py;

	for (i = 1; i < n; i++)
	{
		if (maxX < x[i])maxX = x[i];
		if (maxY < y[i])maxY = y[i];
	}

	//XY axis
	for (i = offset; i < netYaxis; i = i + 1)putpixel(offset, i, BROWN);
	for (i = offset; i < netXaxis; i = i + 1)putpixel(i, offset + AXISLENGTH, BROWN);

	//Grid
	for (i = offset; i < netYaxis; i = i + 30)
	for (j = offset; j <= netXaxis; j = j + 30)putpixel(j, i, LIGHTGRAY);
	for (i = offset; i < netXaxis; i = i + 30)
	for (j = 0; j <= AXISLENGTH; j = j + 30)putpixel(i, offset + AXISLENGTH - j, LIGHTGRAY);

	for (i = 0; i < n; i++)
	{
		tempx[i] = (x[i] / maxX)*100.0;
		tempy[i] = (y[i] / maxY)*100.0;
	}

	//PLOTTING data
	for (i = 0; i <= netXaxis - offset; i++)
	{
		px = i*(double)(beta);

		py = curvefitNPoly(n, ((i) / beta)*(maxX / 100.0));
		py = (py / Bo)*100.0*beta;
		putpixel(offset + i, netYaxis - py, CYAN);//NEVILLE

		py = curvefitCSI(n, ((i) / beta)*(maxX / 100.0));
		py = (py / Bo)*100.0*beta;
		putpixel(offset + i, netYaxis - py, RED); //CubicSpline

		py = curvefitLinear(n, ((i) / beta)*(maxX / 100.0));
		py = (py / Bo)*100.0*beta;
		putpixel(offset + i, netYaxis - py, BLUE); //Linear Fit

		py = curvefitPoly2(n, ((i) / beta)*(maxX / 100.0));
		py = (py / Bo)*100.0*beta;
		putpixel(offset + i, netYaxis - py, DARKGRAY); //Poly2 Fit

		py = curvefitPoly3(n, ((i) / beta)*(maxX / 100.0));
		py = (py / Bo)*100.0*beta;
		putpixel(offset + i, netYaxis - py, WHITE); //Poly3 Fit

		py = curvefitPoly4(n, ((i) / beta)*(maxX / 100.0));
		py = (py / Bo)*100.0*beta;
		putpixel(offset + i, netYaxis - py, YELLOW); //Poly4 Fit
	}

	//Marking xy points
	for (i = 0; i < n; i += 1)
	{
		px = tempx[i] * ((double)(netXaxis - offset) / 100);
		py = tempy[i] * (double)beta;
		putpixel(offset + px, netYaxis - py, WHITE);
		putpixel(offset + px - 1, netYaxis - py - 1, WHITE);
		putpixel(offset + px + 1, netYaxis - py + 1, WHITE);
		putpixel(offset + px - 1, netYaxis - py + 1, WHITE);
		putpixel(offset + px + 1, netYaxis - py - 1, WHITE);
	}

	getch();
}

double curvefitCSI(int num_of_points, double pointX)
{
	double *coeff, max;
	double curvefitY, A, B, C, D, valX;
	int i, splineNum, flag;

	max = Bo;
	coeff = CSIcoeff;
	for (i = 0; i < num_of_points - 1; i++)
	{
		if ((pointX >= x[i] && pointX <= x[i + 1]) || (pointX <= x[i] && pointX >= x[i + 1]))
		{
			flag = 1;
			break;
		}
	}
	if (flag == 0)
	{
		//ERROR 002: y is out of bounds. Cannot Extrapolate with Interpolation formula
		//printf("Error csmain002");
		return -9997;
	}
	splineNum = i;

	A = *(coeff + 0 + splineNum * 4);
	B = *(coeff + 1 + splineNum * 4);
	C = *(coeff + 2 + splineNum * 4);
	D = *(coeff + 3 + splineNum * 4);

	//Si(x) = Ai + Bi(X - Xi) + Ci(X - Xi) ^ 2 + Di(X - Xi) ^ 3;

	valX = pointX - x[splineNum];
	curvefitY = A + B*valX + C*valX*valX + D*valX*valX*valX;

	curvefitY = (curvefitY / 100.0)*max;
	return curvefitY;
}

double curvefitNPoly(int num_of_points, double pointX)
{
	double *coeff, max;
	double curvefitY, A, B, C, D, valX;
	int i, splineNum, flag;

	max = Bo;
	curvefitY = neville_poly(num_of_points, y, x, pointX);
	curvefitY = (curvefitY / 100.0)*max;
	return curvefitY;
}

double backfitCSI(int num_of_points, double pointY)
{
	double backfitX;
	double max;
	int temp;

	max = Bo;
	//Backfitting
	pointY = (pointY / max)*100.0;
	temp = pointY * 100;
	if (temp % 10 >= 5)
	{
		temp -= temp % 10;
		temp += 10;
		pointY = temp / 100.0;
	}
	else
	{
		temp -= temp % 10;
		pointY = temp / 100.0;
	}

	backfitX = implement_CubicSplineandBackFit(x, y, num_of_points, pointY, 0.001);
	return backfitX;
}

double backfitNPoly(int num_of_points, double pointY)
{
	double backfitX;
	double max;
	int temp;

	max = Bo;
	//Backfitting
	pointY = (pointY / max)*100.0;
	temp = pointY * 100;
	if (temp % 10 >= 5)
	{
		temp -= temp % 10;
		temp += 10;
		pointY = temp / 100.0;
	}
	else
	{
		temp -= temp % 10;
		pointY = temp / 100.0;
	}

	backfitX = neville_poly(num_of_points, x, y, pointY);
	return backfitX;
}

void initInterpolation()
{
	int i, temp;
	double max;

	for (i = 0; i < num_of_points; i++)
	{
		if (i == 0)max = oriY[0];
		else if (max < oriY[i])max = oriY[i];
	}

	for (i = 0; i < num_of_points; i++)
	{
		x[i] = oriX[i];
		y[i] = (oriY[i] / max)*100.0;
		temp = y[i] * 100;
		if (temp % 10 >= 5)
		{
			temp -= temp % 10;
			temp += 10;
			y[i] = temp / 100.0;
		}
		else
		{
			temp -= temp % 10;
			y[i] = temp / 100.0;
		}
	}

	Bo = max;
}

double implement_CubicSplineandBackFit(double x[], double y[], int num, double pointY, double accuracy)
{
	double backfit_coeff[MAXSET];
	int i, splineNum, flag = 0;
	double *a;
	double backfitX, initGuess;

	//Complete Implementation of Cubic Spline interpolation with Backfit methods.
	a = CSIcoeff;
	//num is number of points
	if (a == NULL)return -9996;

	//BackFitting

	for (i = 0; i < num - 1; i++)
	{
		if ((pointY >= y[i] && pointY <= y[i + 1]) || (pointY <= y[i] && pointY >= y[i + 1]))
		{
			flag = 1;
			break;
		}
	}

	if (flag == 0)
	{
		//ERROR 002: y is out of bounds. Cannot Extrapolate with Interpolation formula
		//printf("Error csmain002");
		return -9997;
	}

	splineNum = i;
	initGuess = (x[splineNum + 1] - x[splineNum]) / 2;//To be improved... using Newton Raphson's method

	for (i = 0; i < 4; i++)
	{
		backfit_coeff[3 - i] = *(a + i + splineNum * 4);
	}

	backfit_coeff[3] -= pointY;
	backfitX = poly_newtonRaphson(3, backfit_coeff, accuracy, initGuess);

	if (backfitX == -9995)return -9995;
	return backfitX + x[splineNum];
}

void cubicSplineCoeff(double x[], double y[], int num)
{
	double h[MAXSET], b[MAXSET], z[MAXSET], u[MAXSET], v[MAXSET];
	double temp;
	int i, ix;
	int n = num - 1;

	/*
	NOTE:
	u[0] is not required.
	v[0] is not required.

	The Spline expression obtained will be of the form:
	Si(x) = Ai + Bi(X-Xi)+ Ci(X-Xi)^2 + Di(X-Xi)^3;
	CSIcoeff stores A B C D in same order
	*/

	if (num < 4 || num > MAXSET)
	{
		//ERROR 001: Less than four points
		//printf("ERROR csi001");
		return;
	}

	//To find Zi
	for (i = 0; i <= n - 1; i++)
	{
		h[i] = x[i + 1] - x[i];
		b[i] = (y[i + 1] - y[i]) / h[i];
	}

	//Gaussian elimination
	u[1] = 2 * (h[0] + h[1]);
	v[1] = 6 * (b[1] - b[0]);

	for (i = 2; i <= n - 1; i++)
	{
		u[i] = 2 * (h[i - 1] + h[i]) - (h[i - 1] * h[i - 1]) / u[i - 1];
		v[i] = 6 * (b[i] - b[i - 1]) - (h[i - 1] * v[i - 1]) / u[i - 1];
	}

	//Back-substitution
	z[n] = 0;
	for (i = n - 1; i >= 1; i--) z[i] = (v[i] - h[i] * z[i + 1]) / u[i];
	z[0] = 0;

	//abcd values
	for (i = 0; i < n * 4; i = i + 4)
	{
		ix = i / 4;
		CSIcoeff[i] = y[ix];
		CSIcoeff[i + 1] = -(h[ix] / 6.0)*z[ix + 1] - (h[ix] / 3.0)*z[ix] + (y[ix + 1] - y[ix]) / h[ix];
		CSIcoeff[i + 2] = z[ix] / 2.0;
		CSIcoeff[i + 3] = (z[ix + 1] - z[ix]) / (6.0 * h[ix]);
	}

	return;
}

double poly_newtonRaphson(int power, double coeff[], float accuracy, double initApprox)
{
	double xn = initApprox, xp;
	int count = 0;
	/*
	NOTE
	Function f(x) = ax^n + bx^n-1 + ... + c
	n = power;
	coeff[] = {a, b, ... , c}
	example for accuracy: 0.0001
	Number of terms in coeff[] must be power + 1;
	*/

	//Step 1: Approximating to accuracy
	xn = (double)((int)(xn / accuracy) * (double)accuracy);

	//Step 2: Newton-Raphson method to find positive root
	do
	{
		xp = xn;
		xn = xp - (poly_fvalue(power, coeff, xp) / poly_f1dash(power, coeff, xp));
		xn = (double)((int)(xn / accuracy) * (double)accuracy);
		if (xn - xp <= accuracy)break;
		count++;
	} while (xn != xp&&count < 999);
	if (count == 999)return -9995;
	return xn;
}

double poly_fvalue(int power, double coeff[], double x)
{
	double ans = 0, higherpow;
	int i, j;

	/*
	NOTE:
	Function f(x) = ax^n + bx^n-1 + ... + c
	n = power;
	coeff[] = {a, b, ... , c}
	Number of terms in coeff[] must be power + 1;
	*/

	for (i = 0; i <= power; i++)
	{
		higherpow = 1;
		j = power - i;
		if (coeff[i] != 0)while (j--)higherpow *= x;
		ans += coeff[i] * higherpow;
	}

	return ans;
}

double poly_f1dash(int power, double coeff[], double x)
{
	double ans = 0, higherpow;
	int i, j;

	/*
	NOTE:
	Function f(x) = ax^n + bx^n-1 + ... + c
	f'(x) = nax^n-1 + (n-1)bx^n-2 + ... + 0
	n = power;
	coeff[] = {a, b, ... , c}
	Number of terms in coeff[] must be power + 1;
	*/

	for (i = 0; i < power; i++)
	{
		higherpow = 1;
		j = power - i - 1;
		if (coeff[i] != 0)while (j--)higherpow *= x;
		ans += (power - i)*(higherpow * coeff[i]);
	}

	return ans;
}

double neville_poly(int n, double y[], double x[], double eval)
{
	/*
	NOTE:-
	For backfitting, pass X_values array for y[] and Y_values for x[].
	eval is pointY
	For generating graph, pass X_value array for x[] and Y_values for y[].
	eval will be the point on X axis.
	*/
	int i, j;
	double temp[MAXSET];
	double flag = 0;


	for (i = 0; i < n; i++)temp[i] = y[i];
	for (j = 1; j < n; j++)
	for (i = n - 1; i >= j; i--) {
		temp[i] = ((eval - x[i - j])*temp[i] - (eval - x[i])*temp[i - 1]) / (x[i] - x[i - j]);
	}

	return (temp[n - 1]);
}

void calculateCSIinterpolation(int num_of_points)
{
	int i;
	for (i = 0; i < num_of_points; i++)
	{
		CSIinterpolation[i] = curvefitCSI(num_of_points, oriX[i]);
	}
}

void calculateNpolyInterpolation(int num_of_points)
{
	int i;
	for (i = 0; i < num_of_points; i++)
	{
		NpolyInterpolation[i] = curvefitNPoly(num_of_points, oriX[i]);
	}
}

double curvefitLinear(int n, double pointX)
{
	double max = Bo;
	double curvefitY;
	curvefitY = LinearABsolution[0] * pointX + LinearABsolution[1];
	return (curvefitY / 100.0)*max;
}

double backfitLinear(int n, double pointY)
{
	double backfitX;
	double max;
	int temp;

	max = Bo;
	//Backfitting
	pointY = (pointY / max)*100.0;
	temp = pointY * 100;
	if (temp % 10 >= 5)
	{
		temp -= temp % 10;
		temp += 10;
		pointY = temp / 100.0;
	}
	else
	{
		temp -= temp % 10;
		pointY = temp / 100.0;
	}

	if (LinearABsolution[0] == 0)return -9993;
	return (pointY - LinearABsolution[1]) / LinearABsolution[0];
}

int guassJordan(double equation[][6], double solution[], int m, int n)
{
	/*
	This method is used to solve consistent and diagonally dominant matrix equation
	This method assumes that matrix does not have a zero element in its principle diagonal.
	n = m + 1
	*/

	int i, j, k, l;
	double ele;

	//Guass-elimination method.
	for (i = 0; i < m; i++)
	{
		for (j = i; j < m; j++)
		{
			ele = equation[j][i];

			if (ele == 0)return -9994;
			for (k = 0; k < n; k++)
			{
				equation[j][k] /= ele;
			}
		}

		for (j = i + 1; j < m; j++)
		{
			for (k = 0; k < n; k++)
				equation[j][k] -= equation[i][k];
		}
	}


	for (i = m - 1; i >= 0; i--)
	{
		for (j = i; j >= 0; j--)
		{

			ele = equation[j][i];

			if (ele == 0)return -9994;
			for (k = 0; k < n; k++)
			{
				equation[j][k] /= ele;
			}
		}

		for (j = i - 1; j >= 0; j--)
		{
			for (k = 0; k < n; k++)
				equation[j][k] -= equation[i][k];
		}
	}

	for (i = 0; i < m; i++)
	{
		solution[i] = equation[i][n - 1];
	}

	return 1;//success
}

void calculateBestFitPow1(int num_of_points)
{
	int i;
	for (i = 0; i < num_of_points; i++)
	{
		linearFit[i] = curvefitLinear(num_of_points, oriX[i]);
	}
}

void calculateBestFitPow2(int num_of_points)
{
	int i;
	for (i = 0; i < num_of_points; i++)
	{
		Poly2ndOrder[i] = curvefitPoly2(num_of_points, oriX[i]);
	}
}

void calculateBestFitPow3(int num_of_points)
{
	int i;
	for (i = 0; i < num_of_points; i++)
	{
		Poly3rdOrder[i] = curvefitPoly3(num_of_points, oriX[i]);
	}
}

void calculateBestFitPow4(int num_of_points)
{
	int i;
	for (i = 0; i < num_of_points; i++)
	{
		Poly4thOrder[i] = curvefitPoly4(num_of_points, oriX[i]);
	}
}

double curvefitPoly4(int n, double pointX)
{
	double max = Bo;
	double curvefitY;
	curvefitY = PolyABCDEsolution[0] * pointX * pointX * pointX * pointX + PolyABCDEsolution[1] * pointX * pointX * pointX + PolyABCDEsolution[2] * pointX * pointX + PolyABCDEsolution[3] * pointX + PolyABCDEsolution[4];
	return (curvefitY / 100.0)*max;
}

double backfitPoly4(int n, double pointY)
{
	double solution, initGuess, backfit_coeff[5];
	int i, j, num = n, flag = 0;
	double backfitX;
	double max;
	int temp;

	max = Bo;
	//Backfitting
	pointY = (pointY / max)*100.0;
	temp = pointY * 100;
	if (temp % 10 >= 5)
	{
		temp -= temp % 10;
		temp += 10;
		pointY = temp / 100.0;
	}
	else
	{
		temp -= temp % 10;
		pointY = temp / 100.0;
	}

	for (i = 0; i < num - 1; i++)
	{
		if ((pointY >= y[i] && pointY <= y[i + 1]) || (pointY <= y[i] && pointY >= y[i + 1]))
		{
			flag = 1;
			break;
		}
	}

	if (flag == 0)
	{
		//ERROR 002: y is out of bounds. Cannot Extrapolate with Interpolation formula
		//printf("Error csmain002");
		return -9997;
	}

	for (j = 0; j < 5; j++)
	{
		backfit_coeff[j] = PolyABCDEsolution[j];
	}

	backfit_coeff[4] -= pointY;

	initGuess = (x[i + 1] + x[i]) / 2;//To be improved... using Newton Raphson's method

	solution = poly_newtonRaphson(4, backfit_coeff, 0.001, initGuess);
	return solution;

}

double curvefitPoly3(int n, double pointX)
{
	double max = Bo;
	double curvefitY;
	curvefitY = PolyABCDsolution[0] * pointX * pointX * pointX + PolyABCDsolution[1] * pointX * pointX + PolyABCDsolution[2] * pointX + PolyABCDsolution[3];
	return (curvefitY / 100.0)*max;
}

double backfitPoly3(int n, double pointY)
{
	double solution, initGuess, backfit_coeff[MAXSET];
	int i, j, num = n, flag = 0;
	double backfitX;
	double max;
	int temp;

	max = Bo;
	//Backfitting
	pointY = (pointY / max)*100.0;
	temp = pointY * 100;
	if (temp % 10 >= 5)
	{
		temp -= temp % 10;
		temp += 10;
		pointY = temp / 100.0;
	}
	else
	{
		temp -= temp % 10;
		pointY = temp / 100.0;
	}

	for (i = 0; i < num - 1; i++)
	{
		if ((pointY >= y[i] && pointY <= y[i + 1]) || (pointY <= y[i] && pointY >= y[i + 1]))
		{
			flag = 1;
			break;
		}
	}

	if (flag == 0)
	{
		//ERROR 002: y is out of bounds. Cannot Extrapolate with Interpolation formula
		//printf("Error csmain002");
		return -9997;
	}

	for (j = 0; j < 4; j++)
	{
		backfit_coeff[j] = PolyABCDsolution[j];
	}

	backfit_coeff[3] -= pointY;

	initGuess = (x[i + 1] + x[i]) / 2;//To be improved... using Newton Raphson's method

	solution = poly_newtonRaphson(3, backfit_coeff, 0.001, initGuess);
	return solution;
}

double curvefitPoly2(int n, double pointX)
{
	double max = Bo;
	double solution;
	solution = PolyABCsolution[0] * pointX * pointX + PolyABCsolution[1] * pointX + PolyABCsolution[2];
	return (solution / 100.0)*max;
}

double backfitPoly2(int n, double pointY)
{
	double solution, initGuess, backfit_coeff[MAXSET];
	int i, j, num = n, flag = 0;
	double backfitX;
	double max;
	int temp;

	max = Bo;
	//Backfitting
	pointY = (pointY / max)*100.0;
	temp = pointY * 100;
	if (temp % 10 >= 5)
	{
		temp -= temp % 10;
		temp += 10;
		pointY = temp / 100.0;
	}
	else
	{
		temp -= temp % 10;
		pointY = temp / 100.0;
	}

	for (i = 0; i < num - 1; i++)
	{
		if ((pointY >= y[i] && pointY <= y[i + 1]) || (pointY <= y[i] && pointY >= y[i + 1]))
		{
			flag = 1;
			break;
		}
	}

	if (flag == 0)
	{
		//ERROR 002: y is out of bounds. Cannot Extrapolate with Interpolation formula
		//printf("Error csmain002");
		return -9997;
	}

	for (j = 0; j < 3; j++)
	{
		backfit_coeff[j] = PolyABCsolution[j];
	}

	backfit_coeff[2] -= pointY;

	initGuess = (x[i + 1] + x[i]) / 2;//To be improved... using Newton Raphson's method

	solution = poly_newtonRaphson(2, backfit_coeff, 0.001, initGuess);
	return solution;
}

int initCoeffLinear(double x[], double y[], double LinearABsolution[], int n, double PolyCoeff1[][6])
{
	double sy = 0, sx = 0, sxx = 0, syx = 0;
	int i;

	for (i = 0; i < n; i++)
	{
		sx += x[i];
		sy += y[i];
		sxx += x[i] * x[i];
		syx += x[i] * y[i];
	}

	PolyCoeff1[0][0] = sx;
	PolyCoeff1[0][1] = n;
	PolyCoeff1[0][2] = sy;
	PolyCoeff1[1][0] = sxx;
	PolyCoeff1[1][1] = sx;
	PolyCoeff1[1][2] = syx;

	if(guassJordan(PolyCoeff1, LinearABsolution, 2, 3) == -9994) return -9994;
	else return 1;
}

int initCoeffPoly2(double x[], double y[], double PolyABCsolution[], int n, double PolyCoeff2[][6])
{
	double sx = 0, sx2 = 0, sx3 = 0, sx4 = 0;
	double sy = 0, syx = 0, syx2 = 0;
	double eleX, eleY;
	int i;

	for (i = 0; i < n; i++)
	{
		eleX = x[i];
		eleY = y[i];
		sx += eleX;
		sy += eleY;
		syx += eleX * eleY;

		eleX *= x[i];
		syx2 += eleX * eleY;
		sx2 += eleX;

		eleX *= x[i];
		sx3 += eleX;

		eleX *= x[i];
		sx4 += eleX;
	}

	PolyCoeff2[0][0] = sx2;
	PolyCoeff2[0][1] = sx;
	PolyCoeff2[0][2] = n;
	PolyCoeff2[0][3] = sy;

	PolyCoeff2[1][0] = sx3;
	PolyCoeff2[1][1] = sx2;
	PolyCoeff2[1][2] = sx;
	PolyCoeff2[1][3] = syx;

	PolyCoeff2[2][0] = sx4;
	PolyCoeff2[2][1] = sx3;
	PolyCoeff2[2][2] = sx2;
	PolyCoeff2[2][3] = syx2;

	if(guassJordan(PolyCoeff2, PolyABCsolution, 3, 4) == -9994) return -9994;
	else return 1;
}

int initCoeffPoly3(double x[], double y[], double PolyABCDsolution[], int n, double PolyCoeff3[][6])
{
	double sx = 0, sx2 = 0, sx3 = 0, sx4 = 0, sx5 = 0, sx6 = 0;
	double sy = 0, syx = 0, syx2 = 0, syx3 = 0;
	double eleX, eleY;
	int i;

	for (i = 0; i < n; i++)
	{
		eleX = x[i];
		eleY = y[i];
		sx += eleX;
		sy += eleY;
		syx += eleX * eleY;

		eleX *= x[i];
		syx2 += eleX * eleY;
		sx2 += eleX;

		eleX *= x[i];
		sx3 += eleX;
		syx3 += eleY * eleX;

		eleX *= x[i];
		sx4 += eleX;

		eleX *= x[i];
		sx5 += eleX;

		eleX *= x[i];
		sx6 += eleX;
	}


	PolyCoeff3[0][0] = sx3;
	PolyCoeff3[0][1] = sx2;
	PolyCoeff3[0][2] = sx;
	PolyCoeff3[0][3] = n;
	PolyCoeff3[0][4] = sy;

	PolyCoeff3[1][0] = sx4;
	PolyCoeff3[1][1] = sx3;
	PolyCoeff3[1][2] = sx2;
	PolyCoeff3[1][3] = sx;
	PolyCoeff3[1][4] = syx;

	PolyCoeff3[2][0] = sx5;
	PolyCoeff3[2][1] = sx4;
	PolyCoeff3[2][2] = sx3;
	PolyCoeff3[2][3] = sx2;
	PolyCoeff3[2][4] = syx2;

	PolyCoeff3[3][0] = sx6;
	PolyCoeff3[3][1] = sx5;
	PolyCoeff3[3][2] = sx4;
	PolyCoeff3[3][3] = sx3;
	PolyCoeff3[3][4] = syx3;


	if(guassJordan(PolyCoeff3, PolyABCDsolution, 4, 5) == -9994) return -9994;
	else return 1;
}

int initCoeffPoly4(double x[], double y[], double PolyABCDEsolution[], int n, double PolyCoeff4[][6])
{
	double sx = 0, sx2 = 0, sx3 = 0, sx4 = 0, sx5 = 0, sx6 = 0, sx7 = 0, sx8 = 0;
	double sy = 0, syx = 0, syx2 = 0, syx3 = 0, syx4 = 0;
	double eleX, eleY;
	int i;

	for (i = 0; i < n; i++)
	{
		eleX = x[i];
		eleY = y[i];
		sx += eleX;
		sy += eleY;
		syx += eleX * eleY;

		eleX *= x[i];
		syx2 += eleX * eleY;
		sx2 += eleX;

		eleX *= x[i];
		sx3 += eleX;
		syx3 += eleY * eleX;

		eleX *= x[i];
		sx4 += eleX;
		syx4 += eleY * eleX;

		eleX *= x[i];
		sx5 += eleX;

		eleX *= x[i];
		sx6 += eleX;

		eleX *= x[i];
		sx7 += eleX;

		eleX *= x[i];
		sx8 += eleX;
	}

	PolyCoeff4[0][0] = sx4;
	PolyCoeff4[0][1] = sx3;
	PolyCoeff4[0][2] = sx2;
	PolyCoeff4[0][3] = sx;
	PolyCoeff4[0][4] = n;
	PolyCoeff4[0][5] = sy;

	PolyCoeff4[1][0] = sx5;
	PolyCoeff4[1][1] = sx4;
	PolyCoeff4[1][2] = sx3;
	PolyCoeff4[1][3] = sx2;
	PolyCoeff4[1][4] = sx;
	PolyCoeff4[1][5] = syx;

	PolyCoeff4[2][0] = sx6;
	PolyCoeff4[2][1] = sx5;
	PolyCoeff4[2][2] = sx4;
	PolyCoeff4[2][3] = sx3;
	PolyCoeff4[2][4] = sx2;
	PolyCoeff4[2][5] = syx2;

	PolyCoeff4[3][0] = sx7;
	PolyCoeff4[3][1] = sx6;
	PolyCoeff4[3][2] = sx5;
	PolyCoeff4[3][3] = sx4;
	PolyCoeff4[3][4] = sx3;
	PolyCoeff4[3][5] = syx3;

	PolyCoeff4[4][0] = sx8;
	PolyCoeff4[4][1] = sx7;
	PolyCoeff4[4][2] = sx6;
	PolyCoeff4[4][3] = sx5;
	PolyCoeff4[4][4] = sx4;
	PolyCoeff4[4][5] = syx4;

	if(guassJordan(PolyCoeff4, PolyABCDEsolution, 5, 6) == -9994) return -9994;
	else return 1;
}
