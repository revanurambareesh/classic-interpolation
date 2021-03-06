#include <stdio.h>
#define MAXSET 10

//Initialization
void initInterpolation();

//POLY2 FIT
void calculateBestFitPow2(int num_of_points);
double curvefitPoly2(int n, double pointX);
double backfitPoly2(int n, double pointY);
int guassJordan(double equation[][6], double solution[], int m, int n);
double poly_newtonRaphson(int power, double coeff[], float accuracy, double initApprox);
double poly_fvalue(int power, double coeff[], double x);
int initCoeffPoly2(double x[], double y[], double PolyABCsolution[], int n, double PolyCoeff2[][6]);
double poly_f1dash(int power, double coeff[], double x);

//VARIABLES
//INPUTS (undisturbed)
double oriX[MAXSET] = { 0, 1, 2, 3, 4 }, oriY[MAXSET] = { 4, 5, 8, 13, 20 };
int num_of_points = 5;
double pointY = 7, pointX = 4;

//COEFFICIENTS
double PolyCoeff2[3][6];//Best fit polynomial2

//Y-values (interpolated)
double Poly2ndOrder[MAXSET];

//Others essential variables
double Bo; //Max value on Y-axis
double x[MAXSET], y[MAXSET]; //Used variables
double PolyABCsolution[3];

int main()
{
	double backfitX, curvefitY;

	//STEP 1: initialize the values
	initInterpolation();

	//STEP 2: calculate coeff
	if(initCoeffPoly2(x, y, PolyABCsolution, num_of_points, PolyCoeff2) == -9994)printf("Cannot Implement");

	//Step 3: Calculate interpolated values.
	calculateBestFitPow2(num_of_points);

	//Step 4: Backfitting
	//Polynomial 2nd order
	backfitX = backfitPoly2(num_of_points, pointY);
	//Curve fitting
	curvefitY = curvefitPoly2(num_of_points, pointX);

	printf("Backfit (Polynomial 2nd order): %.3lf\n", backfitX);
	printf("Curvefit (Polynomial 2nd order): %.3lf\n", curvefitY);

	return 0;
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

void calculateBestFitPow2(int num_of_points)
{
	int i;
	for (i = 0; i < num_of_points; i++)
	{
		Poly2ndOrder[i] = curvefitPoly2(num_of_points, oriX[i]);
	}
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
