#include <stdio.h>
#define MAXSET 10

//Initialization
void initInterpolation();

//LINEAR FIT
double curvefitLinear(int n, double pointX);
double backfitLinear(int n, double pointY);
int guassJordan(double equation[][6], double solution[], int m, int n);
void calculateBestFitPow1(int num_of_points);
double poly_newtonRaphson(int power, double coeff[], float accuracy, double initApprox);
double poly_fvalue(int power, double coeff[], double x);
int initCoeffLinear(double x[], double y[], double LinearABsolution[], int n, double PolyCoeff1[][6]);
double poly_f1dash(int power, double coeff[], double x);


//VARIABLES
//INPUTS (undisturbed)
double oriX[MAXSET] = { 0, 1, 2, 3, 4 }, oriY[MAXSET] = { 1, 2, 3, 4, 5 };
int num_of_points = 5;
double pointY = 7, pointX = 4;

//COEFFICIENTS
double PolyCoeff1[2][6];//Best fit -Linear

//Y-values (interpolated)
double linearFit[MAXSET];

//Others essential variables
double Bo; //Max value on Y-axis
double x[MAXSET], y[MAXSET]; //Used variables
double LinearABsolution[2];

int main()
{
	double backfitX, curvefitY;

	//STEP 1: initialize the values
	initInterpolation();

	//STEP 2: calculate coeff
	if(initCoeffLinear(x, y, LinearABsolution, num_of_points, PolyCoeff1) == -9994)printf("Cannot Implement");

	//Step 3: Calculate interpolated values.
	calculateBestFitPow1(num_of_points);

	//Step 4: Backfitting

	//LinearCurveFit
	backfitX = backfitLinear(num_of_points, pointY);
	//Curve fitting
	curvefitY = curvefitLinear(num_of_points, pointX);

	printf("Backfit (Linear Best Fit): %.3lf\n", backfitX);
	printf("Curvefit (Linear Best Fit): %.3lf\n", curvefitY);

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

void calculateBestFitPow1(int num_of_points)
{
	int i;
	for (i = 0; i < num_of_points; i++)
	{
		linearFit[i] = curvefitLinear(num_of_points, oriX[i]);
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
