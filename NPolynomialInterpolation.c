#include <stdio.h>

#define MAXSET 10

/*
INTERPOLATION LIST
Neville's polynomial Interpolation

ERROR LIST
-9997 : Point not in range.
-9996 : Number of points less than 4 (or greater than MAXSET)
-9995 : Solution not converging.
*/

//Initialization
void initInterpolation();

//NEVILLE'S POLYNOMIAL INTERPOLATION
double neville_poly(int num_of_points, double y[], double x[], double pointXY);
double curvefitNPoly(int num_of_points, double pointX);
double backfitNPoly(int num_of_points, double pointY);
void calculateNpolyInterpolation(int num_of_points);

//---------------------------------------//
//VARIABLES
//INPUTS (undisturbed)
double oriX[MAXSET] = { 1, 2, 3, 4, 5 }, oriY[MAXSET] = { 100, 200, 300, 400, 500 };
int num_of_points = 5;
double pointY = 360, pointX = 3.56;

//Y-values (interpolated)
double NpolyInterpolation[MAXSET];

//Others essential variables
double Bo; //Max value on Y-axis
double x[MAXSET], y[MAXSET]; //Used variables

int main()
{
	double backfitX, curvefitY;

	//STEP 1: initialize the values
	initInterpolation();

	//No Step 2 (No calculations for coeff)

	//Step 3: Calculate interpolated values.
	calculateNpolyInterpolation(num_of_points);

	//Step 4: Backfitting
	//Neville's polynomial
	//backfitting
	backfitX = backfitNPoly(num_of_points, pointY);
	//Curve fitting
	curvefitY = curvefitNPoly(num_of_points, pointX);

	printf("Backfit (Neville's Polynomial Interpolation): %.3lf\n", backfitX);
	printf("Curvefit (Neville's Polynomial Interpolation): %.3lf\n", curvefitY);

	return 0;
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

void calculateNpolyInterpolation(int num_of_points)
{
	int i;
	for (i = 0; i < num_of_points; i++)
	{
		NpolyInterpolation[i] = curvefitNPoly(num_of_points, oriX[i]);
	}
}
