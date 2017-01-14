#include <stdio.h>

#define MAXSET 10

/*
INTERPOLATION LIST
Cubic Spline Interpolation

ERROR LIST
-9997 : Point not in range.
-9996 : Number of points less than 4 (or greater than MAXSET)
-9995 : Solution not converging.
*/

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

//---------------------------------------//
//VARIABLES
//INPUTS (undisturbed)
double oriX[MAXSET] = { 1, 2, 3, 4, 5 }, oriY[MAXSET] = { 100, 200, 300, 400, 500 };
int num_of_points = 5;
double pointY = 350, pointX = 3.5;

//COEFFICIENTS
double CSIcoeff[MAXSET * 4 - 4]; //cubicSpline

//Y-values (interpolated)
double CSIinterpolation[MAXSET];

//Others essential variables
double Bo; //Max value on Y-axis
double x[MAXSET], y[MAXSET]; //Used variables

int main()
{
	double backfitX, curvefitY;

	//STEP 1: initialize the values
	initInterpolation();

	//STEP 2: calculate coeff (for cubic spline interpolation only)
	cubicSplineCoeff(x, y, num_of_points);

	//Step 3: Calculate interpolated values.
	calculateCSIinterpolation(num_of_points);

	//Step 4: Backfitting
	//Cubic Spline Interpolation
	backfitX = backfitCSI(num_of_points, pointY);
	//Curve fitting
	curvefitY = curvefitCSI(num_of_points, pointX);

	printf("Backfit (Cubic Spline Interpolation): %.3lf\n", backfitX);
	printf("Curvefit (Cubic Spline Interpolation): %.3lf\n", curvefitY);

	return 0;
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
	initGuess = (x[splineNum + 1] - x[splineNum]) / 2;//To be improved...

	for (i = 0; i < 4; i++)
	{
		backfit_coeff[3 - i] = *(a + i + splineNum * 4);
	}

	backfit_coeff[3] -= pointY;
	backfitX = poly_newtonRaphson(3, backfit_coeff, accuracy, initGuess);

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

void calculateCSIinterpolation(int num_of_points)
{
	int i;
	for (i = 0; i < num_of_points; i++)
	{
		CSIinterpolation[i] = curvefitCSI(num_of_points, oriX[i]);
	}
}
