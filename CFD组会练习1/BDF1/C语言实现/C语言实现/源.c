#include<stdio.h>
#include<math.h>
//double power(double x, int n)
//{
//	double num = 1;
//	for (int i = 1; i <= n; i++)
//	{
//		num = x * num;
//	}
//	return num;
//}

int main()
{
	//Pre-proceeding
	const double Pi = 3.1415926;
	const double Unit = 64;
	const double CFL = 0.01;
	const double endtau = 1;
	const double endx = 1;
	const double tol = 1 / pow(10, 10);
	const double nu = 1;
	const double deltax = endx / Unit;
	const int numberx = endx / deltax + 1;
	const double Lr = 1 / (2 * Pi);
	const double Tr = pow(Lr, 2) / nu;
	const double abslambda = sqrt(nu / Tr);
	const double deltatau = CFL * deltax / abslambda;
	const double B1 = 1;
	double C[2][2] = { B1,0,0,B1 / deltax };
	double Mtau[2][2] = { deltax,0,0,1 / deltax };
	double A[2][2] = { abslambda,0,0,abslambda };

	double Uexasolution1[10000];
	double Uexasolution2[10000];


	//Proceeding
	//solve the exasolution
	double x = 0;
	for (int k = 0; k < numberx; k++)
	{
		Uexasolution1[k] = sin(Pi * x);
		//Uexasolution2[k] = Pi * cos(Pi * x);
		x += deltax;
	}
	for (int k = 0; k < numberx; k++)
	{
		printf("%lf\n", Uexasolution1[k]);
	}

	return 0;
}