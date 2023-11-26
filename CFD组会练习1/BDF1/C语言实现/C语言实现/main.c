#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>
#include"calDet.h"

////将数组进行行列计算
//double calDet(double* arr, int n)
//{
//	int i = 0;
//	double sum = 0;
//	if (n == 1)
//	{
//		return *arr;
//	}
//	else
//	{
//		while (i < n)
//		{
//			double* arrTemp = (double*)calloc(1, sizeof(double) * (n - 1) * (n - 1));
//			double* temp = arrTemp;
//			//去除正在递归的该行列的
//			for (int k = 1; k < n; k++)
//			{
//				for (int j = 0; j < n; j++)
//				{
//					if (j != i)
//					{
//						*temp = *((arr + k * n) + j);
//						temp++;
//					}
//				}
//			}
//			sum += pow(-1, i) * (*(arr + i) * calDet(arrTemp, n - 1));//递归调用
//			i++;
//			free(arrTemp);
//		}
//		
//	}
//	return sum;
//}



//求解线性方程组
double* solveEqu(double* A, double* R, int n)
{
	double* value = (double*)malloc((n+1) * sizeof(double));
	//double value[100000000] = { 0 };
	//double* Value= (double*)calloc(1, sizeof(double) * n);
	double* Value = (double*)malloc(n * sizeof(double));
	double* Value1;
	Value1 = Value;
	value[0] = calDet(A , n);
	int i = 1;
	while (i != n + 1)
	{
		double* B= (double*)calloc(1, sizeof(double) * n*n);
		
		for (int k = 0; k < n*n; k++)
		{
			B[k] = A[k];
		}
		for (int j = 0; j < n; j++)
		{
			*(B+(i - 1)+ j*n) = *(R+j);
		}

		value[i] = calDet(B , n);
		i++;
		
	}
	//输出方程的解
	i = 1;
	while (i != n + 1)
	{
		if (value[0] == 0 && value[i] != 0)
		{
			printf("方程无解");
			break;
		}
		else
		{
			*Value = value[i] / value[0];
			Value++;
			i++;
		}
	}
	return Value1;
}


//构建矩阵
double** Make_Matrix(int row, int col)
{

	double** arr = (double**)malloc(sizeof(double*) * row);
	if (arr != NULL)
	{
		for (int i = 0; i < row; i++)
		{
			arr[i] = (double*)malloc(sizeof(double) * col);
		}
	}
	return arr;
}
//释放空间
void free_Matrix(double** src, int row)
{
	assert(src);
	for (int i = 0; i < row; i++)
	{
		free(src[i]);
	}
}

//矩阵初始化
void Init_Matrix(double** arr)
{
	int i;
	int row = (int)_msize(arr) / (int)sizeof(double*);
	for (i = 0; i < row; i++)
	{
		memset((arr[i]), 0, _msize(*arr));
	}
}

//矩阵乘法
double** Matrix_multiple(double** src,double** arr,int row1,int col1, int row2, int col2)
{
	double** c;
	c = Make_Matrix(row1, col2);//以字节存储
	Init_Matrix(c);//初始化，归零
	for (int i = 0; i < row1; i++)
	{
		
		for (int j = 0; j < col2; j++)
		{
			double sum = 0;
			for (int k = 0; k < col1; k++)
			{
				sum = sum + src[i][k] * arr[k][j];
			}
			c[i][j] = sum;
		}
	}
	return c;
}
//矩阵加法
double** Matrix_plus(double** src, double** arr, int row1, int col1, int row2, int col2)
{
	double** c;
	c = Make_Matrix(row1, col2);//以字节存储
	Init_Matrix(c);//初始化，归零
	for (int i = 0; i < row1; i++)
	{

		for (int j = 0; j < col2; j++)
		{
			c[i][j] = src[i][j]+arr[i][j];
		}
	}
	return c;
}

//矩阵转置
double** Matrix_T(double** arr)
{
	if (arr == NULL)exit(-1);
	int row = (int)_msize(arr) / (int)sizeof(double*);
	int col = (int)_msize(*arr) / (int)sizeof(double);
	double** T = (double**)malloc(sizeof(double*) * col);
	int i = 0;
	int j = 0;
	if (T != NULL)
	{
		for (i = 0; i < col; i++)
		{
			T[i] = (double*)malloc(sizeof(double) * row);
		}
	}
	for (i = 0; i < col; i++)
	{
		for (j = 0; j < row; j++)
		{
			T[i][j] = arr[j][i];
		}
	}
	return T;
}

int main()
{
	//计算运行时间
	clock_t start, end;/*首先用clock_t定义两个变量来存储开始与结束的值*/

	start = clock();/*记录开始的值*/
//Pre-proceeding
const double Pi=3.1415926;
const double Unit=64;
const double CFL=0.01;
const double endtau=1;
const double endx=1;
const double tol=1/pow(10,10);
const double nu=1;
const double deltax=endx/Unit;
const int numberx = (int)(endx / deltax)+ 1;
const double Lr = 1 / (2 * Pi);
const double Tr = pow(Lr, 2) / nu;
const double abslambda = sqrt(nu / Tr);
const double deltatau = CFL * deltax / abslambda;
const double B1 = 1;
//double C[2][2] = { B1,0,0,B1 / deltax };
double** C; double** C1; double** C2;
C = Make_Matrix(2, 2); C1 = Make_Matrix(2, 2); C2 = Make_Matrix(2, 2);
C[0][0] = B1; C[0][1] = 0; C[1][0] = 0; C[1][1] = B1 / deltax;
C1[0][0] = abslambda / 2; C1[0][1] = -nu / (2 * deltax); C1[1][0] = -1 / (2 * Tr); C1[1][1] = abslambda / (2 * deltax);
C2[0][0] = -abslambda / 2; C2[0][1] = -nu / (2 * deltax); C2[1][0] = -1 / (2 * Tr); C2[1][1] = -abslambda / (2 * deltax);

double Mtau[2][2] = { deltax,0,0,1 / deltax };
double A[2][2] = { abslambda,0,0,abslambda };

double* Uexasolution1, * q1, * Uexasolution2, * q2;
Uexasolution1 = (double*)malloc(numberx*sizeof(double));//以字节存储
Uexasolution2 = (double*)malloc(numberx * sizeof(double));
//if (Uexasolution1 == NULL)
//{
//	printf("Not able to alloc memory.\n");
//	exit(1);
//}
q1 = Uexasolution1;
q2 = Uexasolution2;



//Proceeding
//solve the exasolution
double x = 0;
for (int i=0;i<numberx ; i++)
{
	*(Uexasolution1+i) = sin(Pi * x);
	*(Uexasolution2 + i) = Pi * cos(Pi * x);
	x += deltax;
}

//构建LHS
double** LHS1, ** LHS2,** LHS3,**LHS;
LHS1 = Make_Matrix(2*Unit, 2*Unit);//以字节存储
LHS2 = Make_Matrix(2 * Unit, 2 * Unit);
LHS3 = Make_Matrix(2 * Unit, 2 * Unit);
LHS = Make_Matrix(2 * Unit, 2 * Unit);
Init_Matrix(LHS1); Init_Matrix(LHS2); Init_Matrix(LHS3); Init_Matrix(LHS);
//构建LHS1
for (int i = 0; i < 2 * Unit; i+=2)
{

	int j = i;
		LHS1[i][j] = deltax / deltatau;
	
}
for (int i = 1; i < 2 * Unit; i += 2)
{
	int j = i;
	LHS1[i][j] = 1 / (deltax * deltatau);
}

//构建Rdomain
for (int i = 1; i < 2 * Unit; i += 2)
{
	int j = i;
	LHS2[i][j] = 1 / (Tr*deltax);
}

//Rboundary
for (int iface = 2; iface <= numberx - 1; iface++)
{
	int ieL = iface - 1;
	int ieR = iface;
	//diag
	for (int i = 2*(ieL - 1); i < 2*ieL; i++)
	{
		for (int j = 2 * (ieL - 1); j < 2 * ieL; j++)
			if (i == j && i == 2*(ieL - 1))
			{
				LHS3[i][j] = LHS3[i][j] + Matrix_multiple(Matrix_T(C), C1, 2, 2, 2, 2)[0][0];
			}
			else if (i == j && i == 2 * ieL - 1)
			{
				LHS3[i][j] = LHS3[i][j] + Matrix_multiple(Matrix_T(C), C1, 2, 2, 2, 2)[1][1];
			}
			else if (i != j)
			{
				switch (i - j)
				{
				case -1:
					LHS3[i][j] = LHS3[i][j] + Matrix_multiple(Matrix_T(C), C1, 2, 2, 2, 2)[0][1];
					break;
				case 1:
					LHS3[i][j] = LHS3[i][j] + Matrix_multiple(Matrix_T(C), C1, 2, 2, 2, 2)[1][0];
					break;
				default:
					break;
				}
			}
			
			
	}

	for (int i = 2 * (ieR - 1); i < 2 * ieR; i++)
	{
		for (int j = 2 * (ieR - 1); j < 2 * ieR; j++)
			if (i == j && i == 2 * (ieR - 1))
			{
				LHS3[i][j] = LHS3[i][j] - Matrix_multiple(Matrix_T(C), C2, 2, 2, 2, 2)[0][0];
			}
			else if (i == j && i == 2 * ieR-1)
			{
				LHS3[i][j] = LHS3[i][j] - Matrix_multiple(Matrix_T(C), C2, 2, 2, 2, 2)[1][1];
			}
			else if (i != j)
			{
				switch (i - j)
				{
				case -1:
					LHS3[i][j] = LHS3[i][j] - Matrix_multiple(Matrix_T(C), C2, 2, 2, 2, 2)[0][1];
					break;
				case 1:
					LHS3[i][j] = LHS3[i][j] - Matrix_multiple(Matrix_T(C), C2, 2, 2, 2, 2)[1][0];
					break;
				default:
					break;
				}
			}


	}
	//upper
	for (int i = 2 * (ieL - 1); i < 2 * ieL; i++)
	{
		for (int j = 2 * (ieR - 1); j < 2 * ieR; j++)
			if (i == (j-2) && i == 2 * (ieL - 1))
			{
				LHS3[i][j] = LHS3[i][j] + Matrix_multiple(Matrix_T(C), C2, 2, 2, 2, 2)[0][0];
			}
			else if (i == (j-2) && i == 2 * ieL - 1)
			{
				LHS3[i][j] = LHS3[i][j] + Matrix_multiple(Matrix_T(C), C2, 2, 2, 2, 2)[1][1];
			}
			else if (i != (j-2))
			{
				switch (i - j)
				{
				case -3:
					LHS3[i][j] = LHS3[i][j] + Matrix_multiple(Matrix_T(C), C2, 2, 2, 2, 2)[0][1];
					break;
				case -1:
					LHS3[i][j] = LHS3[i][j] + Matrix_multiple(Matrix_T(C), C2, 2, 2, 2, 2)[1][0];
					break;
				default:
					break;
				}
			}


	}
	//lower
	for (int i = 2 * (ieR - 1); i < 2 * ieR; i++)
	{
		for (int j = 2 * (ieL - 1); j < 2 * ieL; j++)
			if (i == (j + 2) && i == 2 * (ieR - 1))
			{
				LHS3[i][j] = LHS3[i][j] - Matrix_multiple(Matrix_T(C), C1, 2, 2, 2, 2)[0][0];
			}
			else if (i == (j + 2) && i == 2 * ieR - 1)
			{
				LHS3[i][j] = LHS3[i][j] - Matrix_multiple(Matrix_T(C), C1, 2, 2, 2, 2)[1][1];
			}
			else if (i != (j + 2))
			{
				switch (i - j)
				{
				case 1:
					LHS3[i][j] = LHS3[i][j] - Matrix_multiple(Matrix_T(C), C1, 2, 2, 2, 2)[0][1];
					break;
				case 3:
					LHS3[i][j] = LHS3[i][j] - Matrix_multiple(Matrix_T(C), C1, 2, 2, 2, 2)[1][0];
					break;
				default:
					break;
				}
			}


	}

}
//边界值
double** C3;
C3 = Make_Matrix(2, 2);
C3[0][0] = 0; C3[0][1] = 0; C3[1][0] = 0; C3[1][1] = 1;
double** C4;
C4 = Make_Matrix(2, 2);
C4[0][0] = abslambda / 2; C4[0][1] = -nu / 2 ; C4[1][0] = -1 / (2 * Tr); C4[1][1] = abslambda / 2 ;
double** C5;
C5 = Make_Matrix(2, 2);
C5[0][0] = -abslambda / 2; C5[0][1] = -nu / 2; C5[1][0] = -1 / (2 * Tr); C5[1][1] = -abslambda / 2;

for (int i = 2 * (1 - 1); i < 2 * 1; i++)
{
	for (int j = 2 * (1 - 1); j < 2 * 1; j++)
		if (i == j && i == 2 * (1 - 1))
		{
			LHS3[i][j] = LHS3[i][j] - Matrix_multiple(Matrix_multiple(Matrix_T(C),Matrix_plus(Matrix_multiple(C4,C3,2,2,2,2),C5,2,2,2,2), 2, 2, 2, 2),C,2,2,2,2)[0][0];
		}
		else if (i == j && i == 2 * 1 - 1)
		{
			LHS3[i][j] = LHS3[i][j] - Matrix_multiple(Matrix_multiple(Matrix_T(C), Matrix_plus(Matrix_multiple(C4, C3, 2, 2, 2, 2), C5, 2, 2, 2, 2), 2, 2, 2, 2), C, 2, 2, 2, 2)[1][1];
		}
		else if (i != j)
		{
			switch (i - j)
			{
			case -1:
				LHS3[i][j] = LHS3[i][j] - Matrix_multiple(Matrix_multiple(Matrix_T(C), Matrix_plus(Matrix_multiple(C4, C3, 2, 2, 2, 2), C5, 2, 2, 2, 2), 2, 2, 2, 2), C, 2, 2, 2, 2)[0][1];;
				break;
			case 1:
				LHS3[i][j] = LHS3[i][j] - Matrix_multiple(Matrix_multiple(Matrix_T(C), Matrix_plus(Matrix_multiple(C4, C3, 2, 2, 2, 2), C5, 2, 2, 2, 2), 2, 2, 2, 2), C, 2, 2, 2, 2)[1][0];;
				break;
			default:
				break;
			}
		}


}
for (int i = 2 * (numberx - 1 - 1); i < 2 *(numberx-1); i++)
{
	for (int j = 2 * (numberx - 1 - 1); j < 2 * (numberx - 1); j++)
		if (i == j && i == 2 * (numberx - 1 - 1))
		{
			LHS3[i][j] = LHS3[i][j] + Matrix_multiple(Matrix_multiple(Matrix_T(C),Matrix_plus(C4,Matrix_multiple(C5,C3,2,2,2,2),2,2,2,2),2,2,2,2),C,2,2,2,2)[0][0];
		}
		else if (i == j && i == 2 * (numberx - 1) - 1)
		{
			LHS3[i][j] = LHS3[i][j] + Matrix_multiple(Matrix_multiple(Matrix_T(C), Matrix_plus(C4, Matrix_multiple(C5, C3, 2, 2, 2, 2), 2, 2, 2, 2), 2, 2, 2, 2), C, 2, 2, 2, 2)[1][1];
		}
		else if (i != j)
		{
			switch (i - j)
			{
			case -1:
				LHS3[i][j] = LHS3[i][j] + Matrix_multiple(Matrix_multiple(Matrix_T(C), Matrix_plus(C4, Matrix_multiple(C5, C3, 2, 2, 2, 2), 2, 2, 2, 2), 2, 2, 2, 2), C, 2, 2, 2, 2)[0][1];;
				break;
			case 1:
				LHS3[i][j] = LHS3[i][j] + Matrix_multiple(Matrix_multiple(Matrix_T(C), Matrix_plus(C4, Matrix_multiple(C5, C3, 2, 2, 2, 2), 2, 2, 2, 2), 2, 2, 2, 2), C, 2, 2, 2, 2)[1][0];;
				break;
			default:
				break;
			}
		}


}
//计算LHS
for (int i = 0; i < 2 * Unit; i++)
{
	for (int j = 0; j < 2 * Unit; j++)
	{
		LHS[i][j]=LHS1[i][j]+ LHS2[i][j] + LHS3[i][j];
	}
}

//构建所需要的D
double** D;
D = Make_Matrix(2 * Unit, 2 * Unit);
Init_Matrix(D);
for (int iface = 2; iface <= numberx; iface++)
{
	int ieL = iface - 1;
	D[2 * (ieL - 1)][2 * (ieL - 1)] = LHS[2 * (ieL - 1)][2 * (ieL - 1)];
	D[2 * (ieL - 1)][2 * (ieL - 1)+1] = LHS[2 * (ieL - 1)][2 * (ieL - 1)+1];
	D[2 * (ieL - 1)+1][2 * (ieL - 1)] = LHS[2 * (ieL - 1)+1][2 * (ieL - 1)];
	D[2 * (ieL - 1)+1][2 * (ieL - 1)+1] = LHS[2 * (ieL - 1)+1][2 * (ieL - 1)+1];

}
double** Ucurrent; double** Unext; double** Rd; double** Rb; double** Fn; double** R;
Ucurrent = Make_Matrix(2, numberx - 1); Unext = Make_Matrix(2 * Unit, 1); Rd = Make_Matrix(2, numberx - 1); Rb = Make_Matrix(2, numberx - 1); Fn = Make_Matrix(2, numberx); R = Make_Matrix(2*Unit, 1);
Init_Matrix(Ucurrent); Init_Matrix(Unext); Init_Matrix(Rb); Init_Matrix(Rd); Init_Matrix(R);

 x = 0;
for (int k = 0; k < numberx - 1; k++)
{

	Ucurrent[0][k]= (x + deltax / 2)* (x + deltax / 2) - (x + deltax / 2);
	Ucurrent[1][k]= (2 * (x + deltax / 2) - 1) * deltax;
	x += deltax;

}

//Rdomain
x = 0;
for (int k = 0; k < numberx - 1; k++)
{

	Rd[0][k] = Pi*(cos(Pi*x)-cos(Pi*(x+deltax)));
	Rd[1][k] = -Ucurrent[1][k]/(Tr*deltax);
	x = x + deltax;

}
//Rboundary
for (int iface = 2; iface <= numberx - 1; iface++)
{
	int ieL = iface - 1;
	int ieR = iface;
	Fn[0][iface - 1] = 0.5 * (-nu * Ucurrent[1][ieL - 1] / deltax - nu * Ucurrent[1][ieR - 1] / deltax) - 0.5 * (A[0][0] * (Ucurrent[0][ieR - 1] - Ucurrent[0][ieL - 1]) + A[0][1] * (Ucurrent[1][ieR - 1] / deltax - Ucurrent[1][ieL - 1] / deltax));
	Fn[1][iface - 1] = 0.5 * (-Ucurrent[0][ieL - 1] / Tr - Ucurrent[0][ieR - 1] / Tr) - 0.5 * (A[1][0] * (Ucurrent[0][ieR - 1] - Ucurrent[0][ieL - 1]) + A[1][1] * (Ucurrent[1][ieR - 1] / deltax - Ucurrent[1][ieL - 1] / deltax));
	Rb[0][ieL - 1] = Rb[0][ieL - 1] - Fn[0][iface - 1]; Rb[1][ieL - 1] = Rb[1][ieL - 1] - (1/deltax)*Fn[1][iface - 1];
	Rb[0][ieR - 1] = Rb[0][ieR - 1] + Fn[0][iface - 1]; Rb[1][ieR - 1] = Rb[1][ieR - 1] + (1 / deltax) * Fn[1][iface - 1];

}
Fn[0][0]= 0.5 * (-nu * Ucurrent[1][0] / deltax - nu * Ucurrent[1][0] / deltax) - 0.5 * (A[0][0] * (Ucurrent[0][0] - 0) + A[0][1] * (Ucurrent[1][0] / deltax - Ucurrent[1][0] / deltax));
Fn[1][0] = 0.5 * (0 - Ucurrent[0][0] / Tr) - 0.5 * (A[1][0] * (Ucurrent[0][0] - 0) + A[1][1] * (Ucurrent[1][0] / deltax - Ucurrent[1][0] / deltax));
Fn[0][numberx - 1] = 0.5 * (-nu * Ucurrent[1][numberx-2] / deltax - nu * Ucurrent[1][numberx - 2] / deltax) - 0.5 * (A[0][0] * (0 - Ucurrent[0][numberx-2]) + A[0][1] * (Ucurrent[1][numberx-2] / deltax - Ucurrent[1][numberx-2] / deltax));
Fn[1][numberx - 1] = 0.5 * (-Ucurrent[0][numberx - 2] / Tr - 0) - 0.5 * (A[1][0] * (0 - Ucurrent[0][numberx - 2]) + A[1][1] * (Ucurrent[1][numberx - 2] / deltax - Ucurrent[1][numberx - 2] / deltax));
Rb[0][0] = Rb[0][0] + Fn[0][0];
Rb[1][0] = Rb[1][0] + Fn[1][0]/deltax;
Rb[0][numberx-2] = Rb[0][numberx-2] - Fn[0][numberx-1];
Rb[1][numberx-2] = Rb[1][numberx-2] - Fn[1][numberx-1] / deltax;

//R组装

for (int k = 1; k <= numberx - 1; k++)
{
	R[2 * k - 2][0] = Rd[0][k - 1] + Rb[0][k - 1];
	R[2 * k-1][0] = Rd[1][k - 1] + Rb[1][k - 1];

}
//进行必要的向量转换
for (int k = 1; k <= numberx - 1; k++)
{
	Unext[2 * k - 2][0] = Ucurrent[0][k - 1];
	Unext[2 * k - 1][0] = Ucurrent[1][k - 1];
}




//循环
int flag = 0; double no = deltatau; double Vbreak = 0;

for (double n = deltatau; n <= endtau; n += deltatau)
{
	//求V=D\R
	double* V1; double* V;
	V = (double*)malloc(2 * Unit * sizeof(double));
	for (int i = 0; i < 2 * Unit; i++)
	{

		if (i % 2 == 1)
		{
			double D2[4]; double R2[2];
			D2[0] = D[i - 1][i - 1];
			D2[1] = D[i - 1][i];
			D2[2] = D[i][i - 1];
			D2[3] = D[i][i];
			R2[0] = R[i - 1][0];
			R2[1] = R[i][0];
			V1 = solveEqu(D2, R2, 2);
			V[i - 1] = V1[0];
			V[i] = V1[1];
		}
	}
	Vbreak = 0;
	for (int N = 0; N < 2 * Unit; N++)
	{
		Vbreak += V[N];
		
	}
	if (Vbreak < tol)
	{
		flag = 1;
	}

	if (flag == 1)
	{
		break;
	}

	for (int i = 0; i < 2 * Unit; i++)
	{
		Unext[i][0] += V[i];
	}


	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < numberx - 1; j++)
		{
			Rb[i][j] = 0;
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < numberx - 1; j++)
		{
			Rd[i][j] = 0;
		}
	}
	


	for (int i = 1; i <= numberx - 1; i++)
	{
		 Ucurrent[0][i - 1]=Unext[2 * i - 2][0] ;
		 Ucurrent[1][i - 1]=Unext[2 * i - 1][0] ;
	}

	//Rdomain
	x = 0;
	for (int k = 0; k < numberx - 1; k++)
	{

		Rd[0][k] = Pi * (cos(Pi * x) - cos(Pi * (x + deltax)));
		Rd[1][k] = -Ucurrent[1][k] / (Tr * deltax);
		x = x + deltax;

	}
	//Rboundary
	for (int iface = 2; iface <= numberx - 1; iface++)
	{
		int ieL = iface - 1;
		int ieR = iface;
		Fn[0][iface - 1] = 0.5 * (-nu * Ucurrent[1][ieL - 1] / deltax - nu * Ucurrent[1][ieR - 1] / deltax) - 0.5 * (A[0][0] * (Ucurrent[0][ieR - 1] - Ucurrent[0][ieL - 1]) + A[0][1] * (Ucurrent[1][ieR - 1] / deltax - Ucurrent[1][ieL - 1] / deltax));
		Fn[1][iface - 1] = 0.5 * (-Ucurrent[0][ieL - 1] / Tr - Ucurrent[0][ieR - 1] / Tr) - 0.5 * (A[1][0] * (Ucurrent[0][ieR - 1] - Ucurrent[0][ieL - 1]) + A[1][1] * (Ucurrent[1][ieR - 1] / deltax - Ucurrent[1][ieL - 1] / deltax));
		Rb[0][ieL - 1] = Rb[0][ieL - 1] - Fn[0][iface - 1]; Rb[1][ieL - 1] = Rb[1][ieL - 1] - (1 / deltax) * Fn[1][iface - 1];
		Rb[0][ieR - 1] = Rb[0][ieR - 1] + Fn[0][iface - 1]; Rb[1][ieR - 1] = Rb[1][ieR - 1] + (1 / deltax) * Fn[1][iface - 1];

	}
	Fn[0][0] = 0.5 * (-nu * Ucurrent[1][0] / deltax - nu * Ucurrent[1][0] / deltax) - 0.5 * (A[0][0] * (Ucurrent[0][0] - 0) + A[0][1] * (Ucurrent[1][0] / deltax - Ucurrent[1][0] / deltax));
	Fn[1][0] = 0.5 * (0 - Ucurrent[0][0] / Tr) - 0.5 * (A[1][0] * (Ucurrent[0][0] - 0) + A[1][1] * (Ucurrent[1][0] / deltax - Ucurrent[1][0] / deltax));
	Fn[0][numberx - 1] = 0.5 * (-nu * Ucurrent[1][numberx - 2] / deltax - nu * Ucurrent[1][numberx - 2] / deltax) - 0.5 * (A[0][0] * (0 - Ucurrent[0][numberx - 2]) + A[0][1] * (Ucurrent[1][numberx - 2] / deltax - Ucurrent[1][numberx - 2] / deltax));
	Fn[1][numberx - 1] = 0.5 * (-Ucurrent[0][numberx - 2] / Tr - 0) - 0.5 * (A[1][0] * (0 - Ucurrent[0][numberx - 2]) + A[1][1] * (Ucurrent[1][numberx - 2] / deltax - Ucurrent[1][numberx - 2] / deltax));
	Rb[0][0] = Rb[0][0] + Fn[0][0];
	Rb[1][0] = Rb[1][0] + Fn[1][0] / deltax;
	Rb[0][numberx - 2] = Rb[0][numberx - 2] - Fn[0][numberx - 1];
	Rb[1][numberx - 2] = Rb[1][numberx - 2] - Fn[1][numberx - 1] / deltax;

	//R组装

	for (int k = 1; k <= numberx - 1; k++)
	{
		R[2 * k - 2][0] = Rd[0][k - 1] + Rb[0][k - 1];
		R[2 * k - 1][0] = Rd[1][k - 1] + Rb[1][k - 1];

	}
	//进行必要的向量转换
	for (int k = 1; k <= numberx - 1; k++)
	{
		Unext[2 * k - 2][0] = Ucurrent[0][k - 1];
		Unext[2 * k - 1][0] = Ucurrent[1][k - 1];
	}
	no += deltatau;
}

for (int k = 1; k <= numberx - 1; k++)
{
	 Ucurrent[1][k - 1]= Ucurrent[1][k - 1]/deltax;

}



for (int i = 0; i < numberx-1; i++)
{
	printf("%lf\n", Ucurrent[1][i]);


}


free(q1);
free(q2);
free_Matrix(LHS1, 2 * Unit);
free_Matrix(LHS2, 2 * Unit);
free_Matrix(LHS3, 2 * Unit);
free_Matrix(LHS, 2 * Unit);
free_Matrix(C, 2); free_Matrix(C1, 2); free_Matrix(C2, 2); free_Matrix(C3, 2); free_Matrix(C4, 2); free_Matrix(C5, 2);

//计算运行时间
//clock_t start, end;/*首先用clock_t定义两个变量来存储开始与结束的值*/
//
//start = clock();/*记录开始的值*/

end = clock(); /* 记录结束的值*/

double shijian = ((double)(end - start)) / CLK_TCK;/*用结束时间减去开始时间,因为他是毫秒单位为此除以CLK_TCK来转化为秒*/

printf("\n运行的时间是%lf", shijian);

return 0;
}