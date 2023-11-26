#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>
//将数组进行行列计算
double calDet(double* arr, int n)
{
	int i = 0;
	double sum = 0;
	if (n == 1)
	{
		return *arr;
	}
	else
	{
		while (i < n)
		{
			double* arrTemp = (double*)calloc(1, sizeof(double) * (n - 1) * (n - 1));
			double* temp = arrTemp;
			//去除正在递归的该行列的
			for (int k = 1; k < n; k++)
			{
				for (int j = 0; j < n; j++)
				{
					if (j != i)
					{
						*temp = *((arr + k * n) + j);
						temp++;
					}
				}
			}
			sum += pow(-1, i) * (*(arr + i) * calDet(arrTemp, n - 1));//递归调用
			i++;
			free(arrTemp);
		}

	}
	return sum;
}