#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>
//������������м���
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
			//ȥ�����ڵݹ�ĸ����е�
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
			sum += pow(-1, i) * (*(arr + i) * calDet(arrTemp, n - 1));//�ݹ����
			i++;
			free(arrTemp);
		}

	}
	return sum;
}