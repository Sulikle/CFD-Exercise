#include <stdio.h>
#include <math.h>
#define MAX 20

void swap(double* a, double* b) {
    double tmp = *a;
    *a = *b;
    *b = tmp;
}

void swapArray(double a[MAX], double b[MAX], int n) {
    int i;
    for (i = 0; i < n; i++)
        swap(&a[i], &b[i]);
}

double computeDet(double array[MAX][MAX], int n) {
    // 通过初等行变换，再计算对角元素乘积，来计算行列式的值

    int i, j, k;
    int sign = 0;             //行列式交换一次需要改变符号，此变量记录交换次数 
    double sum = 1.0;         //结果
    double tmp;             //暂存乘积因子
    double zero = 1e-6;       //两浮点数差距小于1e-6视为相等

    for (i = 0; i < n - 1; i++) {
        k = 1;
        while (fabs(array[i][i] - 0.0) < zero && i + k < n) {
            swapArray(array[i], array[i + (k++)], n);
            sign++;
        }

        for (j = i + 1; j < n; j++) {
            if (array[j][i] == 0.0)			//如为0则那一行不用化简 
                continue;
            else {
                tmp = -(double)array[j][i] / array[i][i];		//保存乘积因子 
                for (k = i; k < n; k++)
                    array[j][k] += (tmp * array[i][k]);
            }
        }
    }

    for (i = 0; i < n; i++)
        sum *= array[i][i];
    if (sign % 2 != 0)				//交换偶数次符仍为正 
        sum *= -1;
    return sum;
}

void rowTrans(double array[MAX][MAX], int n) {
    // 做初等行变换
    int i, j, k;
    double tmp;             //暂存乘积因子
    double zero = 1e-6;       //两浮点数差距小于1e-6视为相等

    for (i = 0; i < n - 1; i++) {
        k = 1;
        while (fabs(array[i][i] - 0.0) < zero && i + k < n) {
            swapArray(array[i], array[i + (k++)], n);
        }

        for (j = i + 1; j < n; j++) {
            if (array[j][i] == 0.0)			//如为0则那一行不用化简 
                continue;
            else {
                tmp = -(double)array[j][i] / array[i][i];		//保存乘积因子 
                for (k = i; k < n; k++)
                    array[j][k] += (tmp * array[i][k]);
            }
        }
    }

}

void rowTrans_Ab(double array[MAX][MAX], int n) {
    // 做A、b同时做初等行变换
    int i, j, k;
    double tmp;             //暂存乘积因子
    double zero = 1e-6;       //两浮点数差距小于1e-6视为相等

    for (i = 0; i < n; i++) {
        k = 1;
        while (fabs(array[i][i] - 0.0) < zero && i + k < n) {
            swapArray(array[i], array[i + (k++)], n + 1);
        }

        for (j = i + 1; j < n; j++) {
            if (array[j][i] == 0.0)			//如为0则那一行不用化简 
                continue;
            else {
                tmp = -(double)array[j][i] / array[i][i];		//保存乘积因子 
                for (k = i; k <= n; k++)
                    array[j][k] += (tmp * array[i][k]);
            }
        }
    }

}

void printMatrix(double array[MAX][MAX], int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%lf", array[i][j]);
            if (j == n - 1) {
                printf("\n");
            }
            else {
                printf("\t");
            }
        }
    }
}

void printArray(double x[MAX], int n) {
    int i;
    for (i = 0; i < n; i++) {
        printf("%lf", x[i]);
        if (i == n - 1) {
            printf("\n");
        }
        else {
            printf("\t");
        }
    }
}

void copyMatrix(double a[MAX][MAX], double b[MAX][MAX], int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; i++) {
            b[i][j] = a[i][j];
        }
    }
}

void copyArray(double a[MAX], double b[MAX], int n) {
    int i;
    for (i = 0; i < n; i++) {
        b[i] = a[i];
    }
}

void solveLowerEquations(double A[MAX][MAX], double x[MAX], double b[MAX], int n) {
    //求解系数矩阵为下三角矩阵的线性方程组
    int i, j;
    copyArray(b, x, n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
}

void solveUpperEquations(double A[MAX][MAX], double x[MAX], double b[MAX], int n) {
    //求解系数矩阵为上三角矩阵的线性方程组
    int i, j;
    copyArray(b, x, n);
    for (i = n - 1; i >= 0; i--) {
        for (j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
}