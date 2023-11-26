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
    // ͨ�������б任���ټ���Խ�Ԫ�س˻�������������ʽ��ֵ

    int i, j, k;
    int sign = 0;             //����ʽ����һ����Ҫ�ı���ţ��˱�����¼�������� 
    double sum = 1.0;         //���
    double tmp;             //�ݴ�˻�����
    double zero = 1e-6;       //�����������С��1e-6��Ϊ���

    for (i = 0; i < n - 1; i++) {
        k = 1;
        while (fabs(array[i][i] - 0.0) < zero && i + k < n) {
            swapArray(array[i], array[i + (k++)], n);
            sign++;
        }

        for (j = i + 1; j < n; j++) {
            if (array[j][i] == 0.0)			//��Ϊ0����һ�в��û��� 
                continue;
            else {
                tmp = -(double)array[j][i] / array[i][i];		//����˻����� 
                for (k = i; k < n; k++)
                    array[j][k] += (tmp * array[i][k]);
            }
        }
    }

    for (i = 0; i < n; i++)
        sum *= array[i][i];
    if (sign % 2 != 0)				//����ż���η���Ϊ�� 
        sum *= -1;
    return sum;
}

void rowTrans(double array[MAX][MAX], int n) {
    // �������б任
    int i, j, k;
    double tmp;             //�ݴ�˻�����
    double zero = 1e-6;       //�����������С��1e-6��Ϊ���

    for (i = 0; i < n - 1; i++) {
        k = 1;
        while (fabs(array[i][i] - 0.0) < zero && i + k < n) {
            swapArray(array[i], array[i + (k++)], n);
        }

        for (j = i + 1; j < n; j++) {
            if (array[j][i] == 0.0)			//��Ϊ0����һ�в��û��� 
                continue;
            else {
                tmp = -(double)array[j][i] / array[i][i];		//����˻����� 
                for (k = i; k < n; k++)
                    array[j][k] += (tmp * array[i][k]);
            }
        }
    }

}

void rowTrans_Ab(double array[MAX][MAX], int n) {
    // ��A��bͬʱ�������б任
    int i, j, k;
    double tmp;             //�ݴ�˻�����
    double zero = 1e-6;       //�����������С��1e-6��Ϊ���

    for (i = 0; i < n; i++) {
        k = 1;
        while (fabs(array[i][i] - 0.0) < zero && i + k < n) {
            swapArray(array[i], array[i + (k++)], n + 1);
        }

        for (j = i + 1; j < n; j++) {
            if (array[j][i] == 0.0)			//��Ϊ0����һ�в��û��� 
                continue;
            else {
                tmp = -(double)array[j][i] / array[i][i];		//����˻����� 
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
    //���ϵ������Ϊ�����Ǿ�������Է�����
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
    //���ϵ������Ϊ�����Ǿ�������Է�����
    int i, j;
    copyArray(b, x, n);
    for (i = n - 1; i >= 0; i--) {
        for (j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
}