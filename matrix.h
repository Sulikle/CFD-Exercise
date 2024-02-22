
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
 //定义结构体来存放矩阵的行，列，以及首地址
typedef struct Matrix
{
    int row;
    int col;
    double **data;
} Matrix;

typedef struct Matrix_TD
{
    int row;
    int col;
    int height;
    double ***data;
} Matrix_TD;

//构建3维矩阵，向堆区申请空间
Matrix_TD MakeMatrix_TD(int row, int col,int height)
{
    int i = 0,j=0;
    Matrix_TD arr1 = {0};
    arr1.row = row;
    arr1.col = col;
    arr1.height=height;
    arr1.data = (double ***)malloc(sizeof(double **) * arr1.row);
    if (arr1.data == NULL)
        exit(-1);
    for (i = 0; i < arr1.row; i++)
    {
        arr1.data[i] = (double **)malloc(sizeof(double*) * arr1.col);
        for (j=0;j<arr1.height;j++)
        {
           arr1.data[i][j] = (double *)malloc(sizeof(double) * arr1.height); 
           memset(arr1.data[i][j], 0, sizeof(double) * arr1.height);
        }
        
    }
    return arr1;
}


//构建矩阵，向堆区申请空间
Matrix MakeMatrix(int row, int col)
{
    int i = 0;
    Matrix arr = {0};
    arr.row = row;
    arr.col = col;
    arr.data = (double **)malloc(sizeof(double *) * arr.row);
    if (arr.data == NULL)
        exit(-1);
    for (i = 0; i < arr.row; i++)
    {
        arr.data[i] = (double *)malloc(sizeof(double) * arr.col);
        memset(arr.data[i], 0, sizeof(double) * arr.col);
    }
    return arr;
}

//释放不必要的矩阵内存，避免堆区内存不够
void free_Matrix(Matrix src)
{
    int i;
    for (i = 0; i < src.row; i++)
    {
        free(src.data[i]);
    }
    free(src.data);
}
//打印矩阵
void print_Matrix(Matrix arr)
{
    int i, j;
    putchar('\n');
    for (i = 0; i < arr.row; i++)
    {
        for (j = 0; j < arr.col; j++)
        {
            printf("%10.4lf ", arr.data[i][j]);
        }
        putchar('\n');
    }
}
//矩阵乘法
Matrix Matrix_Mul(const Matrix left, const Matrix right)
{
    if (left.col != right.row)
        exit(-1); // 判断左列是否等于右行
    int i, j, k;
    Matrix res = MakeMatrix(left.row, right.col);
    for (i = 0; i < left.row; i++)
    {
        for (j = 0; j < right.col; j++)
        {
            for (k = 0; k < left.col; k++)
            {
                res.data[i][j] += left.data[i][k] * right.data[k][j]; 
            }
        }
    }
    return res;
}



//矩阵的LU分解
//U分解
Matrix Matrix_U(Matrix src)
{
    assert(src.row == src.col);
    int i, j, k, principal,nv=0;
    Matrix L1, tmp,P;
    double Lsum, Usum, p, Max;
    L1= MakeMatrix(src.row, src.col);
    tmp = MakeMatrix(src.row, src.col);
    P = MakeMatrix(src.row, src.col);
    for (i = 0; i < src.row; i++)
    {
        memcpy(tmp.data[i], src.data[i], sizeof(src.data[0][0]) * src.col);
    }
    // 初始化斜对角线
    for (i = 0; i < L1.row; i++)
    {
        for (j = 0; j < L1.col; j++)
        {
            if (i == j)
                {L1.data[i][j]=1;
                P.data[i][j]=1;
                } // L1主对角线为1
        }
    }
    
    // 计算U
for(i=0;i<src.col;i++)
{
    // if(tmp.data[i][i]==0)
    // {char name[]="需要左乘P";
    //     printf("%s",name);
    //     nv=1;
    //     for(j=i+1;j<src.row;j++)
    // {
    //     if(fabs(tmp.data[j][i])>0)
    //     {
    //         principal=j;
    //         break;
    //     }
    // }
    // P.data[i][i]=0;P.data[principal][i]=1;P.data[i][principal]=1;P.data[principal][principal]=0;
    // for(k=0;k<src.col;k++)
    // {
    //     double a=tmp.data[i][k];
    //     tmp.data[i][k]=tmp.data[principal][k];
    //     tmp.data[principal][k]=a;
    // }
    // }

    for(j=i+1;j<src.row;j++)
    {
        L1.data[j][i]=-tmp.data[j][i]/tmp.data[i][i];
    }
    Matrix D=Matrix_Mul(L1,tmp);
    tmp=D;
    for(j=i+1;j<src.row;j++)
    {
        L1.data[j][i]=0;
    }

}

   // free_Matrix(tmp);
    if(nv==0)//判断是否发生对换
    {return tmp;//没有则返还U
    }
    else
    {
        return P;//否则返还应左乘的正交矩阵P
    }
}

//L分解
Matrix Matrix_L(Matrix src)
{
    assert(src.row == src.col);
    int i, j, k, principal;
    Matrix L, L1,L2,tmp;
    double Lsum, Usum, p, Max;
    L = MakeMatrix(src.row, src.col);
    L1 = MakeMatrix(src.row, src.col);
    L2 = MakeMatrix(src.row, src.col);
    tmp = MakeMatrix(src.row, src.col);
    for (i = 0; i < src.row; i++)
    {
        memcpy(tmp.data[i], src.data[i], sizeof(src.data[0][0]) * src.col);
    }
    // 初始化斜对角线
    for (i = 0; i < L.row; i++)
    {
        for (j = 0; j < L.col; j++)
        {
            if (i == j)
               {L.data[i][j] = 1;
                L1.data[i][j] = 1; 
                L2.data[i][j] = 1;
                }// L主对角线为1
        }
    }

    // 计算L
    for(j=1;j<src.col;j++)
    {
        L.data[j][0]=tmp.data[j][0]/tmp.data[0][0];
    }


for(i=0;i<src.col;i++)
{
    //判断是否需要列主元对换
    if(tmp.data[i][i]==0)
    {
       char name[]="列主元消去需要左乘P";
        printf("%s",name);
        break; 
    }
    //赋值
    for(j=i+1;j<src.row;j++)
    {
        L1.data[j][i]=-tmp.data[j][i]/tmp.data[i][i];
        L2.data[j][i]=tmp.data[j][i]/tmp.data[i][i];
    }
        if(i>0)
{Matrix C=Matrix_Mul(L,L2);
L=C;

}

    Matrix D=Matrix_Mul(L1,tmp);
    tmp=D;


        for(j=i+1;j<src.row;j++)
    {
        L1.data[j][i]=0;
        L2.data[j][i]=0;
    }


}
    free_Matrix(tmp);

    return L;

}

//下三角矩阵求解方程组
Matrix solveEqu_L(Matrix A,Matrix b)
{
    assert(b.col == 1&&A.row==A.col&&A.data[0][1]==0);
    Matrix X;
    X=MakeMatrix(A.row,1);double sum;int i,j,k;
    for(i=0;i<A.row;i++)
    {
        sum=0;
        for(j=0;j<i;j++)
        {
            sum+=A.data[i][j]*X.data[j][0];
        }
        X.data[i][0]=(b.data[i][0]-sum)/A.data[i][i];
    }
    return X;

}
//上三角矩阵求解方程组
Matrix solveEqu_U(Matrix A,Matrix b)
{
    assert(b.col == 1&&A.row==A.col&&fabs(A.data[A.row-1][0])==0);
    Matrix X;
    X=MakeMatrix(A.row,1);double sum;int i,j,k;
    for(i=A.row-1;i>=0;i--)
    {
        sum=0;
        for(j=A.col-1;j>i;j--)
        {
            sum+=A.data[i][j]*X.data[j][0];
        }
        X.data[i][0]=(b.data[i][0]-sum)/A.data[i][i];
    }
    return X;

}

