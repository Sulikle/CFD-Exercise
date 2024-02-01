#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
//#include"matrix.h"
#include"accuracy.h"

int main()
{
//Preproceeding
const double Unit=32;
const double endx=1;
const double deltax=endx/Unit;
const int npoint=Unit+1;int nv=1;
Matrix A=MakeMatrix(npoint,npoint);
Matrix b=MakeMatrix(npoint,1);
Matrix unumsolution=MakeMatrix(npoint,1);
Matrix uxnumsolution=MakeMatrix(npoint,1);
Matrix uxnumsolution0=MakeMatrix(Unit,1);
//计算order
Matrix Acc=MakeMatrix(3,4);
Matrix a1=MakeMatrix(1,4);Matrix a2=MakeMatrix(1,2);
a1.data[0][0]=1/8;a1.data[0][1]=1/16;a1.data[0][2]=1/32;a1.data[0][3]=1/64;
a2.data[0][0]=1/32;a2.data[0][1]=1/64;

//Proceeding
//set A
for(int element=0;element<Unit;element++)
{
    int ip1=element,ip2=element+1;
    A.data[ip1][ip1]=A.data[ip1][ip1]+1/deltax;
    A.data[ip1][ip2]=A.data[ip1][ip2]-1/deltax;
    A.data[ip2][ip1]=A.data[ip2][ip1]-1/deltax;
    A.data[ip2][ip2]=A.data[ip2][ip2]+1/deltax;   
}
//Due to initial setting
for (int ip=1;ip<npoint;ip++)
{
    A.data[0][ip]=0;
    A.data[ip][0]=0;
}
//set b
b.data[0][0]=0;
for(int ip=1;ip<npoint-1;ip++)
{
    b.data[ip][0]=-10*deltax;
}
b.data[npoint-1][0]=-5*deltax;

//solve the equation
Matrix u=Matrix_U(A),l=Matrix_L(A);
Matrix y=MakeMatrix(A.row,1);
y=solveEqu_L(l,b);
unumsolution=solveEqu_U(u,y);

//calculate Ux
for(int ip=1;ip<npoint-1;ip++)
{
    uxnumsolution.data[ip][0]=0.5*(unumsolution.data[ip+1][0]/deltax-unumsolution.data[ip-1][0]/deltax);   
}
uxnumsolution.data[0][0]=unumsolution.data[1][0]/deltax-unumsolution.data[0][0]/deltax;
uxnumsolution.data[npoint-1][0]=0;
for(int element=0;element<Unit;element++)
{
    uxnumsolution0.data[element][0]=unumsolution.data[element+1][0]/deltax-unumsolution.data[element][0]/deltax;
}


if(nv==1)
{
    Acc.data[0][0]=Accuracy(8)[0][0];Acc.data[1][0]=Accuracy(8)[0][1];Acc.data[2][0]=Accuracy(8)[0][2];
    Acc.data[0][1]=Accuracy(16)[0][0];Acc.data[1][1]=Accuracy(16)[0][1];Acc.data[2][1]=Accuracy(16)[0][2];
    Acc.data[0][2]=Accuracy(32)[0][0];Acc.data[1][2]=Accuracy(32)[0][1];Acc.data[2][2]=Accuracy(32)[0][2];
    Acc.data[0][3]=Accuracy(64)[0][0];Acc.data[1][3]=Accuracy(64)[0][1];Acc.data[2][3]=Accuracy(64)[0][2];

}
print_Matrix(Acc);



return 0;
}