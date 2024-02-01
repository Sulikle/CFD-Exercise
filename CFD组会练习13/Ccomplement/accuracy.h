#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include"matrix.h"
double** Accuracy(double Unit)
{
Matrix accuracy=MakeMatrix(1,3);
//Preproceeding
const double endx=1;
const double deltax=endx/Unit;
const int npoint=Unit+1;
Matrix A=MakeMatrix(npoint,npoint);
Matrix b=MakeMatrix(npoint,1);
Matrix unumsolution=MakeMatrix(npoint,1);
Matrix uxnumsolution=MakeMatrix(npoint,1);
Matrix uxnumsolution0=MakeMatrix(Unit,1);


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

//calculate the accuracy of spaceU
double I1=0;double t[3]={-sqrt(15.)/5.,0.,sqrt(15.)/5.};double w[3]={5./9,8./9,5./9};
for(int k=0;k<Unit;k++)
{
    double fi,xi;
    for(int i=0;i<3;i++)
    {
        xi=(deltax/2)*t[i]+0.5*(2*(k+1)-1)*deltax;
        fi=pow(5*xi*(xi-2)-unumsolution.data[k][0]-(unumsolution.data[k+1][0]-unumsolution.data[k][0])/deltax*(xi-(k)*deltax),2);
       I1+=w[i]*fi*0.5*(deltax);
    }
}
accuracy.data[0][0]=sqrt(I1);

//calculate the accuracy of spaceUx
I1=0;
for(int k=0;k<Unit;k++)
{
    double fi,xi;
    for(int i=0;i<3;i++)
    {
        xi=deltax/2*t[i]+0.5*(2*(k+1)-1)*deltax;
        fi=pow(10*xi-10-uxnumsolution.data[k][0]-(uxnumsolution.data[k+1][0]-uxnumsolution.data[k][0])/deltax*(xi-(k)*deltax),2);
       I1=I1+w[i]*fi*0.5*(deltax); 
    }
}
accuracy.data[0][1]=sqrt(I1);

//calculate the accuracy of spaceUx0
I1=0;
for(int k=0;k<Unit;k++)
{
    double fi,xi;
    for(int i=0;i<3;i++)
    {
        xi=deltax/2*t[i]+0.5*(2*(k+1)-1)*deltax;
        fi=pow(10*xi-10-uxnumsolution0.data[k][0],2);
       I1=I1+w[i]*fi*0.5*(deltax); 
    }
}
accuracy.data[0][2]=sqrt(I1);
free_Matrix(A);
free_Matrix(b);
free_Matrix(unumsolution);
free_Matrix(uxnumsolution);
free_Matrix(uxnumsolution0);
return accuracy.data;
}
