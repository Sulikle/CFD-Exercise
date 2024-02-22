#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include"accuracy.h"

int main()
{
//Preproceeding
FILE *fp;
char line[255];
int skiprows,Ndim,Nnode,Nelement,Npoint,Nface;
int column1,column2,column3,column4,column5,column6;
double x_coord,y_coord;
double number;
char *token;

// 打开文件
fp = fopen("/Users/mac/Documents/CFD-Exercise-main/CFDProject1/Proj1_Mesh/feflo.domn.cylinder.coarse", "r");

// 检查文件是否成功打开
if (fp == NULL) {
printf("无法打开文件\n");
return 1;
}

fgets(line, sizeof(line), fp);
sscanf(line, "%d", &skiprows);
for(int i=0;i<skiprows+2;i++)
{
   fgets(line, sizeof(line), fp);

}
//读取维数和单元相邻点数
 token = strtok(line, " "); // 假设数据之间以空格分隔
 Ndim = atoi(token); // 使用atof函数将字符串转换为int类型的数据
 token = strtok(NULL, " ");// 继续处理下一个分割出来的字符串
 Nnode=atoi(token);


//读取单元数，点数，边界数
for(int i=0;i<+2;i++)
{
   fgets(line, sizeof(line), fp);

}
 token = strtok(line, " "); // 假设数据之间以空格分隔
 Nelement = atoi(token); // 使用atof函数将字符串转换为int类型的数据
 token = strtok(NULL, " ");// 继续处理下一个分割出来的字符串
 Npoint = atoi(token); // 使用atof函数将字符串转换为int类型的数据
 token = strtok(NULL, " ");// 继续处理下一个分割出来的字符串
 Nface = atoi(token); // 使用atof函数将字符串转换为int类型的数据
//读取INPOEL
Matrix INPOEL=MakeMatrix(Nelement,Nnode),COORD=MakeMatrix(Npoint,Ndim),BCOND=MakeMatrix(Nface,3);
int row=0;
for(int i=0;i<1;i++)
{
   fgets(line, sizeof(line), fp);

}


    while (fscanf(fp , "%d %d %d %d", &column1, &column2,&column3,&column4) == 4) {
     INPOEL.data[row][0]=(double)column2;
     INPOEL.data[row][1]=(double)column3;
     INPOEL.data[row][2]=(double)column4;
     fgets(line, sizeof(line), fp);
  //   printf("%lf\t%lf\t%lf\n",INPOEL.data[row][0],INPOEL.data[row][1],INPOEL.data[row][2]);
          row++;

    }
 //读取COORD  
 row=0;
for(int i=0;i<1;i++)
{
   fgets(line, sizeof(line), fp);

}
    while (fscanf(fp , "%d %lf %lf", &column1, &x_coord,&y_coord) == 3) {
     COORD.data[row][0]=x_coord;
     COORD.data[row][1]=y_coord;
     fgets(line, sizeof(line), fp);
    // printf("%lf\t%lf\n",COORD.data[row][0],COORD.data[row][1]);
          row++;

    }


//读取BCOND
row=0;
    for(int i=0;i<2+Npoint;i++)
{
   fgets(line, sizeof(line), fp);

}

    while (fscanf(fp , "%d %d %d %d", &column1, &column2,&column3,&column4) == 4) {
     BCOND.data[row][0]=(double)column2;
     BCOND.data[row][1]=(double)column3;
     BCOND.data[row][2]=(double)column4;
     fgets(line, sizeof(line), fp);
 //    printf("%lf\t%lf\t%lf\n",BCOND.data[row][0],BCOND.data[row][1],BCOND.data[row][2]);
          row++;

    }
 // 关闭文件
fclose(fp); 
//elements surrounding points
Matrix ESUP2=MakeMatrix(1,Npoint+1);
for(int ie=0;ie<Nelement;ie++)
{
    for(int in=0;in<Nnode;in++)
    {
        ESUP2.data[0][(int)INPOEL.data[ie][in]]+=1;//斜对面存储该point涉及到的单元数
    }

}

for(int IPOIN=1;IPOIN<Npoint+1;IPOIN++)
{
    ESUP2.data[0][IPOIN]+=ESUP2.data[0][IPOIN-1];
}


int MESUP=(int)ESUP2.data[0][Npoint];
Matrix ESUP1=MakeMatrix(1,MESUP);
int IPOIN,ISTOR;
for(int ie=0;ie<Nelement;ie++)
{
    for(int in=0;in<Nnode;in++)
    {
         IPOIN=INPOEL.data[ie][in];//点的标号
         ISTOR=ESUP2.data[0][IPOIN-1]+1;//改点所用单元的初地址
         ESUP2.data[0][IPOIN-1]=ISTOR;//更新地址
         ESUP1.data[0][ISTOR-1]=ie+1;//记录单元

    }
}
for(IPOIN=Npoint;IPOIN>0;IPOIN--)
{
    ESUP2.data[0][IPOIN]=ESUP2.data[0][IPOIN-1];
}
ESUP2.data[0][0]=0;

//stiffness matrix(LHS)&load vecctor(RHS)
Matrix LHS=MakeMatrix(Npoint,Npoint),RHS=MakeMatrix(Npoint,1);
Matrix v_element=MakeMatrix(Nelement,Ndim),v_point=MakeMatrix(Npoint,Ndim),v_point_scalar=MakeMatrix(Npoint,1);
Matrix_TD HP=MakeMatrix_TD(Nelement,Nnode,Nnode);//help matrix
double Bv[1][2]={1,0};

// Proceeding
//construct help matrix(3D)
Matrix HPcycle=MakeMatrix(1,Nnode+2),c_ie=MakeMatrix(1,3);//c per element
double ipoin1,ipoin2,ipoin3;
for (int ie=0;ie<Nelement;ie++)
{
    ipoin1=INPOEL.data[ie][0];
    ipoin2=INPOEL.data[ie][1];
    ipoin3=INPOEL.data[ie][2];
    HPcycle.data[0][0]=ipoin1;
    HPcycle.data[0][1]=ipoin2;
    HPcycle.data[0][2]=ipoin3;
    HPcycle.data[0][3]=ipoin1;
    HPcycle.data[0][4]=ipoin2;
    
    for(int k=0;k<Nnode;k++)//calculate c
    {

        c_ie.data[0][k]=COORD.data[(int)HPcycle.data[0][k+1]-1][0]*COORD.data[(int)HPcycle.data[0][k+2]-1][1]-COORD.data[(int)HPcycle.data[0][k+2]-1][0]*COORD.data[(int)HPcycle.data[0][k+1]-1][1];
    }
    
    double sum_c_ie=0;
    for(int k=0;k<Nnode;k++)
    {
        sum_c_ie+=c_ie.data[0][k];
    }
    for(int i=0;i<Nnode;i++)
    {
        double a_i=COORD.data[(int)HPcycle.data[0][i+1]-1][1]-COORD.data[(int)HPcycle.data[0][i+2]-1][1];
        double b_i=-(COORD.data[(int)HPcycle.data[0][i+1]-1][0]-COORD.data[(int)HPcycle.data[0][i+2]-1][0]);
        for(int j=0;j<Nnode;j++)
        {
        double a_j=COORD.data[(int)HPcycle.data[0][j+1]-1][1]-COORD.data[(int)HPcycle.data[0][j+2]-1][1];
        double b_j=-(COORD.data[(int)HPcycle.data[0][j+1]-1][0]-COORD.data[(int)HPcycle.data[0][j+2]-1][0]); 
        HP.data[ie][i][j]=(a_i*a_j+b_i*b_j)/(2*sum_c_ie);           
        }
    }
    
}



//construct stiffness matrix(LHS)
for(int ie=0;ie<Nelement;ie++)
{
    for(int i=0;i<Nnode;i++)
    {
        int ipoin=(int)INPOEL.data[ie][i];
        for(int j=0;j<Nnode;j++)
        {
            int jpoin=(int)INPOEL.data[ie][j];
            LHS.data[ipoin-1][jpoin-1]+=HP.data[ie][i][j];
        }
    }
}

//construct load vector
for(int iface=0;iface<Nface;iface++)
{
    int ip1=(int)BCOND.data[iface][0],ip2=(int)BCOND.data[iface][1],flag=(int)BCOND.data[iface][2],ip3=-1;//信息提取
    double gama_e=sqrt(pow(COORD.data[ip1-1][0]-COORD.data[ip2-1][0],(double)2)+pow(COORD.data[ip1-1][1]-COORD.data[ip2-1][1],2));//length of face
    int NESP=ESUP2.data[0][ip1]-ESUP2.data[0][ip1-1];//number of element surrounding ip1
    int je=-1;//element of this boundary
    for(int ie=(int)ESUP2.data[0][ip1-1]+1;ie<=(int)ESUP2.data[0][ip1];ie++)
    {
       for(int in=0;in<Nnode;in++)
        {
            if((int)INPOEL.data[(int)ESUP1.data[0][ie-1]-1][in]==ip2)
            {je=(int)ESUP1.data[0][ie-1];
            break;}
        }
        if(je!=-1)
        {
            break;
        }
       
    }

    for(int in=0;in<Nnode;in++)
    {
        if((int)INPOEL.data[je-1][in]!=ip1&&(int)INPOEL.data[je-1][in]!=ip2)
        {
            ip3=(int)INPOEL.data[je-1][in];//remark the third point of this element
        }
    }

    double Bvector1[1][2]={COORD.data[ip2-1][0]-COORD.data[ip1-1][0],COORD.data[ip2-1][1]-COORD.data[ip1-1][1]};
    if(Bv[0][0]==0&&Bv[0][1]==0)
    {
        printf("%s","Please input a vector which is not zero");
    }
    double dot_product=Bv[0][0]*Bvector1[0][0]+Bv[0][1]*Bvector1[0][1];
    double abs_dot_product=fabs(dot_product);//向量点积的模
    double length_v1=sqrt(Bv[0][0]*Bv[0][0]+Bv[0][1]*Bv[0][1]);
    double length_v2=sqrt(Bvector1[0][0]*Bvector1[0][0]+Bvector1[0][1]*Bvector1[0][1]);
    double length_product=length_v1*length_v2;//向量长度

    double Bvector2[1][2]={COORD.data[ip3-1][0]-COORD.data[ip1-1][0],COORD.data[ip3-1][1]-COORD.data[ip1-1][1]};//该单元内部点与边界点构成的向量
    double V_a=COORD.data[ip2-1][0]-COORD.data[ip1-1][0];
    double V_b=COORD.data[ip2-1][1]-COORD.data[ip1-1][1];
    double g;
    double V_product1[1][2]={0,1};
    double V_product2[1][2]={0,-1};

    if(V_b!=0)
    {
        V_product1[0][0]=-V_b/sqrt(pow(V_b,2)+pow(V_a,2));//法线向量
        V_product1[0][1]=V_a/sqrt(pow(V_b,2)+pow(V_a,2));
        V_product2[0][0]=V_b/sqrt(pow(V_b,2)+pow(V_a,2));
        V_product2[0][1]=-V_a/sqrt(pow(V_b,2)+pow(V_a,2));        
    }

    if(abs_dot_product-length_product==0)//判断速度方向和切线方向是否平行
    {
        g=0;
    }
    else if(abs_dot_product-length_product!=0&&Bvector2[0][0]*V_product1[0][0]+Bvector2[0][1]*V_product1[0][1]<0&&flag==4)//判断外法线向量并判断是否是无穷远处
    {
        g=Bv[0][0]*V_product1[0][0]+Bv[0][1]*V_product1[0][1];
    }
    else if(abs_dot_product-length_product!=0&&Bvector2[0][0]*V_product2[0][0]+Bvector2[0][1]*V_product2[0][1]<0&&flag==4)//判断外法线向量并判断是否是无穷远处
    {
        g=Bv[0][0]*V_product2[0][0]+Bv[0][1]*V_product2[0][1];        
    }
    else if(abs_dot_product-length_product!=0&&flag==2)//判断是否是固定面
    {
        g=0;
    } 
    RHS.data[ip1-1][0]+=0.5*g*gama_e;
    RHS.data[ip2-1][0]+=0.5*g*gama_e;

}

//solve the equation(SGS)
int endtimes=(int)pow(10,7);
double tol=pow(10,-16);
Matrix b=MakeMatrix(Npoint,1),varphi=MakeMatrix(Npoint,1),varphi_start=MakeMatrix(Npoint,1);
for (int i=0;i<Npoint;i++)
{
    b.data[i][0]=RHS.data[i][0];
}


for(int times=1;times<endtimes;times++)
{
    //calculate X forward sweep
    for(int i=0;i<Npoint;i++)
    {
        for(int j=0;j<=i-1;j++)
        {
            RHS.data[i][0]-=LHS.data[i][j]*varphi.data[j][0];
        }
        for(int j=i+1;j<Npoint;j++)
        {
            RHS.data[i][0]-=LHS.data[i][j]*varphi.data[j][0];
        }
        varphi.data[i][0]=RHS.data[i][0]/LHS.data[i][i];
    }
    for(int cycle=0;cycle<Npoint;cycle++)
    {
        RHS.data[cycle][0]=b.data[cycle][0];
    }
   
    //backward sweep
    for(int i=Npoint-1;i>-1;i--)
    {
        for(int j=Npoint-1;j>i;j--)
        {
            RHS.data[i][0]-=LHS.data[i][j]*varphi.data[j][0];
        }
        
        for(int k=i-1;k>-1;k--)
        {
            RHS.data[i][0]-=LHS.data[i][k]*varphi.data[k][0];
        }
        
        varphi.data[i][0]=RHS.data[i][0]/LHS.data[i][i];
    }

    for(int cycle=0;cycle<Npoint;cycle++)
    {
        RHS.data[cycle][0]=b.data[cycle][0];
    }
    Matrix ratio_matrix=MakeMatrix(Npoint,1);

    double ratio=0;
    for(int i=0;i<Npoint;i++)
    {
    ratio_matrix.data[i][0]=fabs(varphi.data[i][0]-varphi_start.data[i][0]);
    if(ratio_matrix.data[i][0]>ratio)
    {
        ratio=ratio_matrix.data[i][0];
    }
    }
    for(int i=0;i<Npoint;i++)
    {
        varphi_start.data[i][0]=varphi.data[i][0];
    }

    if(ratio<tol)
    {
        break;
    }
 

}
//solve v per element

for(int ie=0;ie<Nelement;ie++)
{
    ipoin1=INPOEL.data[ie][0];
    ipoin2=INPOEL.data[ie][1];
    ipoin3=INPOEL.data[ie][2];
    HPcycle.data[0][0]=ipoin1;
    HPcycle.data[0][1]=ipoin2;
    HPcycle.data[0][2]=ipoin3;
    HPcycle.data[0][3]=ipoin1;
    HPcycle.data[0][4]=ipoin2;
    double b_i=0,a_i=0;
        for(int k=0;k<Nnode;k++)
    {
        c_ie.data[0][k]=COORD.data[(int)HPcycle.data[0][k+1]-1][0]*COORD.data[(int)HPcycle.data[0][k+2]-1][1]-COORD.data[(int)HPcycle.data[0][k+2]-1][0]*COORD.data[(int)HPcycle.data[0][k+1]-1][1];
    }
    
    double sum_c_ie=0;
    for(int k=0;k<Nnode;k++)
    {
        sum_c_ie+=c_ie.data[0][k];
    }

        for(int i=0;i<Nnode;i++)
    {
        a_i=COORD.data[(int)HPcycle.data[0][i+1]-1][1]-COORD.data[(int)HPcycle.data[0][i+2]-1][1];
        b_i=-(COORD.data[(int)HPcycle.data[0][i+1]-1][0]-COORD.data[(int)HPcycle.data[0][i+2]-1][0]);
        v_element.data[ie][0]+=varphi.data[(int)INPOEL.data[ie][i]-1][0]*a_i/sum_c_ie;
        v_element.data[ie][1]+=varphi.data[(int)INPOEL.data[ie][i]-1][0]*b_i/sum_c_ie;
    }

}
//solve v in point
for(int IPOIN=0;IPOIN<Npoint;IPOIN++)
{
    double sum_weight_x=0,sum_weight_y=0;
    int NESP=(int)(ESUP2.data[0][IPOIN+1]-ESUP2.data[0][IPOIN]);
    Matrix WEIGHT=MakeMatrix(1,NESP);
        for(int ie=(int)ESUP2.data[0][IPOIN]+1;ie<=(int)ESUP2.data[0][IPOIN+1];ie++)
    {
        ipoin1=INPOEL.data[(int)ESUP1.data[0][ie-1]-1][0];
        ipoin2=INPOEL.data[(int)ESUP1.data[0][ie-1]-1][1];
        ipoin3=INPOEL.data[(int)ESUP1.data[0][ie-1]-1][2];
        HPcycle.data[0][0]=ipoin1;
        HPcycle.data[0][1]=ipoin2;
        HPcycle.data[0][2]=ipoin3;
        HPcycle.data[0][3]=ipoin1;
        HPcycle.data[0][4]=ipoin2;
        for(int k=0;k<Nnode;k++)
        {
            WEIGHT.data[0][ie-1-(int)ESUP2.data[0][IPOIN]]+=COORD.data[(int)HPcycle.data[0][k+1]-1][0]*COORD.data[(int)HPcycle.data[0][k+2]-1][1]-COORD.data[(int)HPcycle.data[0][k+2]-1][0]*COORD.data[(int)HPcycle.data[0][k+1]-1][1];
            
        }
       
    }

    double sum_weight=0;
    for(int i=0;i<NESP;i++)
    {
        sum_weight+=WEIGHT.data[0][i];
    }

        for(int ie=(int)ESUP2.data[0][IPOIN]+1;ie<=(int)ESUP2.data[0][IPOIN+1];ie++)  
    {
        sum_weight_x+=v_element.data[(int)ESUP1.data[0][ie-1]-1][0]*WEIGHT.data[0][ie-1-(int)ESUP2.data[0][IPOIN]];
        sum_weight_y+=v_element.data[(int)ESUP1.data[0][ie-1]-1][1]*WEIGHT.data[0][ie-1-(int)ESUP2.data[0][IPOIN]];
    } 
    v_point.data[IPOIN][0]=sum_weight_x/sum_weight;
    v_point.data[IPOIN][1]=sum_weight_y/sum_weight; 
    free_Matrix(WEIGHT);
}
//calculate the size of v in points
for(int IPOIN=0;IPOIN<Npoint;IPOIN++)
{
    v_point_scalar.data[IPOIN][0]=sqrt(pow(v_point.data[IPOIN][0],2)+pow(v_point.data[IPOIN][1],2));
}


//Postproceeding
 //写入文本

    // 打开文本文件以写入模式
    fp = fopen("Cdata.dat", "w");

    if (fp == NULL) {
        printf("无法打开文件。\n");
        return 1;
    }

    // 写入数据到文本文件
    fprintf(fp,"%s","Variables=x,y,u,v,Potential,Velocity\n");
    fprintf(fp,"Zone n=%d,e=%d,f=fepoint,et=triangle\n",Npoint,Nelement);
    for(int IPOIN=0;IPOIN<Npoint;IPOIN++)
    {
      fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",COORD.data[IPOIN][0],COORD.data[IPOIN][1],v_point.data[IPOIN][0],v_point.data[IPOIN][1],varphi.data[IPOIN][0],v_point_scalar.data[IPOIN][0]); 
    }

    for(int ie=0;ie<Nelement;ie++)
    {
        fprintf(fp,"%d\t%d\t%d\n",(int)INPOEL.data[ie][0],(int)INPOEL.data[ie][1],(int)INPOEL.data[ie][2]);
    }


    // 关闭文件
    fclose(fp);

    printf("数据已成功写入到文件。\n");





return 0;
}