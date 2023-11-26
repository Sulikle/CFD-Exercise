clc
clear all
close all
%% Pre-proceeding
%Some basic paramater
addpath(genpath('D:/CFD课题组/CFD组会练习7/DGP0P1plusDGP0'))
addpath(genpath('D:/CFD课题组/CFD组会练习7/DGP0P2plusDGP1'))
addpath(genpath('D:/CFD课题组/CFD组会练习7/DGP0P2plusrDGP0P1'))
addpath(genpath('D:/CFD课题组/CFD组会练习7/DGP0P3plusrDGP0P2'))
addpath(genpath('D:/CFD课题组/CFD组会练习7/DGP0P3plusrDGP1P2'))
Unit=8;%单元个数
CFL=0.01;
endtau=10;%伪时间阈值
tol=10^(-8);%跳出循环条件
belta=0.05;%网格扰动系数
%VR重构系数
omega0=1;
omega1=0.5;
omega2=0;
omega3=1;
omegab=0;%即令每个伪时间步上的VR重构Dirichlet边界都为数值边界(也可以理解为不考虑边界)
norder=0;%阶数DG（P0Pnorder+1)+DG(Pnorder)
nreconstruct=2;%重构阶数DG（P0,Pnorder+nreconstruct+1)+rDG(Pnorder,Pnorder+nreconstruct)
nexplicit=1;%显式=1或者隐式=0
nsdv=1;%表示所使用的显式或隐式方法
%% Proceeding
fprintf('1D\n Unit=%d,CFL=%0.3f,belta(网格扰动系数)=%0.3f\n omega0=%0.3f,omega1=%0.3f,omega2=%0.3f,omega3=%0.3f,omegab=%0.3f的情况下\n\n',Unit,CFL,belta,omega0,omega1,omega2,omega3,omegab);
[n,UL2errors,Vl2errors]=Judges(Unit,CFL,endtau,tol,belta,norder,nreconstruct,nexplicit,nsdv,omega0,omega1,omega2,omega3,omegab);
