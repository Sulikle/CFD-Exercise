clc
clear all
close all
%% Pre-proceeding
%Some basic paramater
addpath(genpath('D:/CFD������/CFD�����ϰ7/DGP0P1plusDGP0'))
addpath(genpath('D:/CFD������/CFD�����ϰ7/DGP0P2plusDGP1'))
addpath(genpath('D:/CFD������/CFD�����ϰ7/DGP0P2plusrDGP0P1'))
addpath(genpath('D:/CFD������/CFD�����ϰ7/DGP0P3plusrDGP0P2'))
addpath(genpath('D:/CFD������/CFD�����ϰ7/DGP0P3plusrDGP1P2'))
Unit=8;%��Ԫ����
CFL=0.01;
endtau=10;%αʱ����ֵ
tol=10^(-8);%����ѭ������
belta=0.05;%�����Ŷ�ϵ��
%VR�ع�ϵ��
omega0=1;
omega1=0.5;
omega2=0;
omega3=1;
omegab=0;%����ÿ��αʱ�䲽�ϵ�VR�ع�Dirichlet�߽綼Ϊ��ֵ�߽�(Ҳ�������Ϊ�����Ǳ߽�)
norder=0;%����DG��P0Pnorder+1)+DG(Pnorder)
nreconstruct=2;%�ع�����DG��P0,Pnorder+nreconstruct+1)+rDG(Pnorder,Pnorder+nreconstruct)
nexplicit=1;%��ʽ=1������ʽ=0
nsdv=1;%��ʾ��ʹ�õ���ʽ����ʽ����
%% Proceeding
fprintf('1D\n Unit=%d,CFL=%0.3f,belta(�����Ŷ�ϵ��)=%0.3f\n omega0=%0.3f,omega1=%0.3f,omega2=%0.3f,omega3=%0.3f,omegab=%0.3f�������\n\n',Unit,CFL,belta,omega0,omega1,omega2,omega3,omegab);
[n,UL2errors,Vl2errors]=Judges(Unit,CFL,endtau,tol,belta,norder,nreconstruct,nexplicit,nsdv,omega0,omega1,omega2,omega3,omegab);
