clc
clear all
close all
%% Pre-processing
Unit=64;CFL=0.01;endtau=2;endx=1;deltax=1/Unit;numberx=endx/deltax+1;
Uexasolution=zeros(2,numberx);
UDGP0P1plusDGP0=zeros(2,numberx-1);
Unumsolution1=zeros(1,2);
Unumsolution2=zeros(2,numberx-1);
Acc=zeros(3,4);a1=[1/8,1/16,1/32,1/64];a2=[1/8,1/16];
%% solve the question
 %solve the exasolution
x=0;
for k=1:numberx
    Uexasolution(1,k)=sin(pi*x);
    Uexasolution(2,k)=pi*cos(pi*x);
    x=x+deltax;
end

[UDGP0P1plusDGP0,endtau01plus0]=TVDsubDGP0P1plusDGP0(Unit,CFL,endtau);%acquire u and ux

%% post-processing
%plot the u
%plot the DGP0plusDGP0
figure
 x=0*deltax:deltax:1*deltax;
 Unumsolution1(1,1)=UDGP0P1plusDGP0(1,1)-UDGP0P1plusDGP0(2,1)*deltax/2;Unumsolution1(1,2)=UDGP0P1plusDGP0(1,1)+UDGP0P1plusDGP0(2,1)*deltax/2;
 plot(x,Unumsolution1,'-c.','linewidth',1.5);hold on
 H1=plot(x,Unumsolution1,'-c.','linewidth',1.5);hold on
for i=2:numberx-1
    x=(i-1)*deltax:deltax:i*deltax;
    Unumsolution1(1,1)=UDGP0P1plusDGP0(1,i)-UDGP0P1plusDGP0(2,i)*deltax/2;Unumsolution1(1,2)=UDGP0P1plusDGP0(1,i)+UDGP0P1plusDGP0(2,i)*deltax/2;
    plot(x,Unumsolution1,'-c.','linewidth',1.5)
end

%plot the exact
y=0:deltax:endx;
plot(y,Uexasolution(1,:),'-b*','linewidth',1.5)
H2=plot(y,Uexasolution(1,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DGP0P1+DGP0','解析解');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
title('面循环TVDRK3数值解与解析解(U)，单元格数为64','fontsize',16)
hold off


%plot the ux
%DGP0plusDGP0
 figure
 x=0*deltax:deltax:1*deltax;
 Unumsolution1(1,1)=UDGP0P1plusDGP0(2,1);Unumsolution1(1,2)=UDGP0P1plusDGP0(2,1);
 plot(x,Unumsolution1,'-c.','linewidth',1.5);hold on
 H1=plot(x,Unumsolution1,'-c.','linewidth',1.5);hold on
for i=2:numberx-1
    x=(i-1)*deltax:deltax:i*deltax;
    Unumsolution1(1,1)=UDGP0P1plusDGP0(2,i);Unumsolution1(1,2)=UDGP0P1plusDGP0(2,i);
    plot(x,Unumsolution1,'-c.','linewidth',1.5)
end

%exact
y=0:deltax:endx;
plot(y,Uexasolution(2,:),'-b*','linewidth',1.5)
H2=plot(y,Uexasolution(2,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DGP0P1+DGP0','解析解');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
title('面循环TVDRK3数值解与解析解(Ux)，单元格数为64','fontsize',16)
hold off

%determine the accuracy of space U
Acc(1,1)=AccuracyU(8,subDGP0P1plusDGP0(8,CFL,endtau));
Acc(1,2)=AccuracyU(16,subDGP0P1plusDGP0(16,CFL,endtau));
Acc(1,3)=AccuracyU(32,subDGP0P1plusDGP0(32,CFL,endtau));
Acc(1,4)=AccuracyU(64,subDGP0P1plusDGP0(64,CFL,endtau));

figure
hold on
plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5)
H1=plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5);

H2=plot(log10(a2),1*log10(a2),'--','linewidth',1.5);
plot(log10(a2),2*log10(a2),'--','linewidth',1.5)
H3=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
lgd=legend([H1,H2,H3],'DGP0P1+DGP0','斜率为1','斜率为2');
lgd.FontSize=12;
xlabel('Log(deltax)','fontsize',14)
ylabel('Log(episilo)','fontsize',14)
title('U空间精度分析','fontsize',16)

%determine the accuracy of space Ux
%DGP0
Acc(1,1)=AccuracyUx(8,subDGP0P1plusDGP0(8,CFL,endtau));
Acc(1,2)=AccuracyUx(16,subDGP0P1plusDGP0(16,CFL,endtau));
Acc(1,3)=AccuracyUx(32,subDGP0P1plusDGP0(32,CFL,endtau));
Acc(1,4)=AccuracyUx(64,subDGP0P1plusDGP0(64,CFL,endtau));

figure
hold on
plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5)
H1=plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5);

plot(log10(a2),1*log10(a2),'--','linewidth',1.5)
H2=plot(log10(a2),1*log10(a2),'--','linewidth',1.5);
plot(log10(a2),2*log10(a2),'--','linewidth',1.5)
H3=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
lgd=legend([H1,H2,H3],'DGP0P1+DGP0','斜率为1','斜率为2');
lgd.FontSize=12;
xlabel('Log(deltax)','fontsize',14)
ylabel('Log(episilo)','fontsize',14)
title('Ux空间精度分析','fontsize',16)

