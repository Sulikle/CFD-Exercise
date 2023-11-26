clc
clear all
close all
%% Pre-processing
Unit=8;CFL=0.01;endtau=2;endx=1;deltax=1/Unit;numberx=endx/deltax+1;
% %记录内点位置
% Grid=zeros(1,numberx);
% for i=2:numberx-1
%     Grid(1,i)=(i-1)*deltax+(0.1*rand(1)-0.05)*deltax;
% end
% Grid(1,numberx)=endx;
Uexasolution=zeros(2,numberx);
UDGP0P1plusDGP0=zeros(2,numberx-1);
Unumsolution1=zeros(1,2);
Unumsolution2=zeros(2,numberx-1);
Acc=zeros(3,4);a1=[1/8,1/16,1/32,1/64];a2=[1/8,1/16];
%% solve the question
%solve the numsolution
[UDGP0P1plusDGP0,Grid,endtau01plus0]=BDF1subDGP0P1plusDGP0(Unit,CFL,endtau);%acquire u and ux

 %solve the exasolution
x=0;
for k=1:numberx
    x=Grid(1,k);
    Uexasolution(1,k)=sin(pi*x);
    Uexasolution(2,k)=pi*cos(pi*x);
end
%% post-processing
%plot the u
%plot the DGP0P1plusDGP0
figure
 Unumsolution1(1,1)=UDGP0P1plusDGP0(1,1)-UDGP0P1plusDGP0(2,1)*(Grid(2)-Grid(1))/2;Unumsolution1(1,2)=UDGP0P1plusDGP0(1,1)+UDGP0P1plusDGP0(2,1)*(Grid(2)-Grid(1))/2;
 plot([Grid(1,1),Grid(1,2)],Unumsolution1,'-c.','linewidth',1.5);hold on
 H1=plot([Grid(1,1),Grid(1,2)],Unumsolution1,'-c.','linewidth',1.5);hold on
for i=2:numberx-1
    Unumsolution1(1,1)=UDGP0P1plusDGP0(1,i)-UDGP0P1plusDGP0(2,i)*(Grid(i+1)-Grid(i))/2;Unumsolution1(1,2)=UDGP0P1plusDGP0(1,i)+UDGP0P1plusDGP0(2,i)*(Grid(i+1)-Grid(i))/2;
    plot([Grid(1,i),Grid(1,i+1)],Unumsolution1,'-c.','linewidth',1.5)
end

%plot the exact
y=Grid;
plot(y,Uexasolution(1,:),'-b*','linewidth',1.5)
H2=plot(y,Uexasolution(1,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DGP0P1+DGP0','Exact solution');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
title('Nonuniform grid-BDF1(U)，8 elements','fontsize',16)
hold off


%plot the ux
%DGP0P1plusDGP0
 figure
 Unumsolution1(1,1)=UDGP0P1plusDGP0(2,1);Unumsolution1(1,2)=UDGP0P1plusDGP0(2,1);
 plot([Grid(1,1),Grid(1,2)],Unumsolution1,'-c.','linewidth',1.5);hold on
 H1=plot([Grid(1,1),Grid(1,2)],Unumsolution1,'-c.','linewidth',1.5);hold on
for i=2:numberx-1
    Unumsolution1(1,1)=UDGP0P1plusDGP0(2,i);Unumsolution1(1,2)=UDGP0P1plusDGP0(2,i);
    plot([Grid(1,i),Grid(1,i+1)],Unumsolution1,'-c.','linewidth',1.5)
end

%exact
y=Grid;
plot(y,Uexasolution(2,:),'-b*','linewidth',1.5)
H2=plot(y,Uexasolution(2,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DGP0P1+DGP0','Exact solution');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
title('Nonuniform grid-BDF1(Ux)，8 elements','fontsize',16)
hold off

%determine the accuracy of space U
[U,G,~]=BDF1subDGP0P1plusDGP0(8,CFL,endtau);
Acc(1,1)=AccuracyU(8,U,G);
[U,G,~]=BDF1subDGP0P1plusDGP0(16,CFL,endtau);
Acc(1,2)=AccuracyU(16,U,G);
[U,G,~]=BDF1subDGP0P1plusDGP0(32,CFL,endtau);
Acc(1,3)=AccuracyU(32,U,G);
[U,G,~]=BDF1subDGP0P1plusDGP0(64,CFL,endtau);
Acc(1,4)=AccuracyU(64,U,G);

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
%DGP0P1
[U,G,~]=BDF1subDGP0P1plusDGP0(8,CFL,endtau);
Acc(1,1)=AccuracyUx(8,U,G);
[U,G,~]=BDF1subDGP0P1plusDGP0(16,CFL,endtau);
Acc(1,2)=AccuracyUx(16,U,G);
[U,G,~]=BDF1subDGP0P1plusDGP0(32,CFL,endtau);
Acc(1,3)=AccuracyUx(32,U,G);
[U,G,~]=BDF1subDGP0P1plusDGP0(64,CFL,endtau);
Acc(1,4)=AccuracyUx(64,U,G);

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

