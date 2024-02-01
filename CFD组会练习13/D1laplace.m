clc
clear
close all
%% Preproceeding
Unit=32;endx=1;deltax=endx/Unit;npoint=Unit+1;
A=zeros(npoint,npoint);b=zeros(npoint,1);%matrix default setting
unumsolution=zeros(npoint,1);
uxnumsolution=zeros(npoint,1);
uxnumsolution0=zeros(Unit,1);
%计算order
nv=0;
Acc=zeros(3,4);a1=[1/8,1/16,1/32,1/64];a2=[1/32,1/64];
Ue=@(x) 5*x.*(x-2);
Uex=@(x)10*x-10;



%% Proceeding
%set A
for element=1:Unit
    ip1=element;
    ip2=element+1;
    A(ip1,ip1)=A(ip1,ip1)+1/deltax;
    A(ip1,ip2)=A(ip1,ip2)-1/deltax;
    A(ip2,ip1)=A(ip2,ip1)-1/deltax;
    A(ip2,ip2)=A(ip2,ip2)+1/deltax;
end
%Due to initial setting
for ip=2:npoint
    A(1,ip)=0;
    A(ip,1)=0;
end

%set b
b(1,1)=0;%initial setting
for ip=2:npoint-1
    b(ip,1)=-10*deltax;
end
b(npoint,1)=-5*deltax;

%solve the equation
[L,U]=lu(A);
y=L\b;
unumsolution=U\y;

%calculate Ux
for ip=2:npoint-1
    uxnumsolution(ip,1)=0.5*(unumsolution(ip+1,1)/deltax-unumsolution(ip-1,1)/deltax);
end
uxnumsolution(1,1)=unumsolution(2,1)/deltax-unumsolution(1,1)/deltax;
uxnumsolution(npoint,1)=0;

%Ux是常值
for element=1:Unit
    uxnumsolution0(element,1)=unumsolution(element+1,1)/deltax-unumsolution(element,1)/deltax;
end

%% Postproceeding
%U
figure
x=0:deltax:1;
plot(x,unumsolution,'-r^','linewidth',1.5);
H1=plot(x,unumsolution,'-r^','linewidth',1.5);
hold on
plot(x,Ue(x),'-b*','linewidth',1.5);
H2=plot(x,Ue(x),'-b*','linewidth',1.5);
lgd=legend([H1,H2],'FEM数值解','精确解');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
 title('Laplace方程 FEM数值解&精确解(U-8)','fontsize',16)
grid on
 hold off

%Ux
%1次函数的情况
figure
x=0:deltax:1;
plot(x,uxnumsolution,'-r^','linewidth',1.5);
H1=plot(x,uxnumsolution,'-r^','linewidth',1.5);
hold on
plot(x,Uex(x),'-b*','linewidth',1.5);
H2=plot(x,Uex(x),'-b*','linewidth',1.5);
lgd=legend([H1,H2],'FEM数值解','精确解');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
 title('速度-FEM数值解&精确解(U-8)','fontsize',16)
grid on
 hold off
%0次函数的情况

figure
x=0:deltax:deltax;
y=[uxnumsolution0(1,1),uxnumsolution0(1,1)];
plot(x,y,'-r^','linewidth',1.5);
H1=plot(x,y,'-r^','linewidth',1.5);
hold on
for element=2:Unit
x=(element-1)*deltax:deltax:element*deltax;
y=[uxnumsolution0(element,1),uxnumsolution0(element,1)];
plot(x,y,'-r^','linewidth',1.5);
end

x=0:deltax:1;
plot(x,Uex(x),'-b*','linewidth',1.5);
H2=plot(x,Uex(x),'-b*','linewidth',1.5);
lgd=legend([H1,H2],'FEM数值解','精确解');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
 title('速度-FEM数值解&精确解(U-8)','fontsize',16)
grid on
 hold off
 
 
 
 
 
 
 if nv==1

%计算L2误差
[Acc(1,1),Acc(2,1),Acc(3,1)]=Accuracy(8);
[Acc(1,2),Acc(2,2),Acc(3,2)]=Accuracy(16);
[Acc(1,3),Acc(2,3),Acc(3,3)]=Accuracy(32);
[Acc(1,4),Acc(2,4),Acc(3,4)]=Accuracy(64);

%计算order
AccuracyU=zeros(1,3);
AccuracyUx=zeros(1,3);
AccuracyUx0=zeros(1,3);
for k=1:3
AccuracyU(k)=(log10(Acc(1,k+1))-log10(Acc(1,k)))./(log10(a1(1,k+1))-log10(a1(1,k)));
end
for k=1:3
AccuracyUx(k)=(log10(Acc(2,k+1))-log10(Acc(2,k)))./(log10(a1(1,k+1))-log10(a1(1,k)));
end
for k=1:3
AccuracyUx0(k)=(log10(Acc(3,k+1))-log10(Acc(3,k)))./(log10(a1(1,k+1))-log10(a1(1,k)));
end

%U精度
figure
hold on
plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5)
H1=plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5);

H2=plot(log10(a2),1*log10(a2),'--','linewidth',1.5);
 plot(log10(a2),2*log10(a2),'--','linewidth',1.5)
H3=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
lgd=legend([H1,H2,H3],'FEM','Slope=1','Slope=2');
lgd.FontSize=12;
xlabel('Log(1/DOF)','fontsize',14)
ylabel('Log(episilo)','fontsize',14)
title('FEM 精度分析(U)','fontsize',16)
grid on
hold off

%Ux精度
figure
hold on
plot(log10(a1),log10(Acc(2,:)),'-c*','linewidth',1.5)
H1=plot(log10(a1),log10(Acc(2,:)),'-c*','linewidth',1.5);

H2=plot(log10(a2),1*log10(a2),'--','linewidth',1.5);
 plot(log10(a2),2*log10(a2),'--','linewidth',1.5)
H3=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
lgd=legend([H1,H2,H3],'FEM','Slope=1','Slope=2');
lgd.FontSize=12;
xlabel('Log(1/DOF)','fontsize',14)
ylabel('Log(episilo)','fontsize',14)
title('FEM速度精度分析(Ux-1次)','fontsize',16)
grid on
hold off

%Ux精度
figure
hold on
plot(log10(a1),log10(Acc(3,:)),'-c*','linewidth',1.5)
H1=plot(log10(a1),log10(Acc(3,:)),'-c*','linewidth',1.5);

H2=plot(log10(a2),1*log10(a2),'--','linewidth',1.5);
 plot(log10(a2),2*log10(a2),'--','linewidth',1.5)
H3=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
lgd=legend([H1,H2,H3],'FEM','Slope=1','Slope=2');
lgd.FontSize=12;
xlabel('Log(1/DOF)','fontsize',14)
ylabel('Log(episilo)','fontsize',14)
title('FEM速度精度分析(Ux-0次)','fontsize',16)
grid on
hold off

 end
