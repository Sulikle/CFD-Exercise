clc
clear all
close all
%% Pre-proceeding
Unit=8;endx=1;deltax=endx/Unit;numberx=Unit+1;
%记录内点位置,上下浮动不超过百分之5
Grid=zeros(1,numberx);
for i=2:numberx-1
    Grid(1,i)=(i-1)*deltax+(0.1*rand(1)-0.05)*deltax;
end
Grid(1,numberx)=endx;
f=@(x)x.^2+x+1;F=@(x)2*x+1;
g=@(y)sin(pi*y);G=@(y)pi*cos(pi*y);
Unumsolution=zeros(2,Unit);
Ureconstruct=zeros(1,Unit);
Un=zeros(1,numberx);
Unumsolution1=zeros(1,2);
Unumsolution2=zeros(2,numberx-1);
% x=0:deltax:endx;
% plot(x,f(x))
Acc=zeros(3,4);a1=[1/8,1/16,1/32,1/64];a2=[1/8,1/16];
%% Proceeding
%对f
for k=1:Unit
    Unumsolution(1,k)=f(0.5*(Grid(k)+Grid(k+1)));
    Unumsolution(2,k)=F(0.5*(Grid(k)+Grid(k+1)));
end

%Reconstruct Uxx
for iface=2:numberx-1
    ieL=iface-1;
    ieR=iface;
    Un(iface)=((Grid(ieL+1)-Grid(ieL))*Unumsolution(2,ieL)+(Grid(ieR+1)-Grid(ieR))*Unumsolution(2,ieR))/((Grid(ieL+1)-Grid(ieL))+(Grid(ieR+1)-Grid(ieR)));
    Ureconstruct(ieL)=Ureconstruct(ieL)+Un(iface)/deltax;
    Ureconstruct(ieR)=Ureconstruct(ieR)-Un(iface)/deltax;
end
    Un(1)=Unumsolution(2,1);
    Ureconstruct(1)=Ureconstruct(1)-Un(1)/deltax;
    Un(numberx)=0.5*(Unumsolution(2,numberx-1)+Unumsolution(2,numberx-1));
    Ureconstruct(numberx-1)=Ureconstruct(numberx-1)+Un(numberx)/deltax;
%Post-proceeding
figure
 k=1;
 x=Grid(k):0.2*(Grid(k+1)-Grid(k)):Grid(k+1);
 p=[0.5*Ureconstruct(k),Unumsolution(2,k)-0.5*(Grid(k)+Grid(k+1))*Ureconstruct(k),0.5*Ureconstruct(k)*(0.5*(Grid(k)+Grid(k+1)))^2-0.5*(Grid(k)+Grid(k+1))*Unumsolution(2,k)+Unumsolution(1,k)];
 y=polyval(p,x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on

 for k=2:numberx-1
     x=Grid(k):0.2*(Grid(k+1)-Grid(k)):Grid(k+1);
     p=[0.5*Ureconstruct(k),Unumsolution(2,k)-0.5*(Grid(k)+Grid(k+1))*Ureconstruct(k),0.5*Ureconstruct(k)*(0.5*(Grid(k)+Grid(k+1)))^2-0.5*(Grid(k)+Grid(k+1))*Unumsolution(2,k)+Unumsolution(1,k)];
 y=polyval(p,x);
 plot(x,y,'-r^','linewidth',1.5);
 end
hold on
x=Grid(1):0.01*(Grid(numberx)-Grid(1)):Grid(numberx);
plot(x,f(x),'-b','linewidth',1.5);
H2=plot(x,f(x),'-b','linewidth',1.5);
lgd=legend([H1,H2],'Reconstruct','Exact solution');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
title('Ununiform grid reconstruct f=x^2+x+1','fontsize',16)
hold off

%计算精度
Acc(1,1)=Accuracy(8);
Acc(1,2)=Accuracy(16);
Acc(1,3)=Accuracy(32);
Acc(1,4)=Accuracy(64);
figure
hold on
plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5)
H1=plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5);

H2=plot(log10(a2),1*log10(a2),'--','linewidth',1.5);
plot(log10(a2),2*log10(a2),'--','linewidth',1.5)
H3=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
lgd=legend([H1,H2,H3],'f reconstruct','斜率为1','斜率为2');
lgd.FontSize=12;
xlabel('Log(deltax)','fontsize',14)
ylabel('Log(episilo)','fontsize',14)
title('f精度分析','fontsize',16)

%对g
for k=1:Unit
    Unumsolution(1,k)=g(0.5*(Grid(k)+Grid(k+1)));
    Unumsolution(2,k)=G(0.5*(Grid(k)+Grid(k+1)));
end

%Reconstruct Uxx
for iface=2:numberx-1
    ieL=iface-1;
    ieR=iface;
    Un(iface)=((Grid(ieL+1)-Grid(ieL))*Unumsolution(2,ieL)+(Grid(ieR+1)-Grid(ieR))*Unumsolution(2,ieR))/((Grid(ieL+1)-Grid(ieL))+(Grid(ieR+1)-Grid(ieR)));
    Ureconstruct(ieL)=Ureconstruct(ieL)+Un(iface)/deltax;
    Ureconstruct(ieR)=Ureconstruct(ieR)-Un(iface)/deltax;
end
    Un(1)=Unumsolution(2,1);
    Ureconstruct(1)=Ureconstruct(1)-Un(1)/deltax;
    Un(numberx)=0.5*(Unumsolution(2,numberx-1)+Unumsolution(2,numberx-1));
    Ureconstruct(numberx-1)=Ureconstruct(numberx-1)+Un(numberx)/deltax;
%Post-proceeding
figure
 k=1;
 x=Grid(k):0.2*(Grid(k+1)-Grid(k)):Grid(k+1);
 p=[0.5*Ureconstruct(k),Unumsolution(2,k)-0.5*(Grid(k)+Grid(k+1))*Ureconstruct(k),0.5*Ureconstruct(k)*(0.5*(Grid(k)+Grid(k+1)))^2-0.5*(Grid(k)+Grid(k+1))*Unumsolution(2,k)+Unumsolution(1,k)];
 y=polyval(p,x);
 plot(x,y,'-ro','linewidth',1.5);hold on
 H1=plot(x,y,'-ro','linewidth',1.5);hold on

 for k=2:numberx-1
     x=Grid(k):0.2*(Grid(k+1)-Grid(k)):Grid(k+1);
     p=[0.5*Ureconstruct(k),Unumsolution(2,k)-0.5*(Grid(k)+Grid(k+1))*Ureconstruct(k),0.5*Ureconstruct(k)*(0.5*(Grid(k)+Grid(k+1)))^2-0.5*(Grid(k)+Grid(k+1))*Unumsolution(2,k)+Unumsolution(1,k)];
 y=polyval(p,x);
 plot(x,y,'-ro','linewidth',1.5);
 end
hold on
x=Grid(1):0.01*(Grid(numberx)-Grid(1)):Grid(numberx);
plot(x,g(x),'-b','linewidth',1.5);
H2=plot(x,g(x),'-b','linewidth',1.5);
lgd=legend([H1,H2],'Reconstruct','Exact solution');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
title('Ununiform grid reconstruct g=sin(pi*x) ','fontsize',16)
hold off

%计算精度
Acc(1,1)=Accuracyg(8);
Acc(1,2)=Accuracyg(16);
Acc(1,3)=Accuracyg(32);
Acc(1,4)=Accuracyg(64);
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
title('g精度分析','fontsize',16)










