clc
clear all
close all
%% Pre-proceeding
Unit=8;endx=1;deltax=endx/Unit;numberx=Unit+1;
f=@(x)x.^2+x+1;F=@(x)2*x+1;
g=@(y)sin(pi*y);G=@(y)pi*cos(pi*y);
Unumsolution=zeros(2,Unit);
Ureconstruct=zeros(1,Unit);
Un=zeros(1,numberx);
Unumsolution1=zeros(1,2);
Unumsolution2=zeros(2,numberx-1);
% x=0:deltax:endx;
% plot(x,f(x))

%% Proceeding
%对f
for k=1:Unit
    Unumsolution(1,k)=f((k-0.5)*deltax);
    Unumsolution(2,k)=F((k-0.5)*deltax);
end

%Reconstruct Uxx
for iface=2:numberx-1
    ieL=iface-1;
    ieR=iface;
    Un(iface)=0.5*(Unumsolution(2,ieL)+Unumsolution(2,ieR));
    Ureconstruct(ieL)=Ureconstruct(ieL)+Un(iface)/deltax;
    Ureconstruct(ieR)=Ureconstruct(ieR)-Un(iface)/deltax;
end
    Un(1)=0.5*(Unumsolution(2,1)+Unumsolution(2,1));
    Ureconstruct(1)=Ureconstruct(1)-Un(1)/deltax;
    Un(numberx)=0.5*(Unumsolution(2,numberx-1)+Unumsolution(2,numberx-1));
    Ureconstruct(numberx-1)=Ureconstruct(numberx-1)+Un(numberx)/deltax;
%Post-proceeding
figure
 k=1;
 x=0*deltax:0.2*deltax:1*deltax;
 p=[0.5*Ureconstruct(k),Unumsolution(2,k)-(k-0.5)*deltax*Ureconstruct(k),0.5*Ureconstruct(k)*((k-0.5)*deltax)^2-(k-0.5)*deltax*Unumsolution(2,k)+Unumsolution(1,k)];
 y=polyval(p,x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on

 for k=2:numberx-1
     x=(k-1)*deltax:0.2*deltax:k*deltax;
     p=[0.5*Ureconstruct(k),Unumsolution(2,k)-(k-0.5)*deltax*Ureconstruct(k),0.5*Ureconstruct(k)*((k-0.5)*deltax)^2-(k-0.5)*deltax*Unumsolution(2,k)+Unumsolution(1,k)];
 y=polyval(p,x);
 plot(x,y,'-r^','linewidth',1.5);
 end
hold on
x=0:deltax:endx;
plot(x,f(x),'-b*','linewidth',1.5);
H2=plot(x,f(x),'-b*','linewidth',1.5);
lgd=legend([H1,H2],'Reconstruct','Exact solution');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
title('Reconstruct f=x^2+x+1','fontsize',16)
hold off

%对g
%对f
for k=1:Unit
    Unumsolution(1,k)=g((k-0.5)*deltax);
    Unumsolution(2,k)=G((k-0.5)*deltax);
end

%Reconstruct Uxx
for iface=2:numberx-1
    ieL=iface-1;
    ieR=iface;
    Un(iface)=0.5*(Unumsolution(2,ieL)+Unumsolution(2,ieR));
    Ureconstruct(ieL)=Ureconstruct(ieL)+Un(iface)/deltax;
    Ureconstruct(ieR)=Ureconstruct(ieR)-Un(iface)/deltax;
end
    Un(1)=0.5*(Unumsolution(2,1)+Unumsolution(2,1));
    Ureconstruct(1)=Ureconstruct(1)-Un(1)/deltax;
    Un(numberx)=0.5*(Unumsolution(2,numberx-1)+Unumsolution(2,numberx-1));
    Ureconstruct(numberx-1)=Ureconstruct(numberx-1)+Un(numberx)/deltax;
%Post-proceeding
figure
 k=1;
 x=0*deltax:0.2*deltax:1*deltax;
 p=[0.5*Ureconstruct(k),Unumsolution(2,k)-(k-0.5)*deltax*Ureconstruct(k),0.5*Ureconstruct(k)*((k-0.5)*deltax)^2-(k-0.5)*deltax*Unumsolution(2,k)+Unumsolution(1,k)];
 y=polyval(p,x);
 plot(x,y,'-ro','linewidth',1.5);hold on
 H1=plot(x,y,'-ro','linewidth',1.5);hold on

 for k=2:numberx-1
     x=(k-1)*deltax:0.2*deltax:k*deltax;
     p=[0.5*Ureconstruct(k),Unumsolution(2,k)-(k-0.5)*deltax*Ureconstruct(k),0.5*Ureconstruct(k)*((k-0.5)*deltax)^2-(k-0.5)*deltax*Unumsolution(2,k)+Unumsolution(1,k)];
 y=polyval(p,x);
 plot(x,y,'-ro','linewidth',1.5);
 end
hold on
x=0:0.1*deltax:endx;
plot(x,g(x),'-b','linewidth',1.5);
H2=plot(x,g(x),'-b','linewidth',1.5);
lgd=legend([H1,H2],'Reconstruct','Exact solution');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
title('Reconstruct g=sin(pi*x)','fontsize',16)
hold off