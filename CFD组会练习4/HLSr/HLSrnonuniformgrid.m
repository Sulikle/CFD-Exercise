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
Unumsolution1=zeros(1,2);
Unumsolution2=zeros(2,numberx-1);
Acc=zeros(3,4);a1=[1/8,1/16,1/32,1/64];a2=[1/8,1/16];
%% Proceeding
%对f
for k=1:Unit
    Unumsolution(1,k)=(Grid(k+1)-Grid(k)+0.5*(Grid(k+1)^2-Grid(k)^2)+(Grid(k+1)^3-Grid(k)^3)/3)/(Grid(k+1)-Grid(k));
    Unumsolution(2,k)=F(0.5*(Grid(k)+Grid(k+1)))*(Grid(k+1)-Grid(k));%store Uxc*deltax
end

%Reconstruct Uxx
for k=2:numberx-2
    xci=0.5*(Grid(k)+Grid(k+1));%xci
    xci1=0.5*(Grid(k+1)+Grid(k+2));%xci+1
    xci2=0.5*(Grid(k-1)+Grid(k));%xci-1
    A=[(xci1-xci)^2/(2*(Grid(k+1)-Grid(k))^2)+1/24*((Grid(k+2)-Grid(k+1))^2/(Grid(k+1)-Grid(k))^2-1);
        (xci1-xci)/(Grid(k+1)-Grid(k))^2;
        (xci2-xci)^2/(2*(Grid(k+1)-Grid(k))^2)+1/24*((Grid(k)-Grid(k-1))^2/(Grid(k+1)-Grid(k))^2-1);
        (xci2-xci)/(Grid(k+1)-Grid(k))^2];
    b=[Unumsolution(1,k+1)-Unumsolution(1,k)-Unumsolution(2,k)*(xci1-xci)/(Grid(k+1)-Grid(k));
        Unumsolution(2,k+1)/(Grid(k+2)-Grid(k+1))-Unumsolution(2,k)/(Grid(k+1)-Grid(k));
        Unumsolution(1,k-1)-Unumsolution(1,k)-Unumsolution(2,k)*(xci2-xci)/(Grid(k+1)-Grid(k));
        Unumsolution(2,k-1)/(Grid(k)-Grid(k-1))-Unumsolution(2,k)/(Grid(k+1)-Grid(k))];
    Ureconstruct(k)=A\b;
end
%boundary
%Left
    xci=0.5*(Grid(1)+Grid(2));%xci
    xci1=0.5*(Grid(2)+Grid(3));%xci+1
A=[(xci1-xci)^2/(2*(Grid(2)-Grid(1))^2)+1/24*((Grid(3)-Grid(2))^2/(Grid(2)-Grid(1))^2-1);
        (xci1-xci)/(Grid(1+1)-Grid(1))^2];
b=[Unumsolution(1,2)-Unumsolution(1,1)-Unumsolution(2,1)*(xci1-xci)/(Grid(2)-Grid(1));
        Unumsolution(2,1+1)/(Grid(1+2)-Grid(1+1))-Unumsolution(2,1)/(Grid(1+1)-Grid(1))];    
    Ureconstruct(1)=A\b;
%Right    
    xci=0.5*(Grid(numberx-1)+Grid(numberx-1+1));%xci
    xci2=0.5*(Grid(numberx-1-1)+Grid(numberx-1));%xci-1    
 A=[(xci2-xci)^2/(2*(Grid(numberx-1+1)-Grid(numberx-1))^2)+1/24*((Grid(numberx-1)-Grid(numberx-1-1))^2/(Grid(numberx-1+1)-Grid(numberx-1))^2-1);
        (xci2-xci)/(Grid(numberx-1+1)-Grid(numberx-1))^2];
 b=[Unumsolution(1,numberx-1-1)-Unumsolution(1,numberx-1)-Unumsolution(2,numberx-1)*(xci2-xci)/(Grid(numberx-1+1)-Grid(numberx-1));
        Unumsolution(2,numberx-1-1)/(Grid(numberx-1)-Grid(numberx-1-1))-Unumsolution(2,numberx-1)/(Grid(numberx-1+1)-Grid(numberx-1))];   
       Ureconstruct(numberx-1)=A\b;

%Post-proceeding
figure
 k=1;
 x=Grid(k):0.2*(Grid(k+1)-Grid(k)):Grid(k+1);
 p=[Ureconstruct(k)/(2*(Grid(k+1)-Grid(k))^2),Unumsolution(2,k)/(Grid(k+1)-Grid(k))-Ureconstruct(k)*0.5*(Grid(k+1)+Grid(k))/(Grid(k+1)-Grid(k))^2,Unumsolution(1,k)+Ureconstruct(k)*(0.5*(Grid(k+1)+Grid(k)))^2/(2*(Grid(k+1)-Grid(k))^2)-Unumsolution(2,k)*0.5*(Grid(k+1)+Grid(k))/(Grid(k+1)-Grid(k))-Ureconstruct(k)/24];
 y=polyval(p,x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on

 for k=2:numberx-1
     x=Grid(k):0.2*(Grid(k+1)-Grid(k)):Grid(k+1);
 p=[Ureconstruct(k)/(2*(Grid(k+1)-Grid(k))^2),Unumsolution(2,k)/(Grid(k+1)-Grid(k))-Ureconstruct(k)*0.5*(Grid(k+1)+Grid(k))/(Grid(k+1)-Grid(k))^2,Unumsolution(1,k)+Ureconstruct(k)*(0.5*(Grid(k+1)+Grid(k)))^2/(2*(Grid(k+1)-Grid(k))^2)-Unumsolution(2,k)*0.5*(Grid(k+1)+Grid(k))/(Grid(k+1)-Grid(k))-Ureconstruct(k)/24];
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
title('Nonuniform grid reconstruct f=x^2+x+1','fontsize',16)
hold off

%计算精度
Acc(1,1)=Accuracy(8);
Acc(1,2)=Accuracy(16);
Acc(1,3)=Accuracy(32);
Acc(1,4)=Accuracy(64);
for k=1:3
accuracyf(k)=(log10(Acc(1,k+1))-log10(Acc(1,k)))./(log10(a1(1,k+1))-log10(a1(1,k)));
end
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
    Unumsolution(1,k)=(cos(pi*Grid(k))-cos(pi*Grid(k+1)))/(pi*(Grid(k+1)-Grid(k)));
    Unumsolution(2,k)=G(0.5*(Grid(k)+Grid(k+1)))*(Grid(k+1)-Grid(k));%store Uxc*deltax
end

%Reconstruct Uxx
for k=2:numberx-2
    xci=0.5*(Grid(k)+Grid(k+1));%xci
    xci1=0.5*(Grid(k+1)+Grid(k+2));%xci+1
    xci2=0.5*(Grid(k-1)+Grid(k));%xci-1
    A=[(xci1-xci)^2/(2*(Grid(k+1)-Grid(k))^2)+1/24*((Grid(k+2)-Grid(k+1))^2/(Grid(k+1)-Grid(k))^2-1);
        (xci1-xci)/(Grid(k+1)-Grid(k))^2;
        (xci2-xci)^2/(2*(Grid(k+1)-Grid(k))^2)+1/24*((Grid(k)-Grid(k-1))^2/(Grid(k+1)-Grid(k))^2-1);
        (xci2-xci)/(Grid(k+1)-Grid(k))^2];
    b=[Unumsolution(1,k+1)-Unumsolution(1,k)-Unumsolution(2,k)*(xci1-xci)/(Grid(k+1)-Grid(k));
        Unumsolution(2,k+1)/(Grid(k+2)-Grid(k+1))-Unumsolution(2,k)/(Grid(k+1)-Grid(k));
        Unumsolution(1,k-1)-Unumsolution(1,k)-Unumsolution(2,k)*(xci2-xci)/(Grid(k+1)-Grid(k));
        Unumsolution(2,k-1)/(Grid(k)-Grid(k-1))-Unumsolution(2,k)/(Grid(k+1)-Grid(k))];
    Ureconstruct(k)=A\b;
end
%boundary
%Left
    xci=0.5*(Grid(1)+Grid(2));%xci
    xci1=0.5*(Grid(2)+Grid(3));%xci+1
A=[(xci1-xci)^2/(2*(Grid(2)-Grid(1))^2)+1/24*((Grid(3)-Grid(2))^2/(Grid(2)-Grid(1))^2-1);
        (xci1-xci)/(Grid(1+1)-Grid(1))^2];
b=[Unumsolution(1,2)-Unumsolution(1,1)-Unumsolution(2,1)*(xci1-xci)/(Grid(2)-Grid(1));
        Unumsolution(2,1+1)/(Grid(1+2)-Grid(1+1))-Unumsolution(2,1)/(Grid(1+1)-Grid(1))];    
    Ureconstruct(1)=A\b;
%Right    
    xci=0.5*(Grid(numberx-1)+Grid(numberx-1+1));%xci
    xci2=0.5*(Grid(numberx-1-1)+Grid(numberx-1));%xci-1    
 A=[(xci2-xci)^2/(2*(Grid(numberx-1+1)-Grid(numberx-1))^2)+1/24*((Grid(numberx-1)-Grid(numberx-1-1))^2/(Grid(numberx-1+1)-Grid(numberx-1))^2-1);
        (xci2-xci)/(Grid(numberx-1+1)-Grid(numberx-1))^2];
 b=[Unumsolution(1,numberx-1-1)-Unumsolution(1,numberx-1)-Unumsolution(2,numberx-1)*(xci2-xci)/(Grid(numberx-1+1)-Grid(numberx-1));
        Unumsolution(2,numberx-1-1)/(Grid(numberx-1)-Grid(numberx-1-1))-Unumsolution(2,numberx-1)/(Grid(numberx-1+1)-Grid(numberx-1))];   
       Ureconstruct(numberx-1)=A\b;    
%Post-proceeding
figure
 k=1;
 x=Grid(k):0.2*(Grid(k+1)-Grid(k)):Grid(k+1);
 p=[Ureconstruct(k)/(2*(Grid(k+1)-Grid(k))^2),Unumsolution(2,k)/(Grid(k+1)-Grid(k))-Ureconstruct(k)*0.5*(Grid(k+1)+Grid(k))/(Grid(k+1)-Grid(k))^2,Unumsolution(1,k)+Ureconstruct(k)*(0.5*(Grid(k+1)+Grid(k)))^2/(2*(Grid(k+1)-Grid(k))^2)-Unumsolution(2,k)*0.5*(Grid(k+1)+Grid(k))/(Grid(k+1)-Grid(k))-Ureconstruct(k)/24];
 y=polyval(p,x);
 plot(x,y,'-ro','linewidth',1.5);hold on
 H1=plot(x,y,'-ro','linewidth',1.5);hold on

 for k=2:numberx-1
     x=Grid(k):0.2*(Grid(k+1)-Grid(k)):Grid(k+1);
 p=[Ureconstruct(k)/(2*(Grid(k+1)-Grid(k))^2),Unumsolution(2,k)/(Grid(k+1)-Grid(k))-Ureconstruct(k)*0.5*(Grid(k+1)+Grid(k))/(Grid(k+1)-Grid(k))^2,Unumsolution(1,k)+Ureconstruct(k)*(0.5*(Grid(k+1)+Grid(k)))^2/(2*(Grid(k+1)-Grid(k))^2)-Unumsolution(2,k)*0.5*(Grid(k+1)+Grid(k))/(Grid(k+1)-Grid(k))-Ureconstruct(k)/24];
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
title('Nonuniform grid reconstruct g=sin(pi*x) ','fontsize',16)
hold off

%计算精度
Acc(1,1)=Accuracyg(8);
Acc(1,2)=Accuracyg(16);
Acc(1,3)=Accuracyg(32);
Acc(1,4)=Accuracyg(64);
for k=1:3
accuracyg(k)=(log10(Acc(1,k+1))-log10(Acc(1,k)))./(log10(a1(1,k+1))-log10(a1(1,k)));
end
figure
hold on
plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5)
H1=plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5);

H2=plot(log10(a2),1*log10(a2),'--','linewidth',1.5);
plot(log10(a2),2*log10(a2),'--','linewidth',1.5)
H3=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
lgd=legend([H1,H2,H3],'g reconstruct','斜率为1','斜率为2');
lgd.FontSize=12;
xlabel('Log(deltax)','fontsize',14)
ylabel('Log(episilo)','fontsize',14)
title('g精度分析','fontsize',16)










