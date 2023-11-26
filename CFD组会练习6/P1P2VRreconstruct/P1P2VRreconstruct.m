clc
clear all
close all
%% Pre-proceeding
Unit=8;endx=1;deltax=endx/Unit;numberx=Unit+1;omega0=1;omega1=1;omega2=1;omegab=omega0/sqrt(2);
%记录内点位置,上下浮动不超过百分之5
Grid=zeros(1,numberx);
Deltax=zeros(1,Unit);
for i=2:numberx-1
    Grid(1,i)=(i-1)*deltax+(0.1*rand(1)-0.05)*deltax;
end
Grid(1,numberx)=endx;
for i=2:numberx
        Deltax(i-1)=Grid(1,i)-Grid(1,i-1);%记录每个单元的区间长度
end
f=@(x)1+x+x.^2;F=@(x)2*x+1;
h=@(x)sin(pi*x);H=@(x)pi*cos(pi*x);
Unumsolution=zeros(2,Unit);
Ureconstruct=zeros(Unit,1);
A=sparse(1:Unit,1:Unit,0,Unit,Unit);
R=zeros(Unit,1);
Unumsolution1=zeros(1,2);
Unumsolution2=zeros(2,numberx-1);
Acc=zeros(3,4);a1=[1/8,1/16,1/32,1/64];a2=[1/8,1/16];

%% Proceeding
%对f
%initial 
for i=1:numberx-1
    Unumsolution(1,i)=(Grid(i+1)-Grid(i)+(Grid(i+1)^2-Grid(i)^2)/2+(Grid(i+1)^3-Grid(i)^3)/3)/(Grid(i+1)-Grid(i));
    Unumsolution(2,i)=(2*0.5*(Grid(i+1)+Grid(i))+1)*(Grid(i+1)-Grid(i));
end

%构建大型稀疏矩阵
for iface=2:numberx-1
    ieL=iface-1;xciL=0.5*(Grid(ieL)+Grid(ieL+1));
    ieR=iface;xciR=0.5*(Grid(ieR)+Grid(ieR+1));
    dLR=0.5*(Deltax(ieR)+Deltax(ieL));
    B2L=(Grid(iface)-xciL)/Deltax(ieL);B2R=(Grid(iface)-xciR)/Deltax(ieR);
    B3L=0.5*B2L^2-1/24;B3R=0.5*B2R^2-1/24;
    
    A(ieL,ieL)=A(ieL,ieL)+2*(omega0^2*B3L^2+omega1^2*dLR^2*B2L^2/Deltax(ieL)^2+omega2^2*dLR^4/Deltax(ieL)^4)/dLR;
    A(ieL,ieL+1)=A(ieL,ieL+1)-2*(omega0^2*B3R*B3L+omega1^2*dLR^2*B2L*B2R/(Deltax(ieL)*Deltax(ieR))+omega2^2*dLR^4/(Deltax(ieL)^2*Deltax(ieR)^2))/dLR;
    R(ieL)=R(ieL)-2*(omega0^2*(Unumsolution(1,ieL)+Unumsolution(2,ieL)*B2L-Unumsolution(1,ieR)-Unumsolution(2,ieR)*B2R)*B3L+omega1^2*B2L*dLR^2*(Unumsolution(2,ieL)/Deltax(ieL)-Unumsolution(2,ieR)/Deltax(ieR))/Deltax(ieL))/dLR;
    
    A(ieR,ieR)=A(ieR,ieR)+2*(omega0^2*B3R^2+omega1^2*dLR^2*B2R^2/Deltax(ieR)^2+omega2^2*dLR^4/Deltax(ieR)^4)/dLR;
    A(ieR,ieR-1)=A(ieR,ieR-1)-2*(omega0^2*B3R*B3L+omega1^2*dLR^2*B2L*B2R/(Deltax(ieL)*Deltax(ieR))+omega2^2*dLR^4/(Deltax(ieL)^2*Deltax(ieR)^2))/dLR;
    R(ieR)=R(ieR)+2*(omega0^2*(Unumsolution(1,ieL)+Unumsolution(2,ieL)*B2L-Unumsolution(1,ieR)-Unumsolution(2,ieR)*B2R)*B3R+omega1^2*B2R*dLR^2*(Unumsolution(2,ieL)/Deltax(ieL)-Unumsolution(2,ieR)/Deltax(ieR))/Deltax(ieR))/dLR;
end

%B.C
%left
iface=1;
ieR=iface;
xciR=0.5*(Grid(ieR)+Grid(ieR+1));
B2R=(Grid(iface)-xciR)/Deltax(ieR);
B3R=0.5*B2R^2-1/24;
A(ieR,ieR)=A(ieR,ieR)+4*omegab^2*B3R^2/Deltax(ieR);
R(ieR)=R(ieR)+4*omegab^2*(1-(Unumsolution(1,ieR)+Unumsolution(2,ieR)*B2R))*B3R/Deltax(ieR);

%Right
iface=numberx;
ieL=iface-1;
xciL=0.5*(Grid(ieL)+Grid(ieL+1));
B2L=(Grid(iface)-xciL)/Deltax(ieL);
B3L=0.5*B2L^2-1/24;
A(ieL,ieL)=A(ieL,ieL)+4*omegab^2*B3L^2/Deltax(ieL);
R(ieL)=R(ieL)-4*omegab^2*(Unumsolution(1,ieL)+Unumsolution(2,ieL)*B2L-3)*B3L/Deltax(ieL);

% Ureconstruct=A\R;

%Thomas 解三对角矩阵
L=zeros(1,Unit);U=zeros(1,Unit);C=zeros(1,Unit);
U(1)=A(1,1);
for i=2:numberx-1
    L(i)=A(i,i-1)/U(i-1);
    U(i)=A(i,i)-L(i)*A(i-1,i);
end
Y=zeros(Unit,1);
Y(1)=R(1);
for i=2:numberx-1
    Y(i)=R(i)-L(i)*Y(i-1);
end
Ureconstruct(numberx-1)=Y(numberx-1)/U(numberx-1);
for i=numberx-2:-1:1
    Ureconstruct(i)=(Y(i)-A(i,i+1)*Ureconstruct(i+1))/U(i);
end

    
    
%% Post-proceeding
figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/Deltax(k)+Ureconstruct(k,1)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on

 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/Deltax(k)+Ureconstruct(k,1)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);
 end
 
hold on
x=Grid(1):0.01*(Grid(numberx)-Grid(1)):Grid(numberx);
plot(x,f(x),'-b','linewidth',1.5);
H2=plot(x,f(x),'-b','linewidth',1.5);
lgd=legend([H1,H2],'Reconstruct','Exact solution');
lgd.FontSize=12;
xlabel('Position x','fontsize',14)
ylabel('Numerical value','fontsize',14)
title('VR(P1P2,NUG) F(x)=1+x+x^2','fontsize',16)
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

H2=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
plot(log10(a2),3*log10(a2),'--','linewidth',1.5)
H3=plot(log10(a2),3*log10(a2),'--','linewidth',1.5);
lgd=legend([H1,H2,H3],'F reconstruct','Slope=2','Slope=3');
lgd.FontSize=12;
xlabel('Log(1/DOF)','fontsize',14)
ylabel('Log(episilo)','fontsize',14)
title('F(x)精度分析 P1P2','fontsize',16)

%% Proceeding
%对h
%initial 
for i=1:numberx-1
    Unumsolution(1,i)=(cos(pi*Grid(i))-cos(pi*Grid(i+1)))/(pi*(Grid(i+1)-Grid(i)));
    Unumsolution(2,i)=pi*cos(pi*0.5*(Grid(i+1)+Grid(i)))*(Grid(i+1)-Grid(i));
end

%构建大型稀疏矩阵
for iface=2:numberx-1
    ieL=iface-1;xciL=0.5*(Grid(ieL)+Grid(ieL+1));
    ieR=iface;xciR=0.5*(Grid(ieR)+Grid(ieR+1));
    dLR=0.5*(Deltax(ieR)+Deltax(ieL));
    B2L=(Grid(iface)-xciL)/Deltax(ieL);B2R=(Grid(iface)-xciR)/Deltax(ieR);
    B3L=0.5*B2L^2-1/24;B3R=0.5*B2R^2-1/24;
    
    A(ieL,ieL)=A(ieL,ieL)+2*(omega0^2*B3L^2+omega1^2*dLR^2*B2L^2/Deltax(ieL)^2+omega2^2*dLR^4/Deltax(ieL)^4)/dLR;
    A(ieL,ieL+1)=A(ieL,ieL+1)-2*(omega0^2*B3R*B3L+omega1^2*dLR^2*B2L*B2R/(Deltax(ieL)*Deltax(ieR))+omega2^2*dLR^4/(Deltax(ieL)^2*Deltax(ieR)^2))/dLR;
    R(ieL)=R(ieL)-2*(omega0^2*(Unumsolution(1,ieL)+Unumsolution(2,ieL)*B2L-Unumsolution(1,ieR)-Unumsolution(2,ieR)*B2R)*B3L+omega1^2*B2L*dLR^2*(Unumsolution(2,ieL)/Deltax(ieL)-Unumsolution(2,ieR)/Deltax(ieR))/Deltax(ieL))/dLR;
    
    A(ieR,ieR)=A(ieR,ieR)+2*(omega0^2*B3R^2+omega1^2*dLR^2*B2R^2/Deltax(ieR)^2+omega2^2*dLR^4/Deltax(ieR)^4)/dLR;
    A(ieR,ieR-1)=A(ieR,ieR-1)-2*(omega0^2*B3R*B3L+omega1^2*dLR^2*B2L*B2R/(Deltax(ieL)*Deltax(ieR))+omega2^2*dLR^4/(Deltax(ieL)^2*Deltax(ieR)^2))/dLR;
    R(ieR)=R(ieR)+2*(omega0^2*(Unumsolution(1,ieL)+Unumsolution(2,ieL)*B2L-Unumsolution(1,ieR)-Unumsolution(2,ieR)*B2R)*B3R+omega1^2*B2R*dLR^2*(Unumsolution(2,ieL)/Deltax(ieL)-Unumsolution(2,ieR)/Deltax(ieR))/Deltax(ieR))/dLR;
end

%B.C
%left
iface=1;
ieR=iface;
xciR=0.5*(Grid(ieR)+Grid(ieR+1));
B2R=(Grid(iface)-xciR)/Deltax(ieR);
B3R=0.5*B2R^2-1/24;
A(ieR,ieR)=A(ieR,ieR)+4*omegab^2*B3R^2/Deltax(ieR);
R(ieR)=R(ieR)+4*omegab^2*(0-(Unumsolution(1,ieR)+Unumsolution(2,ieR)*B2R))*B3R/Deltax(ieR);

%Right
iface=numberx;
ieL=iface-1;
xciL=0.5*(Grid(ieL)+Grid(ieL+1));
B2L=(Grid(iface)-xciL)/Deltax(ieL);
B3L=0.5*B2L^2-1/24;
A(ieL,ieL)=A(ieL,ieL)+4*omegab^2*B3L^2/Deltax(ieL);
R(ieL)=R(ieL)-4*omegab^2*(Unumsolution(1,ieL)+Unumsolution(2,ieL)*B2L-0)*B3L/Deltax(ieL);

%Ureconstruct=A\R;

%Thomas 解三对角矩阵
L=zeros(1,Unit);U=zeros(1,Unit);C=zeros(1,Unit);
U(1)=A(1,1);
for i=2:numberx-1
    L(i)=A(i,i-1)/U(i-1);
    U(i)=A(i,i)-L(i)*A(i-1,i);
end
Y=zeros(Unit,1);
Y(1)=R(1);
for i=2:numberx-1
    Y(i)=R(i)-L(i)*Y(i-1);
end
Ureconstruct(numberx-1)=Y(numberx-1)/U(numberx-1);
for i=numberx-2:-1:1
    Ureconstruct(i)=(Y(i)-A(i,i+1)*Ureconstruct(i+1))/U(i);
end



    
    
%% Post-proceeding
figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/Deltax(k)+Ureconstruct(k)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-rh','linewidth',1.5);hold on
 H1=plot(x,y,'-rh','linewidth',1.5);hold on

 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/Deltax(k)+Ureconstruct(k)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-rh','linewidth',1.5);
 end
 
hold on
x=Grid(1):0.01*(Grid(numberx)-Grid(1)):Grid(numberx);
plot(x,h(x),'-b','linewidth',1.5);
H2=plot(x,h(x),'-b','linewidth',1.5);
lgd=legend([H1,H2],'Reconstruct','Exact solution');
lgd.FontSize=12;
xlabel('Position x','fontsize',14)
ylabel('Numerical value','fontsize',14)
title('VR(P1P2,NUG) G(x)=sin(pi*x)','fontsize',16)
hold off

%计算精度
Acc(1,1)=Accuracyh(8);
Acc(1,2)=Accuracyh(16);
Acc(1,3)=Accuracyh(32);
Acc(1,4)=Accuracyh(64);
for k=1:3
accuracyh(k)=(log10(Acc(1,k+1))-log10(Acc(1,k)))./(log10(a1(1,k+1))-log10(a1(1,k)));
end

figure
hold on
plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5)
H1=plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5);

H2=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
plot(log10(a2),3*log10(a2),'--','linewidth',1.5)
H3=plot(log10(a2),3*log10(a2),'--','linewidth',1.5);
lgd=legend([H1,H2,H3],'G reconstruct','Slope=2','Slope=3');
lgd.FontSize=12;
xlabel('Log(1/DOF)','fontsize',14)
ylabel('Log(episilo)','fontsize',14)
title('G(x)精度分析 P1P2','fontsize',16)












