clc
clear all
close all
%% Preproceeding
%Some basic paramater
Unit=8;%单元个数
option=1;
afa=0.06;varphi0=50;
Lr=1/(2*pi);Tr=Lr^2/afa;abslambda=sqrt(afa/Tr);A=[abslambda,0;0,abslambda];
CFLt=0.01;CFLtau=0.01;
endt=10;
endtau=10;%伪时间阈值
tol=10^(-5);%跳出循环条件
belta=0.05;%网格扰动系数
endx=1;deltax=endx/Unit;numberx=endx/deltax+1;
deltat=CFLt*deltax/abslambda;
Vcurrent=zeros(2,numberx-1);
Vn=zeros(2,numberx-1);
Vm=zeros(2,numberx-1);
Vm1=zeros(2,numberx-1);
%RHS
R=zeros(2,Unit);
Rd=zeros(2,Unit);
Rb=zeros(2,Unit);
%记录内点位置,上下浮动不超过belta
Grid=zeros(1,numberx);
Deltax=zeros(1,Unit);
for i=2:numberx-1
    Grid(1,i)=(i-1)*deltax+(2*rand(1)-1)*belta*deltax;
end
Grid(1,numberx)=endx;
%记录每个单元的区间长度
for i=2:numberx
        Deltax(i-1)=Grid(1,i)-Grid(1,i-1);
end

%伪时间上的local time stepping
deltatau=zeros(1,numberx-1);%伪时间变量
for i=1:numberx-1
deltatau(i)=CFLtau*(Grid(1,i+1)-Grid(1,i))/abslambda;%伪时间变量
end
%记录每个物理时间步上的伪时间终止时刻
n=zeros(1,floor(endt/deltat));
%VR重构所需量
Ar=zeros(Unit,Unit);
Rr=zeros(Unit,1);
omega0=1;
omega1=0.5;
omega2=0;
omegab=0;

 %solve the exasolution
 Uexasolution=zeros(2,numberx); Vnumsolution=zeros(2,numberx-1);Vreconstruct=zeros(1,numberx-1);
for k=1:numberx
    Uexasolution(1,k)=varphi0*sin(pi*Grid(k))*exp(-afa*pi^2*endt);
    Uexasolution(2,k)=pi*varphi0*cos(pi*Grid(k))*exp(-afa*pi^2*endt);
end

%% Proceeding
%initial  condition set up
for k=1:numberx-1
    Vcurrent(1,k)=-varphi0*(cos(pi*Grid(k+1))-cos(pi*Grid(k)))/(pi*Deltax(k));
    Vcurrent(2,k)=varphi0*(sin(pi*Grid(k+1))-sin(pi*Grid(k)));
%     Vcurrent(3,k)=-varphi0*pi^2*sin(pi*0.5*(Grid(k+1)+Grid(k)))*Deltax(k)^2;
end



%对t=0时刻赋值
Vn=Vcurrent;
in=1;
for itime=0:deltat:endt%物理时间迭代
    Vm=Vn;
     %VR重构varphixxDeltax^2

%构建大型稀疏矩阵
for iface=2:numberx-1
    ieL=iface-1;xciL=0.5*(Grid(ieL)+Grid(ieL+1));
    ieR=iface;xciR=0.5*(Grid(ieR)+Grid(ieR+1));
    dLR=0.5*(Deltax(ieR)+Deltax(ieL));
    B2L=(Grid(iface)-xciL)/Deltax(ieL);B2R=(Grid(iface)-xciR)/Deltax(ieR);
    B3L=0.5*B2L^2-1/24;B3R=0.5*B2R^2-1/24;
    
    Ar(ieL,ieL)=Ar(ieL,ieL)+2*(omega0^2*B3L^2+omega1^2*dLR^2*B2L^2/Deltax(ieL)^2+omega2^2*dLR^4/Deltax(ieL)^4)/dLR;
    Ar(ieL,ieL+1)=Ar(ieL,ieL+1)-2*(omega0^2*B3R*B3L+omega1^2*dLR^2*B2L*B2R/(Deltax(ieL)*Deltax(ieR))+omega2^2*dLR^4/(Deltax(ieL)^2*Deltax(ieR)^2))/dLR;
    Rr(ieL)=Rr(ieL)-2*(omega0^2*(Vm(1,ieL)+Vm(2,ieL)*B2L-Vm(1,ieR)-Vm(2,ieR)*B2R)*B3L+omega1^2*B2L*dLR^2*(Vm(2,ieL)/Deltax(ieL)-Vm(2,ieR)/Deltax(ieR))/Deltax(ieL))/dLR;
    
    Ar(ieR,ieR)=Ar(ieR,ieR)+2*(omega0^2*B3R^2+omega1^2*dLR^2*B2R^2/Deltax(ieR)^2+omega2^2*dLR^4/Deltax(ieR)^4)/dLR;
    Ar(ieR,ieR-1)=Ar(ieR,ieR-1)-2*(omega0^2*B3R*B3L+omega1^2*dLR^2*B2L*B2R/(Deltax(ieL)*Deltax(ieR))+omega2^2*dLR^4/(Deltax(ieL)^2*Deltax(ieR)^2))/dLR;
    Rr(ieR)=Rr(ieR)+2*(omega0^2*(Vm(1,ieL)+Vm(2,ieL)*B2L-Vm(1,ieR)-Vm(2,ieR)*B2R)*B3R+omega1^2*B2R*dLR^2*(Vm(2,ieL)/Deltax(ieL)-Vm(2,ieR)/Deltax(ieR))/Deltax(ieR))/dLR;
end

%B.C
%left
iface=1;
ieR=iface;
xciR=0.5*(Grid(ieR)+Grid(ieR+1));
B2R=(Grid(iface)-xciR)/Deltax(ieR);
B3R=0.5*B2R^2-1/24;
Ar(ieR,ieR)=Ar(ieR,ieR)+4*omegab^2*B3R^2/Deltax(ieR);
Rr(ieR)=Rr(ieR)+4*omegab^2*(0-(Vm(1,ieR)+Vm(2,ieR)*B2R))*B3R/Deltax(ieR);

%Right
iface=numberx;
ieL=iface-1;
xciL=0.5*(Grid(ieL)+Grid(ieL+1));
B2L=(Grid(iface)-xciL)/Deltax(ieL);
B3L=0.5*B2L^2-1/24;
Ar(ieL,ieL)=Ar(ieL,ieL)+4*omegab^2*B3L^2/Deltax(ieL);
Rr(ieL)=Rr(ieL)-4*omegab^2*(Vm(1,ieL)+Vm(2,ieL)*B2L-0)*B3L/Deltax(ieL);

% Ureconstruct=A\R;

%Thomas 解三对角矩阵
L=zeros(1,Unit);U=zeros(1,Unit);C=zeros(1,Unit);
U(1)=Ar(1,1);
for i=2:numberx-1
    L(i)=Ar(i,i-1)/U(i-1);
    U(i)=Ar(i,i)-L(i)*Ar(i-1,i);
end
Y=zeros(Unit,1);
Y(1)=Rr(1);
for i=2:numberx-1
    Y(i)=Rr(i)-L(i)*Y(i-1);
end
Vreconstruct(numberx-1)=Y(numberx-1)/U(numberx-1);
for i=numberx-2:-1:1
    Vreconstruct(i)=(Y(i)-Ar(i,i+1)*Vreconstruct(i+1))/U(i);
end
    
    for itau=0:min(deltatau):endtau%伪时间迭代
        
        %组装RHS
        %利用Gauss积分计算Rdomain
        t=[-sqrt(15)/5,0,sqrt(15)/5];
        W=[5/9,8/9,5/9];
       for ie=1:Unit
           xci=(Grid(ie+1)+Grid(ie))/2;
           for ig=1:3
               xig=Deltax(ie)/2*t(ig)+xci;
               B2g=(xig-xci)/Deltax(ie);
               B3g=0.5*B2g^2-1/24;
               varphig=Vm(1,ie)+B2g*Vm(2,ie)+B3g*Vreconstruct(1,ie);
               Vg=Vm(2,ie)/Deltax(ie)+Vreconstruct(1,ie)*B2g/Deltax(ie);
               S1g=0;
               S2g=-Vg/Tr;
               F1g=-afa*Vg;
               F2g=-varphig/Tr;
               const=W(ig)*0.5*Deltax(ie);
               Rd(2*(ie-1)+1)=Rd(2*(ie-1)+1)+S1g*const;
               Rd(2*(ie-1)+2)=Rd(2*(ie-1)+2)+(S1g*B2g+S2g/Deltax(ie)+F1g/Deltax(ie))*const;
%                Rd(3*(ie-1)+3)=Rd(3*(ie-1)+3)+(S1g*B3g+S2g*B2g/Deltax(ie)+F1g*B2g/Deltax(ie)+F2g/Deltax(ie)^2)*const;
              
           end
       end



        %Rboundary
           for iface=2:numberx-1
               ieL=iface-1;
               ieR=iface;
               B2L=1/2;B3L=0.5*B2L^2-1/24;
               B2R=-1/2;B3R=0.5*B2R^2-1/24;
               varphiL=Vm(1,ieL)+B2L*Vm(2,ieL)+B3L*Vreconstruct(1,ieL);
               varphiR=Vm(1,ieR)+B2R*Vm(2,ieR)+B3R*Vreconstruct(1,ieR);
               VL=Vm(2,ieL)/Deltax(ieL)+Vreconstruct(1,ieL)*B2L/Deltax(ieL);
               VR=Vm(2,ieR)/Deltax(ieR)+Vreconstruct(1,ieR)*B2R/Deltax(ieR);
    
               Fn(:,iface)=0.5*([-afa*VL;-varphiL/Tr]+[-afa*VR;-varphiR/Tr])-0.5*A*([varphiR;VR]-[varphiL;VL]);
               Rb(2*(ieL-1)+1)=Rb(2*(ieL-1)+1)-Fn(1,iface);
               Rb(2*(ieL-1)+2)=Rb(2*(ieL-1)+2)-Fn(1,iface)*B2L-Fn(2,iface)/Deltax(ieL);
%                Rb(3*(ieL-1)+3)=Rb(3*(ieL-1)+3)-Fn(1,iface)*B3L-Fn(2,iface)*B2L/Deltax(ieL);
               Rb(2*(ieR-1)+1)=Rb(2*(ieR-1)+1)+Fn(1,iface);
               Rb(2*(ieR-1)+2)=Rb(2*(ieR-1)+2)+Fn(1,iface)*B2R+Fn(2,iface)/Deltax(ieR);
%                Rb(3*(ieR-1)+3)=Rb(3*(ieR-1)+3)+Fn(1,iface)*B3R+Fn(2,iface)*B2R/Deltax(ieR);
            end
%边界：左
iface=1;
ieR=iface;
B2R=-1/2;B3R=0.5*B2R^2-1/24;
varphiR=Vm(1,ieR)+B2R*Vm(2,ieR)+B3R*Vreconstruct(1,ieR);
VR=Vm(2,ieR)/Deltax(ieR)+Vreconstruct(1,ieR)*B2R/Deltax(ieR);
Fn(:,iface)=0.5*([-afa*VR;0]+[-afa*VR;-varphiR/Tr])-0.5*A*([varphiR;VR]-[0;VR]);
Rb(2*(ieR-1)+1)=Rb(2*(ieR-1)+1)+Fn(1,iface);
Rb(2*(ieR-1)+2)=Rb(2*(ieR-1)+2)+Fn(1,iface)*B2R+Fn(2,iface)/Deltax(ieR);
% Rb(3*(ieR-1)+3)=Rb(3*(ieR-1)+3)+Fn(1,iface)*B3R+Fn(2,iface)*B2R/Deltax(ieR);

%边界：右
iface=numberx;
ieL=iface-1;
B2L=1/2;B3L=0.5*B2L^2-1/24;
varphiL=Vm(1,ieL)+B2L*Vm(2,ieL)+B3L*Vreconstruct(1,ieL);
VL=Vm(2,ieL)/Deltax(ieL)+Vreconstruct(1,ieL)*B2L/Deltax(ieL);
Fn(:,iface)=0.5*([-afa*VL;-varphiL/Tr]+[-afa*VL;0])-0.5*A*([0;VL]-[varphiL;VL]);
Rb(2*(ieL-1)+1)=Rb(2*(ieL-1)+1)-Fn(1,iface);
Rb(2*(ieL-1)+2)=Rb(2*(ieL-1)+2)-Fn(1,iface)*B2L-Fn(2,iface)/Deltax(ieL);
% Rb(3*(ieL-1)+3)=Rb(3*(ieL-1)+3)-Fn(1,iface)*B3L-Fn(2,iface)*B2L/Deltax(ieL);

%R组装
R=Rd+Rb;

if option==1
    for k=1:numberx-1
        Mtau=[Deltax(k),0;0,Deltax(k)/12+1/Deltax(k)];
        Mt=[Deltax(k),0;0,Deltax(k)/12];
        Vm1(:,k)=Vm(:,k)+Mtau\(R(:,k)-Mt*(Vm(:,k)-Vn(:,k))/deltat)*deltatau(k);
    end
end
if option==2
    for k=1:numberx-1
        Mtau=[Deltax(k),0;0,Deltax(k)/12+1/Deltax(k)];
        Mt=[Deltax(k),0;0,Deltax(k)/12];
        Vm1(:,k)=(Mtau/deltatau(k)+Mt/deltat)\(R(:,k)+Mt*Vn(:,k)/deltat+Mtau*Vm(:,k)/deltatau(k));
    end
end


        if max(max(abs(Vm1-Vm)))<tol &&itau>=min(deltatau)
            break;
        end
        
        Vm=Vm1;
        Rd=zeros(2,Unit);
        Rb=zeros(2,Unit);
        Ar=zeros(Unit,Unit);
        Rr=zeros(Unit,1);
        
    end
    n(in)=itau;in=in+1;    
    Vn=Vm1;    
end
Vnumsolution=Vn;

%% Post-proceeding
%U
figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Vnumsolution(1,k)+Vnumsolution(2,k)*(x-xci)/Deltax(k)+Vreconstruct(1,k)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on

 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Vnumsolution(1,k)+Vnumsolution(2,k)*(x-xci)/Deltax(k)+Vreconstruct(1,k)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);
 end
 
%plot the exact
plot(Grid,Uexasolution(1,:),'-b*','linewidth',1.5)
H2=plot(Grid,Uexasolution(1,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DG(P0P2)+rDG(P0P1)','解析解');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
 title('DG(P0P2)+rDG(P0P1) Explicit Euler数值解与解析解(U)','fontsize',16)
hold off

% Ux
 figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Vnumsolution(2,k)/Deltax(k)+Vreconstruct(1,k)*(x-xci)/Deltax(k)^2;
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on
 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
      xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Vnumsolution(2,k)/Deltax(k)+Vreconstruct(1,k)*(x-xci)/Deltax(k)^2;
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);
 end

%exact
plot(Grid,Uexasolution(2,:),'-b*','linewidth',1.5)
H2=plot(Grid,Uexasolution(2,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DG(P0P2)+rDG(P0P1)','解析解');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
title('DG(P0P2)+rDG(P0P1) Explicit Euler数值解与解析解(Ux)','fontsize',16)
hold off






