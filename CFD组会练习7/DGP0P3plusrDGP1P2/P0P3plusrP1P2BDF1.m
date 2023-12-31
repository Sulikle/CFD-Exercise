function [Unumsolution,n]=P0P3plusrP1P2BDF1(Unit,CFL,endtau,tol,Grid,Deltax,nsdv,omega0,omega1,omega2,omega3,omegab)
dimention=3;
%%Pre-proceeding
endx=1;deltax=endx/Unit;numberx=endx/deltax+1;
nu=1;Lr=1/(2*pi);Tr=Lr^2/nu;
abslambda=sqrt(nu/Tr);
deltatau=zeros(1,numberx-1);%伪时间变量
for i=1:numberx-1
deltatau(i)=CFL*(Grid(1,i+1)-Grid(1,i))/abslambda;%伪时间变量
end
B1=1;
%Mtau=[deltax,0;0,deltax/12+1/deltax];dimention=2
%Mtau=[deltax,0,0;0,deltax/12+1/deltax,0;0,0,1/720*deltax+1/(12*deltax)];dimention=3
LHS1=zeros(dimention*Unit,dimention*Unit);
LHS2=zeros(dimention*Unit,dimention*Unit);
LHS3=zeros(dimention*Unit,dimention*Unit);
Rr=zeros(Unit,1);
Ar=zeros(Unit,Unit);
A=[abslambda,0;0,abslambda];
R=zeros(dimention*Unit,1);
Rd=zeros(dimention*Unit,1);
Rb=zeros(dimention*Unit,1);
Fn=zeros(2,numberx);
Ureconstruct=zeros(1,numberx-1);
%为循环所预设的一些量
Ucurrent=zeros(dimention,numberx-1);
Unext=zeros(dimention*Unit,1);

%构建LHS
%Mtau/deltatau
for i=1:Unit
    LHS1(dimention*(i-1)+1,dimention*(i-1)+1)=Deltax(i)/deltatau(i);
    LHS1(dimention*(i-1)+2,dimention*(i-1)+2)=(Deltax(i)/12+1/Deltax(i))/deltatau(i);
    LHS1(dimention*i,dimention*i)=(Deltax(i)/720+1/(12*Deltax(i)))/deltatau(i);
end
%Rdomain
for i=1:Unit
    LHS2(dimention*(i-1)+2,dimention*(i-1)+2)=(1/Tr+nu)/Deltax(i);
    LHS2(dimention*i,dimention*(i-1)+1)=1/(Tr*Deltax(i));
    LHS2(dimention*i,dimention*i)=(1/Tr+nu)/(12*Deltax(i));
end

%Rboundary
for iface=2:numberx-1
    ieL=iface-1;
    ieR=iface;
    CL=[B1,1/2,0.5*(1/2)^2-1/24;0,B1/Deltax(ieL),(1/2)/Deltax(ieL)];
    CR=[B1,-1/2,0.5*(-1/2)^2-1/24;0,B1/Deltax(ieR),(-1/2)/Deltax(ieR)];
    %diag
    LHS3(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieL-1)+1:dimention*ieL)=LHS3(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieL-1)+1:dimention*ieL)+CL'*[abslambda/2,-nu/2;-1/(2*Tr),abslambda/2]*CL;
    LHS3(dimention*(ieR-1)+1:dimention*ieR,dimention*(ieR-1)+1:dimention*ieR)=LHS3(dimention*(ieR-1)+1:dimention*ieR,dimention*(ieR-1)+1:dimention*ieR)-CR'*[-abslambda/2,-nu/2;-1/(2*Tr),-abslambda/2]*CR;
    %upper
    LHS3(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieR-1)+1:dimention*ieR)=LHS3(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieR-1)+1:dimention*ieR)+CL'*[-abslambda/2,-nu/2;-1/(2*Tr),-abslambda/2]*CR;
    %lower
    LHS3(dimention*(ieR-1)+1:dimention*ieR,dimention*(ieL-1)+1:dimention*ieL)=LHS3(dimention*(ieR-1)+1:dimention*ieR,dimention*(ieL-1)+1:dimention*ieL)-CR'*[abslambda/2,-nu/2;-1/(2*Tr),abslambda/2]*CL;
end

%边界：左
iface=1;
ieR=iface;
CR=[B1,-1/2,0.5*(-1/2)^2-1/24;0,B1/Deltax(ieR),(-1/2)/Deltax(ieR)];
LHS3(dimention*(ieR-1)+1:dimention*ieR,dimention*(ieR-1)+1:dimention*ieR)=LHS3(dimention*(ieR-1)+1:dimention*ieR,dimention*(ieR-1)+1:dimention*ieR)-CR'*([abslambda/2,-nu/2;-1/(2*Tr),abslambda/2]*[0,0;0,1]*CR+[-abslambda/2,-nu/2;-1/(2*Tr),-abslambda/2]*CR);
%边界：右
iface=numberx;
ieL=iface-1;
CL=[B1,1/2,0.5*(1/2)^2-1/24;0,B1/Deltax(ieL),(1/2)/Deltax(ieL)];
LHS3(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieL-1)+1:dimention*ieL)=LHS3(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieL-1)+1:dimention*ieL)+CL'*([abslambda/2,-nu/2;-1/(2*Tr),abslambda/2]+[-abslambda/2,-nu/2;-1/(2*Tr),-abslambda/2]*[0,0;0,1])*CL;

%组装LHS
LHS=LHS1+LHS2+LHS3;

% %取出我们所需要的D
% D=zeros(dimention*Unit,dimention*Unit);
% for iface=2:numberx
%     ieL=iface-1;
%     D(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieL-1)+1:dimention*ieL)=LHS(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieL-1)+1:dimention*ieL);
% end


%initial  condition set up

for k=1:numberx-1
    Ucurrent(1,k)=((Grid(k+1)^3-Grid(k)^3)/3-(Grid(k+1)^2-Grid(k)^2)/2)/Deltax(k);
    Ucurrent(2,k)=(Grid(k+1)^2-Grid(k)^2)-(Grid(k+1)-Grid(k));
    Ucurrent(3,k)=2*Deltax(k)^2;
end

%进行必要的向量等价转变
for k=1:numberx-1
    Unext(dimention*(k-1)+1:dimention*k,1)=Ucurrent(:,k);
end


%循环迭代
for n=0:min(deltatau):endtau

%VR重构varphixxxDeltax^3 

%构建大型稀疏矩阵
for iface=2:numberx-1
    ieL=iface-1;xciL=0.5*(Grid(ieL)+Grid(ieL+1));
    ieR=iface;xciR=0.5*(Grid(ieR)+Grid(ieR+1));
    dLR=0.5*(Deltax(ieR)+Deltax(ieL));
    B2L=(Grid(iface)-xciL)/Deltax(ieL);B2R=(Grid(iface)-xciR)/Deltax(ieR);
    B3L=0.5*B2L^2-1/24;B3R=0.5*B2R^2-1/24;
    B4L=(Grid(iface)-xciL)*((Grid(iface)-xciL)^2-Deltax(ieL)^2/4)/(6*Deltax(ieL)^3);
    B4R=(Grid(iface)-xciR)*((Grid(iface)-xciR)^2-Deltax(ieR)^2/4)/(6*Deltax(ieR)^3);
    
    Ar(ieL,ieL)=Ar(ieL,ieL)+2*(omega0^2*B4L^2+omega1^2*dLR^2*B3L^2/Deltax(ieL)^2+omega2^2*B2L^2*dLR^4/Deltax(ieL)^4+omega3^2*dLR^6/Deltax(ieL)^6)/dLR;
    Ar(ieL,ieL+1)=Ar(ieL,ieL+1)-2*(omega0^2*B4R*B4L+omega1^2*dLR^2*B3R*B3L/(Deltax(ieL)*Deltax(ieR))+omega2^2*B2L*B2R*dLR^4/(Deltax(ieL)^2*Deltax(ieR)^2)+omega3^2*dLR^6/(Deltax(ieL)^3*Deltax(ieR)^3))/dLR;
    Rr(ieL)=Rr(ieL)-2*(omega0^2*(Ucurrent(1,ieL)+Ucurrent(2,ieL)*B2L+Ucurrent(3,ieL)*B3L-(Ucurrent(1,ieR)+Ucurrent(2,ieR)*B2R+Ucurrent(3,ieR)*B3R))*B4L+omega1^2*(Ucurrent(2,ieL)/Deltax(ieL)+Ucurrent(3,ieL)*B2L/Deltax(ieL)-Ucurrent(2,ieR)/Deltax(ieR)-Ucurrent(3,ieR)*B2R/Deltax(ieR))*B3L*dLR^2/Deltax(ieL)+omega2^2*B2L*dLR^4*(Ucurrent(3,ieL)/Deltax(ieL)^2-Ucurrent(3,ieR)/Deltax(ieR)^2)/Deltax(ieL)^2)/dLR;
    
    Ar(ieR,ieR)=Ar(ieR,ieR)+2*(omega0^2*B4R^2+omega1^2*dLR^2*B3R^2/Deltax(ieR)^2+omega2^2*B2R^2*dLR^4/Deltax(ieR)^4+omega3^2*dLR^6/Deltax(ieR)^6)/dLR;
    Ar(ieR,ieR-1)=Ar(ieR,ieR-1)-2*(omega0^2*B4L*B4R+omega1^2*dLR^2*B3R*B3L/(Deltax(ieR)*Deltax(ieL))+omega2^2*dLR^4*B2L*B2R/(Deltax(ieL)^2*Deltax(ieR)^2)+omega3^2*dLR^6/(Deltax(ieL)^3*Deltax(ieR)^3))/dLR;
    Rr(ieR)=Rr(ieR)+2*(omega0^2*(Ucurrent(1,ieL)+Ucurrent(2,ieL)*B2L+Ucurrent(3,ieL)*B3L-(Ucurrent(1,ieR)+Ucurrent(2,ieR)*B2R+Ucurrent(3,ieR)*B3R))*B4R+omega1^2*(Ucurrent(2,ieL)/Deltax(ieL)+Ucurrent(3,ieL)*B2L/Deltax(ieL)-Ucurrent(2,ieR)/Deltax(ieR)-Ucurrent(3,ieR)*B2R/Deltax(ieR))*B3R*dLR^2/Deltax(ieR)+omega2^2*B2R*dLR^4*(Ucurrent(3,ieL)/Deltax(ieL)^2-Ucurrent(3,ieR)/Deltax(ieR)^2)/Deltax(ieR)^2)/dLR;
end

%B.C
%left
iface=1;
ieR=iface;
xciR=0.5*(Grid(ieR)+Grid(ieR+1));
B2R=(Grid(iface)-xciR)/Deltax(ieR);
B3R=0.5*B2R^2-1/24;
B4R=(Grid(iface)-xciR)*((Grid(iface)-xciR)^2-Deltax(ieR)^2/4)/(6*Deltax(ieR)^3);

Ar(ieR,ieR)=Ar(ieR,ieR)+4*omegab^2*B4R^2/Deltax(ieR);
Rr(ieR)=Rr(ieR)+4*omegab^2*(0-(Ucurrent(1,ieR)+Ucurrent(2,ieR)*B2R+Ucurrent(3,ieR)*B3R))*B4R/Deltax(ieR);

%Right
iface=numberx;
ieL=iface-1;
xciL=0.5*(Grid(ieL)+Grid(ieL+1));
B2L=(Grid(iface)-xciL)/Deltax(ieL);
B3L=0.5*B2L^2-1/24;
B4L=(Grid(iface)-xciL)*((Grid(iface)-xciL)^2-Deltax(ieL)^2/4)/(6*Deltax(ieL)^3);

Ar(ieL,ieL)=Ar(ieL,ieL)+4*omegab^2*B4L^2/Deltax(ieL);
Rr(ieL)=Rr(ieL)-4*omegab^2*(Ucurrent(1,ieL)+Ucurrent(2,ieL)*B2L+Ucurrent(3,ieL)*B3L-0)*B4L/Deltax(ieL);

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
Ureconstruct(numberx-1)=Y(numberx-1)/U(numberx-1);
for i=numberx-2:-1:1
    Ureconstruct(i)=(Y(i)-Ar(i,i+1)*Ureconstruct(i+1))/U(i);
end

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
        B4g=(xig-xci)*((xig-xci)^2-Deltax(ie)^2/4)/(6*Deltax(ie)^3);
        varphig=Ucurrent(1,ie)+B2g*Ucurrent(2,ie)+B3g*Ucurrent(3,ie)+B4g*Ureconstruct(1,ie);
        Vg=Ucurrent(2,ie)/Deltax(ie)+Ucurrent(3,ie)*B2g/Deltax(ie)+Ureconstruct(1,ie)*B3g/Deltax(ie);
        S1g=pi^2*sin(pi*xig);
        S2g=-Vg/Tr;
        F1g=-nu*Vg;
        F2g=-varphig/Tr;
        const=W(ig)*0.5*Deltax(ie);
        Rd(dimention*(ie-1)+1)=Rd(dimention*(ie-1)+1)+S1g*const;
        Rd(dimention*(ie-1)+2)=Rd(dimention*(ie-1)+2)+(S1g*B2g+S2g/Deltax(ie)+F1g/Deltax(ie))*const;
        Rd(dimention*(ie-1)+3)=Rd(dimention*(ie-1)+3)+(S1g*B3g+S2g*B2g/Deltax(ie)+F1g*B2g/Deltax(ie)+F2g/Deltax(ie)^2)*const;
              
      end
  end



%Rboundary
  for iface=2:numberx-1
    ieL=iface-1;xciL=0.5*(Grid(ieL)+Grid(ieL+1));
    ieR=iface;xciR=0.5*(Grid(ieR)+Grid(ieR+1));
    B2L=1/2;B3L=0.5*B2L^2-1/24;
    B2R=-1/2;B3R=0.5*B2R^2-1/24;
    B4L=(Grid(iface)-xciL)*((Grid(iface)-xciL)^2-Deltax(ieL)^2/4)/(6*Deltax(ieL)^3);
    B4R=(Grid(iface)-xciR)*((Grid(iface)-xciR)^2-Deltax(ieR)^2/4)/(6*Deltax(ieR)^3);
    varphiL=Ucurrent(1,ieL)+B2L*Ucurrent(2,ieL)+B3L*Ucurrent(3,ieL)+B4L*Ureconstruct(1,ieL);
    varphiR=Ucurrent(1,ieR)+B2R*Ucurrent(2,ieR)+B3R*Ucurrent(3,ieR)+B4R*Ureconstruct(1,ieR);
    VL=Ucurrent(2,ieL)/Deltax(ieL)+Ucurrent(3,ieL)*B2L/Deltax(ieL)+Ureconstruct(1,ieL)*B3L/Deltax(ieL);
    VR=Ucurrent(2,ieR)/Deltax(ieR)+Ucurrent(3,ieR)*B2R/Deltax(ieR)+Ureconstruct(1,ieR)*B3R/Deltax(ieR);
    
    Fn(:,iface)=0.5*([-nu*VL;-varphiL/Tr]+[-nu*VR;-varphiR/Tr])-0.5*A*([varphiR;VR]-[varphiL;VL]);
    Rb(dimention*(ieL-1)+1)=Rb(dimention*(ieL-1)+1)-Fn(1,iface);
    Rb(dimention*(ieL-1)+2)=Rb(dimention*(ieL-1)+2)-Fn(1,iface)*B2L-Fn(2,iface)/Deltax(ieL);
    Rb(dimention*(ieL-1)+3)=Rb(dimention*(ieL-1)+3)-Fn(1,iface)*B3L-Fn(2,iface)*B2L/Deltax(ieL);
    Rb(dimention*(ieR-1)+1)=Rb(dimention*(ieR-1)+1)+Fn(1,iface);
    Rb(dimention*(ieR-1)+2)=Rb(dimention*(ieR-1)+2)+Fn(1,iface)*B2R+Fn(2,iface)/Deltax(ieR);
    Rb(dimention*(ieR-1)+3)=Rb(dimention*(ieR-1)+3)+Fn(1,iface)*B3R+Fn(2,iface)*B2R/Deltax(ieR);
  end
%边界：左
iface=1;
ieR=iface;xciR=0.5*(Grid(ieR)+Grid(ieR+1));
B2R=-1/2;B3R=0.5*B2R^2-1/24;
B4R=(Grid(iface)-xciR)*((Grid(iface)-xciR)^2-Deltax(ieR)^2/4)/(6*Deltax(ieR)^3);
varphiR=Ucurrent(1,ieR)+B2R*Ucurrent(2,ieR)+B3R*Ucurrent(3,ieR)+B4R*Ureconstruct(1,ieR);
VR=Ucurrent(2,ieR)/Deltax(ieR)+Ucurrent(3,ieR)*B2R/Deltax(ieR)+Ureconstruct(1,ieR)*B3R/Deltax(ieR);
Fn(:,iface)=0.5*([-nu*VR;0]+[-nu*VR;-varphiR/Tr])-0.5*A*([varphiR;VR]-[0;VR]);
Rb(dimention*(ieR-1)+1)=Rb(dimention*(ieR-1)+1)+Fn(1,iface);
Rb(dimention*(ieR-1)+2)=Rb(dimention*(ieR-1)+2)+Fn(1,iface)*B2R+Fn(2,iface)/Deltax(ieR);
Rb(dimention*(ieR-1)+3)=Rb(dimention*(ieR-1)+3)+Fn(1,iface)*B3R+Fn(2,iface)*B2R/Deltax(ieR);

%边界：右
iface=numberx;
 ieL=iface-1;xciL=0.5*(Grid(ieL)+Grid(ieL+1));
B2L=1/2;B3L=0.5*B2L^2-1/24;B4L=(Grid(iface)-xciL)*((Grid(iface)-xciL)^2-Deltax(ieL)^2/4)/(6*Deltax(ieL)^3);
varphiL=Ucurrent(1,ieL)+B2L*Ucurrent(2,ieL)+B3L*Ucurrent(3,ieL)+B4L*Ureconstruct(1,ieL);
VL=Ucurrent(2,ieL)/Deltax(ieL)+Ucurrent(3,ieL)*B2L/Deltax(ieL)+Ureconstruct(1,ieL)*B3L/Deltax(ieL);
Fn(:,iface)=0.5*([-nu*VL;-varphiL/Tr]+[-nu*VL;0])-0.5*A*([0;VL]-[varphiL;VL]);
Rb(dimention*(ieL-1)+1)=Rb(dimention*(ieL-1)+1)-Fn(1,iface);
Rb(dimention*(ieL-1)+2)=Rb(dimention*(ieL-1)+2)-Fn(1,iface)*B2L-Fn(2,iface)/Deltax(ieL);
 Rb(dimention*(ieL-1)+3)=Rb(dimention*(ieL-1)+3)-Fn(1,iface)*B3L-Fn(2,iface)*B2L/Deltax(ieL);

%R组装
R=Rd+Rb;

%求解
% X=D\R;
if nsdv==1
X=Jacobi(LHS,R,Unit);
elseif nsdv==2
    X=LUSGS(LHS,R,Unit);
end
    if max(X)<tol
        break
    end
    Unext=Unext+X;
Rd=zeros(dimention*Unit,1);
Rb=zeros(dimention*Unit,1);
Rr=zeros(Unit,1);
Ar=zeros(Unit,Unit);

for k=1:numberx-1
    Ucurrent(:,k)=Unext(dimention*(k-1)+1:dimention*k,1);
end


end
Unumsolution(1,:)=Ucurrent(1,:);Unumsolution(2,:)=Ucurrent(2,:);Unumsolution(3,:)=Ucurrent(3,:);Unumsolution(4,:)=Ureconstruct(1,:);

if nsdv==1
fprintf('DG(P0P3)+rDG(P1P2) Jacobi 达到稳态所需要的时间是：%f 秒\n',n)
elseif nsdv==2
    fprintf('DG(P0P3)+rDG(P1P2) LUSGS 达到稳态所需要的时间是：%f 秒\n',n)
end
end
