function A2=Accuracy(Unit)
%% Pre-processing
endx=1;deltax=endx/Unit;numberx=Unit+1;omega0=1;omega1=1;omega2=1;omegab=omega0/sqrt(2);
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
f=@(x)1+x+x.^2;F=@(x)2.*x+1;
Unumsolution=zeros(1,Unit);
Ureconstruct=zeros(2*Unit,1);
A=sparse(1:2*Unit,1:2*Unit,0,2*Unit,2*Unit);
R=zeros(2*Unit,1);


%% Proceeding
%对f
%initial 
for i=1:numberx-1
    Unumsolution(1,i)=(Grid(i+1)-Grid(i)+(Grid(i+1)^2-Grid(i)^2)/2+(Grid(i+1)^3-Grid(i)^3)/3)/(Grid(i+1)-Grid(i));
end

%构建大型分块稀疏矩阵
for iface=2:numberx-1
    ieL=iface-1;xciL=0.5*(Grid(ieL)+Grid(ieL+1));
    ieR=iface;xciR=0.5*(Grid(ieR)+Grid(ieR+1));
    dLR=0.5*(Deltax(ieR)+Deltax(ieL));
    B2L=(Grid(iface)-xciL)/Deltax(ieL);B2R=(Grid(iface)-xciR)/Deltax(ieR);
    B3L=0.5*B2L^2-1/24;B3R=0.5*B2R^2-1/24;
    %diag
    A(2*ieL-1,2*ieL-1)=A(2*ieL-1,2*ieL-1)+2*(omega0^2*B2L^2+omega1^2*dLR^2/Deltax(ieL)^2)/dLR;
    A(2*ieL-1,2*ieL)=A(2*ieL-1,2*ieL)+2*(omega0^2*B3L*B2L+omega1^2*dLR^2*B2L/Deltax(ieL)^2)/dLR;
    A(2*ieL,2*ieL-1)=A(2*ieL,2*ieL-1)+2*(omega0^2*B2L*B3L+omega1^2*dLR^2*B2L/Deltax(ieL)^2)/dLR;
    A(2*ieL,2*ieL)=A(2*ieL,2*ieL)+2*(omega0^2*B3L^2+omega1^2*dLR^2*B2L^2/Deltax(ieL)^2+omega2^2*dLR^4/Deltax(ieL)^4)/dLR;
    A(2*ieR-1,2*ieR-1)=A(2*ieR-1,2*ieR-1)+2*(omega0^2*B2R^2+omega1^2*dLR^2/Deltax(ieR)^2)/dLR;
    A(2*ieR-1,2*ieR)=A(2*ieR-1,2*ieR)+2*(omega0^2*B3R*B2R+omega1^2*dLR^2*B2R/Deltax(ieR)^2)/dLR;
    A(2*ieR,2*ieR-1)=A(2*ieR,2*ieR-1)+2*(omega0^2*B2R*B3R+omega1^2*dLR^2*B2R/Deltax(ieR)^2)/dLR;
    A(2*ieR,2*ieR)=A(2*ieR,2*ieR)+2*(omega0^2*B3R^2+omega1^2*dLR^2*B2R^2/Deltax(ieR)^2+omega2^2*dLR^4/Deltax(ieR)^4)/dLR;
     
    %upper
    A(2*ieL-1,2*ieR-1)=A(2*ieL-1,2*ieR-1)-2*(omega0^2*B2R*B2L+omega1^2*dLR^2/(Deltax(ieL)*Deltax(ieR)))/dLR;
    A(2*ieL-1,2*ieR)=A(2*ieL-1,2*ieR)-2*(omega0^2*B3R*B2L+omega1^2*dLR^2*B2R/(Deltax(ieL)*Deltax(ieR)))/dLR;
    A(2*ieL,2*ieR-1)=A(2*ieL,2*ieR-1)-2*(omega0^2*B2R*B3L+omega1^2*dLR^2*B2L/(Deltax(ieL)*Deltax(ieR)))/dLR;
    A(2*ieL,2*ieR)=A(2*ieL,2*ieR)-2*(omega0^2*B3R*B3L+omega1^2*dLR^2*B2R*B2L/(Deltax(ieL)*Deltax(ieR))+omega2^2*dLR^4/(Deltax(ieL)^2*Deltax(ieR)^2))/dLR;
    %lower
    A(2*ieR-1,2*ieL-1)=A(2*ieR-1,2*ieL-1)-2*(omega0^2*B2L*B2R+omega1^2*dLR^2/(Deltax(ieL)*Deltax(ieR)))/dLR;
    A(2*ieR-1,2*ieL)=A(2*ieR-1,2*ieL)-2*(omega0^2*B3L*B2R+omega1^2*dLR^2*B2L/(Deltax(ieR)*Deltax(ieL)))/dLR;
    A(2*ieR,2*ieL-1)=A(2*ieR,2*ieL-1)-2*(omega0^2*B2L*B3R+omega1^2*dLR^2*B2R/(Deltax(ieR)*Deltax(ieL)))/dLR;
    A(2*ieR,2*ieL)=A(2*ieR,2*ieL)-2*(omega0^2*B3L*B3R+omega1^2*dLR^2*B2L*B2R/(Deltax(ieR)*Deltax(ieL))+omega2^2*dLR^4/(Deltax(ieR)^2*Deltax(ieL)^2))/dLR;
    
    %RHS
    R(2*ieL-1)=R(2*ieL-1)-2*omega0^2*(Unumsolution(ieL)-Unumsolution(ieR))*B2L/dLR;
    R(2*ieL)=R(2*ieL)-2*omega0^2*(Unumsolution(ieL)-Unumsolution(ieR))*B3L/dLR;
    R(2*ieR-1)=R(2*ieR-1)+2*omega0^2*(Unumsolution(ieL)-Unumsolution(ieR))*B2R/dLR;
    R(2*ieR)=R(2*ieR)+2*omega0^2*(Unumsolution(ieL)-Unumsolution(ieR))*B3R/dLR;    
end

%B.C
%left
iface=1;
ieR=iface;
xciR=0.5*(Grid(ieR)+Grid(ieR+1));
B2R=(Grid(iface)-xciR)/Deltax(ieR);
B3R=0.5*B2R^2-1/24;
    A(2*ieR-1,2*ieR-1)=A(2*ieR-1,2*ieR-1)+4*omegab^2*B2R^2/Deltax(ieR);
    A(2*ieR-1,2*ieR)=A(2*ieR-1,2*ieR)+4*omegab^2*B2R*B3R/Deltax(ieR);
    A(2*ieR,2*ieR-1)=A(2*ieR,2*ieR-1)+4*omegab^2*B2R*B3R/Deltax(ieR);
    A(2*ieR,2*ieR)=A(2*ieR,2*ieR)+4*omegab^2*B3R^2/Deltax(ieR);
    
    R(2*ieR-1)=R(2*ieR-1)+4*omegab^2*(1-Unumsolution(ieR))*B2R/Deltax(ieR);
    R(2*ieR)=R(2*ieR)+4*omegab^2*(1-Unumsolution(ieR))*B3R/Deltax(ieR);  

%Right
iface=numberx;
ieL=iface-1;
xciL=0.5*(Grid(ieL)+Grid(ieL+1));
B2L=(Grid(iface)-xciL)/Deltax(ieL);
B3L=0.5*B2L^2-1/24;

    A(2*ieL-1,2*ieL-1)=A(2*ieL-1,2*ieL-1)+4*omegab^2*B2L^2/Deltax(ieL);
    A(2*ieL-1,2*ieL)=A(2*ieL-1,2*ieL)+4*omegab^2*B2L*B3L/Deltax(ieL);
    A(2*ieL,2*ieL-1)=A(2*ieL,2*ieL-1)+4*omegab^2*B2L*B3L/Deltax(ieL);
    A(2*ieL,2*ieL)=A(2*ieL,2*ieL)+4*omegab^2*B3L^2/Deltax(ieL);
    
    R(2*ieL-1)=R(2*ieL-1)-4*omegab^2*(Unumsolution(ieL)-3)*B2L/Deltax(ieL);
    R(2*ieL)=R(2*ieL)-4*omegab^2*(Unumsolution(ieL)-3)*B3L/Deltax(ieL);

%Ureconstruct=A\R;

%LU-SGS 解三对角矩阵
%取出我们所需要的D
D=zeros(2*Unit,2*Unit);
for iface=2:numberx
    ieL=iface-1;
    D(2*ieL-1:2*ieL,2*ieL-1:2*ieL)=A(2*ieL-1:2*ieL,2*ieL-1:2*ieL);
end
%取出我们所需要的L
L=zeros(2*Unit,2*Unit);
for iface=2:numberx-1
    ieR=iface;
    ieL=iface-1;
    L(2*ieR-1:2*ieR,2*ieL-1:2*ieL)=A(2*ieR-1:2*ieR,2*ieL-1:2*ieL);
end

%取出我们所需要的U
U=zeros(2*Unit,2*Unit);
for iface=2:numberx-1
    ieR=iface;
    ieL=iface-1;
    U(2*ieL-1:2*ieL,2*ieR-1:2*ieR)= A(2*ieL-1:2*ieL,2*ieR-1:2*ieR);
end
b=R;%用来存储最初的rhs
Ureconstruct0=zeros(2*Unit,1);
for k=1:30%k表示SGS(k)
    %Forward sweep
    ie=1;
    Ureconstruct(ie:ie+1,1)=D(ie:ie+1,ie:ie+1)\R(ie:ie+1,1);
    for ie=2:numberx-1
        R(2*ie-1:2*ie,1)=R(2*ie-1:2*ie,1)-L(2*ie-1:2*ie,2*(ie-1)-1:2*(ie-1))*Ureconstruct(2*(ie-1)-1:2*(ie-1),1);
        Ureconstruct(2*ie-1:2*ie,1)=D(2*ie-1:2*ie,2*ie-1:2*ie)\R(2*ie-1:2*ie,1);
    end
    %Backward sweep
    for ie=1:numberx-1
        R(2*ie-1:2*ie,1)= D(2*ie-1:2*ie,2*ie-1:2*ie)*Ureconstruct(2*ie-1:2*ie,1);
    end
    
    ie=numberx-1;
    Ureconstruct(2*ie-1:2*ie)=D(2*ie-1:2*ie,2*ie-1:2*ie)\R(2*ie-1:2*ie,1);
    for ie=numberx-2:-1:1
        R(2*ie-1:2*ie,1)=R(2*ie-1:2*ie,1)-U(2*ie-1:2*ie,2*(ie+1)-1:2*(ie+1))*Ureconstruct(2*(ie+1)-1:2*(ie+1),1);
        Ureconstruct(2*ie-1:2*ie,1)=D(2*ie-1:2*ie,2*ie-1:2*ie)\R(2*ie-1:2*ie,1);
    end
    deltaUreconstruct=Ureconstruct;
    Ureconstruct=Ureconstruct0+Ureconstruct;
    if max(deltaUreconstruct)<10^(-10)
        break;
    end

    %重新整理R
    R=b-A*Ureconstruct;
    Ureconstruct0=Ureconstruct;
end  

%calculate the accuracy of space
I1=0;
t=[-sqrt(15)/5,0,sqrt(15)/5];
% t=[-1/sqrt(5),0,1/sqrt(5)];
W=[5/9,8/9,5/9];
k=1;%determine the correctness of the program
for K=1:numberx-1
    xci=(Grid(K+1)+Grid(K))/2;
    deltaxi=Grid(K+1)-Grid(K);
   for i=1:3
       xi=(Grid(K+1)-Grid(K))/2*t(i)+0.5*(Grid(K+1)+Grid(K));
     %p=[Ureconstruct(k)/(2*(Grid(k+1)-Grid(k))^2),Unumsolution(2,k)/(Grid(k+1)-Grid(k))-Ureconstruct(k)*0.5*(Grid(k+1)+Grid(k))/(Grid(k+1)-Grid(k))^2,Unumsolution(1,k)+Ureconstruct(k)*(0.5*(Grid(k+1)+Grid(k)))^2/(2*(Grid(k+1)-Grid(k))^2)-Unumsolution(2,k)*0.5*(Grid(k+1)+Grid(k))/(Grid(k+1)-Grid(k))-Ureconstruct(k)/24];
      %y=polyval(p,xi);
      %fi=(xi^3+xi^2+xi+1-y)^2;
      %对v
       %fi=(3*xi^2+2*xi+1-(Unumsolution(2,K)/deltaxi+Unumsolution(3,K)*(xi-xci)/deltaxi^2+Ureconstruct(K)*((xi-xci)^2/(2*deltaxi^2)-1/24)/deltaxi))^2;
      %对phi
       fi=(xi^2+xi+1-(Unumsolution(1,K)+Ureconstruct(2*K-1,1)*(xi-xci)/deltaxi+Ureconstruct(2*K,1)*(0.5*((xi-xci)/deltaxi)^2-1/24)))^2;
       I1=I1+W(i)*fi*0.5*(Grid(K+1)-Grid(K));        
   end
   k=k+1;
end
A2=sqrt(I1);
end
