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
Unumsolution=zeros(2,Unit);
Ureconstruct=zeros(Unit,1);
A=sparse(1:Unit,1:Unit,0,Unit,Unit);
R=zeros(Unit,1);


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
       fi=(xi^2+xi+1-(Unumsolution(1,K)+Unumsolution(2,K)*(xi-xci)/deltaxi+Ureconstruct(K,1)*(0.5*((xi-xci)/deltaxi)^2-1/24)))^2;
       I1=I1+W(i)*fi*0.5*(Grid(K+1)-Grid(K));        
   end
   k=k+1;
end
A2=sqrt(I1);
end
