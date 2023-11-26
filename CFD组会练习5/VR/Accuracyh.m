function A2=Accuracyh(Unit)
%% Pre-processing
endx=1;deltax=endx/Unit;numberx=Unit+1;omega0=0.5;omega1=0.5;omegab=0;
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
f=@(x)1+x;F=@(x)1;
h=@(x)sin(pi*x);H=@(x)pi*cos(pi*x);
Unumsolution=zeros(1,Unit);
Ureconstruct=zeros(Unit,1);
A=sparse(1:Unit,1:Unit,0,Unit,Unit);
R=zeros(Unit,1);
Unumsolution1=zeros(1,2);
Unumsolution2=zeros(2,numberx-1);
Acc=zeros(3,4);a1=[1/(2*8),1/(2*16),1/(2*32),1/(2*64)];a2=[1/(2*8),1/(2*16)];

%% Proceeding
%对h
%initial 
for i=1:numberx-1
    Unumsolution(1,i)=(cos(pi*Grid(i))-cos(pi*Grid(i+1)))/(pi*(Grid(i+1)-Grid(i)));
end

%构建大型稀疏矩阵
for iface=2:numberx-1
    ieL=iface-1;xciL=0.5*(Grid(ieL)+Grid(ieL+1));
    ieR=iface;xciR=0.5*(Grid(ieR)+Grid(ieR+1));
    dLR=0.5*(Deltax(ieR)+Deltax(ieL));
    A(ieL,ieL)=A(ieL,ieL)+2*(omega0^2*((Grid(iface)-xciL)/Deltax(ieL))^2+omega1^2*dLR^2/Deltax(ieL)^2)/dLR;
    A(ieL,ieL+1)=A(ieL,ieL+1)-2*(omega0^2*((Grid(iface)-xciL)/Deltax(ieL))*((Grid(iface)-xciR)/Deltax(ieR))+omega1^2*dLR^2/(Deltax(ieR)*Deltax(ieL)))/dLR;
    R(ieL)=R(ieL)-2*omega0^2*(Unumsolution(ieL)-Unumsolution(ieR))*((Grid(iface)-xciL)/Deltax(ieL))/dLR;
    
    A(ieR,ieR)=A(ieR,ieR)+2*(omega0^2*((Grid(iface)-xciR)/Deltax(ieR))^2+omega1^2*dLR^2/Deltax(ieR)^2)/dLR;
    A(ieR,ieR-1)=A(ieR,ieR-1)-2*(omega0^2*((Grid(iface)-xciL)/Deltax(ieL))*((Grid(iface)-xciR)/Deltax(ieR))+omega1^2*dLR^2/(Deltax(ieR)*Deltax(ieL)))/dLR;
    R(ieR)=R(ieR)+2*omega0^2*(Unumsolution(ieL)-Unumsolution(ieR))*((Grid(iface)-xciR)/Deltax(ieR))/dLR;
end

%B.C
%left
iface=1;
ieR=iface;
xciR=0.5*(Grid(ieR)+Grid(ieR+1));
A(ieR,ieR)=A(ieR,ieR)+4*omegab^2*((Grid(iface)-xciR)/Deltax(ieR))^2/Deltax(ieR);
R(ieR)=R(ieR)+4*omegab^2*(0-Unumsolution(ieR))*((Grid(iface)-xciR)/Deltax(ieR))/Deltax(ieR);

%Right
iface=numberx;
ieL=iface-1;
xciL=0.5*(Grid(ieL)+Grid(ieL+1));
A(ieL,ieL)=A(ieL,ieL)+4*omegab^2*((Grid(iface)-xciL)/Deltax(ieL))^2/Deltax(ieL);
R(ieL)=R(ieL)+4*omegab^2*(sin(pi)-Unumsolution(ieL))*((Grid(iface)-xciL)/Deltax(ieL))/Deltax(ieL);

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
I1=0;t=[-1/sqrt(5),0,1/sqrt(5)];W=[5/9,8/9,5/9];
k=1;%determine the correctness of the program
for K=1:numberx-1
    xci=(Grid(K+1)+Grid(K))/2;
    deltaxi=Grid(K+1)-Grid(K);
   for i=1:3
       xi=(Grid(K+1)-Grid(K))/2*t(i)+0.5*(Grid(K+1)+Grid(K));
       %对v
%       fi=(pi*cos(pi*xi)-(Unumsolution(2,K)/deltaxi+Unumsolution(3,K)*(xi-xci)/deltaxi^2+Ureconstruct(K)*((xi-xci)^2/(2*deltaxi^2)-1/24)/deltaxi))^2;
      %对phi
      fi=(sin(pi*xi)-(Unumsolution(1,K)+Ureconstruct(K)*(xi-xci)/deltaxi))^2;

       I1=I1+W(i)*fi*0.5*(Grid(K+1)-Grid(K));        
   end
   k=k+1;
end
A2=sqrt(I1);
end
