function [Unumsolution,n]=subDGP0P1plusDGP0(Unit,CFL,endtau)
%Some basic paramater
endx=1;deltax=endx/Unit;numberx=endx/deltax+1;
tol=10^(-10);
nu=1;Lr=1/(2*pi);Tr=Lr^2/nu;
abslambda=sqrt(nu/Tr);deltatau=CFL*deltax/abslambda;%伪时间变量
B1=1;
%C=[B1,0;0,B1/deltax];
Mtau=[deltax,0;0,deltax/12+1/deltax];%此为推导出的U=CV中的C
A=[abslambda,0;0,abslambda];
R=zeros(2*Unit,1);
Rd=zeros(2,numberx-1);
Rb=zeros(2,numberx-1);
Fn=zeros(2,numberx);

%构建LHS
%Mtau/deltatau
LHS1=sparse(1:2:2*Unit-1,1:2:2*Unit-1,deltax/deltatau,2*Unit,2*Unit);
LHS1=LHS1+sparse(2:2:2*Unit,2:2:2*Unit,(deltax/12+1/deltax)/deltatau,2*Unit,2*Unit);
%Rdomain
LHS2=-sparse(2:2:2*Unit,2:2:2*Unit,-(nu+1/Tr)/deltax,2*Unit,2*Unit);
%Rboundary
LHS3=zeros(2*Unit,2*Unit);
for iface=2:numberx-1
    ieL=iface-1;
    ieR=iface;
    CL=[B1,1/2;0,B1/deltax];
    CR=[B1,-1/2;0,B1/deltax];
    %diag
    LHS3(2*ieL-1:2*ieL,2*ieL-1:2*ieL)=LHS3(2*ieL-1:2*ieL,2*ieL-1:2*ieL)+CL'*[abslambda/2,-nu/2;-1/(2*Tr),abslambda/2]*CL;
    LHS3(2*ieR-1:2*ieR,2*ieR-1:2*ieR)=LHS3(2*ieR-1:2*ieR,2*ieR-1:2*ieR)-CR'*[-abslambda/2,-nu/2;-1/(2*Tr),-abslambda/2]*CR;
    %upper
    LHS3(2*ieL-1:2*ieL,2*ieR-1:2*ieR)=LHS3(2*ieL-1:2*ieL,2*ieR-1:2*ieR)+CL'*[-abslambda/2,-nu/2;-1/(2*Tr),-abslambda/2]*CR;
    %lower
    LHS3(2*ieR-1:2*ieR,2*ieL-1:2*ieL)=LHS3(2*ieR-1:2*ieR,2*ieL-1:2*ieL)-CR'*[abslambda/2,-nu/2;-1/(2*Tr),abslambda/2]*CL;
end
    LHS3(2*1-1:2*1,2*1-1:2*1)=LHS3(2*1-1:2*1,2*1-1:2*1)-CR'*([abslambda/2,-nu/2;-1/(2*Tr),abslambda/2]*[0,0;0,1/deltax]+[-abslambda/2,-nu/2;-1/(2*Tr),-abslambda/2]*CR);
    LHS3(2*(numberx-1)-1:2*(numberx-1),2*(numberx-1)-1:2*(numberx-1))=LHS3(2*(numberx-1)-1:2*(numberx-1),2*(numberx-1)-1:2*(numberx-1))+CL'*([abslambda/2,-nu/2;-1/(2*Tr),abslambda/2]+[-abslambda/2,-nu/2;-1/(2*Tr),-abslambda/2]*[0,0;0,1])*CL;

LHS=LHS1+LHS2+LHS3;

%取出我们所需要的D
D=zeros(2*Unit,2*Unit);
for iface=2:numberx
    ieL=iface-1;
    D(2*ieL-1:2*ieL,2*ieL-1:2*ieL)=LHS(2*ieL-1:2*ieL,2*ieL-1:2*ieL);
end
%为循环所预设的一些量
Ucurrent=zeros(2,numberx-1);
Unext=zeros(2*Unit,1);
%initial  condition set up
x=0;
for k=1:numberx-1
    Ucurrent(1,k)=(x+deltax/2)^2-(x+deltax/2);
    Ucurrent(2,k)=(2*(x+deltax/2)-1)*deltax;
    x=x+deltax;
end

%Rdomain
x=0;
for k=1:numberx-1
    Rd(1,k)=gauss1((k-1)*deltax,k*deltax);
    %Rd(1,k)=pi*(cos(pi*x)-cos(pi*(x+deltax)));
    %Rd(2,k)=-nu*Ucurrent(2,k)/deltax+(-pi/deltax*(deltax/2*cos(pi*(x+deltax))-(-deltax/2)*cos(pi*x)-1/pi*(sin(pi*(x+deltax))-sin(pi*x))))-Ucurrent(2,k)/(Tr*deltax);
    Rd(2,k)=-nu*Ucurrent(2,k)/deltax-Ucurrent(2,k)/(Tr*deltax)+gauss2((k-1)*deltax,k*deltax,k);
    x=x+deltax;
end
%Rboundary
for iface=2:numberx-1
    ieL=iface-1;
    ieR=iface;
    B2L=1/2;
    B2R=-1/2;
    Fn(:,iface)=0.5*([-nu*Ucurrent(2,ieL)/deltax;-(Ucurrent(1,ieL)+B2L*Ucurrent(2,ieL))/Tr]+[-nu*Ucurrent(2,ieR)/deltax;-(Ucurrent(1,ieR)+B2R*Ucurrent(2,ieR))/Tr])-0.5*A*([Ucurrent(1,ieR)+B2R*Ucurrent(2,ieR);Ucurrent(2,ieR)/deltax]-[Ucurrent(1,ieL)+B2L*Ucurrent(2,ieL);Ucurrent(2,ieL)/deltax]);
    Rb(:,ieL)=Rb(:,ieL)-[1,B2L;0,1/deltax]'*Fn(:,iface);
    Rb(:,ieR)=Rb(:,ieR)+[1,B2R;0,1/deltax]'*Fn(:,iface);
end
Fn(:,1)=0.5*([-nu*Ucurrent(2,1)/deltax;0]+[-nu*Ucurrent(2,1)/deltax;-(Ucurrent(1,1)+B2R*Ucurrent(2,1))/Tr])-0.5*A*([Ucurrent(1,1)+B2R*Ucurrent(2,1);Ucurrent(2,1)/deltax]-[0;Ucurrent(2,1)/deltax]);
Fn(:,numberx)=0.5*([-nu*Ucurrent(2,numberx-1)/deltax;-(Ucurrent(1,numberx-1)+B2L*Ucurrent(2,numberx-1))/Tr]+[-nu*Ucurrent(2,numberx-1)/deltax;0])-0.5*A*([0;Ucurrent(2,numberx-1)/deltax]-[Ucurrent(1,numberx-1)+B2L*Ucurrent(2,numberx-1);Ucurrent(2,numberx-1)/deltax]);
Rb(:,1)=Rb(:,1)+[1,B2R;0,1/deltax]'*Fn(:,1);
Rb(:,numberx-1)=Rb(:,numberx-1)-[1,B2L;0,1/deltax]'*Fn(:,numberx);

%R组装
for k=1:numberx-1
    R(2*k-1:2*k,1)=Rd(:,k)+Rb(:,k);
end

%进行必要的向量等价转变
for k=1:numberx-1
    Unext(2*k-1:2*k,1)=Ucurrent(:,k);
end
%循环迭代
for n=deltatau:deltatau:endtau
    X=D\R;
    if max(X)<tol
        break
    end
    Unext=Unext+X;
    Rd=zeros(2,numberx-1);
    Rb=zeros(2,numberx-1);    
for k=1:numberx-1
    Ucurrent(:,k)=Unext(2*k-1:2*k,1);
end
%Rdomain
x=0;
for k=1:numberx-1
    Rd(1,k)=gauss1((k-1)*deltax,k*deltax);
    %Rd(1,k)=pi*(cos(pi*x)-cos(pi*(x+deltax)));
    %Rd(2,k)=-nu*Ucurrent(2,k)/deltax+(-pi/deltax*(deltax/2*cos(pi*(x+deltax))-(-deltax/2)*cos(pi*x)-1/pi*(sin(pi*(x+deltax))-sin(pi*x))))-Ucurrent(2,k)/(Tr*deltax);
    Rd(2,k)=-nu*Ucurrent(2,k)/deltax-Ucurrent(2,k)/(Tr*deltax)+gauss2((k-1)*deltax,k*deltax,k);
    x=x+deltax;
end
%Rboundary
for iface=2:numberx-1
    ieL=iface-1;
    ieR=iface;
    B2L=1/2;
    B2R=-1/2;
    Fn(:,iface)=0.5*([-nu*Ucurrent(2,ieL)/deltax;-(Ucurrent(1,ieL)+B2L*Ucurrent(2,ieL))/Tr]+[-nu*Ucurrent(2,ieR)/deltax;-(Ucurrent(1,ieR)+B2R*Ucurrent(2,ieR))/Tr])-0.5*A*([Ucurrent(1,ieR)+B2R*Ucurrent(2,ieR);Ucurrent(2,ieR)/deltax]-[Ucurrent(1,ieL)+B2L*Ucurrent(2,ieL);Ucurrent(2,ieL)/deltax]);
    Rb(:,ieL)=Rb(:,ieL)-[1,B2L;0,1/deltax]'*Fn(:,iface);
    Rb(:,ieR)=Rb(:,ieR)+[1,B2R;0,1/deltax]'*Fn(:,iface);
end
Fn(:,1)=0.5*([-nu*Ucurrent(2,1)/deltax;0]+[-nu*Ucurrent(2,1)/deltax;-(Ucurrent(1,1)+B2R*Ucurrent(2,1))/Tr])-0.5*A*([Ucurrent(1,1)+B2R*Ucurrent(2,1);Ucurrent(2,1)/deltax]-[0;Ucurrent(2,1)/deltax]);
Fn(:,numberx)=0.5*([-nu*Ucurrent(2,numberx-1)/deltax;-(Ucurrent(1,numberx-1)+B2L*Ucurrent(2,numberx-1))/Tr]+[-nu*Ucurrent(2,numberx-1)/deltax;0])-0.5*A*([0;Ucurrent(2,numberx-1)/deltax]-[Ucurrent(1,numberx-1)+B2L*Ucurrent(2,numberx-1);Ucurrent(2,numberx-1)/deltax]);
Rb(:,1)=Rb(:,1)+[1,B2R;0,1/deltax]'*Fn(:,1);
Rb(:,numberx-1)=Rb(:,numberx-1)-[1,B2L;0,1/deltax]'*Fn(:,numberx);

%R组装
for k=1:numberx-1
    R(2*k-1:2*k,1)=Rd(:,k)+Rb(:,k);
end

for k=1:numberx-1
    Unext(2*k-1:2*k,1)=Ucurrent(:,k);
end
end
Unumsolution(1,:)=Ucurrent(1,:);Unumsolution(2,:)=Ucurrent(2,:)/deltax;
end
