function [Unumsolution,n]=subDGP0P1plusDGP0(Unit,CFL,endtau)
deltax=1/Unit;tol=0.0001;
nu=1;Lr=1/(2*pi);Tr=Lr^2/nu;
abslambda=sqrt(nu/Tr);deltatau=CFL*deltax/abslambda;%伪时间变量
endx=1;
numberx=endx/deltax+1;
Ucurrent=zeros(2,numberx-1);
Unext=zeros(2,numberx-1);
B1=1;
%C=[B1,0;0,B1/deltax];
Mtau=[deltax,0;0,deltax/12+1/deltax];%此为推导出的U=CV中的C
A=[abslambda,0;0,abslambda];
R=zeros(2,numberx-1);
Rd=zeros(2,numberx-1);
Rb=zeros(2,numberx-1);
Fn=zeros(2,numberx);
V0=zeros(1,numberx-1);
%% solve the question
%initial  condition set up
x=0;
for k=1:numberx-1
    Ucurrent(1,k)=(x+deltax/2)^2-(x+deltax/2);
    x=x+deltax;
end
V0=Ucurrent(1,:);

x=0;
for k=1:numberx-1
    Ucurrent(2,k)=(2*(x+deltax/2)-1)*deltax;%Ucurrent(2,:)存储的是Ux*deltx
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
    R(:,k)=Rd(:,k)+Rb(:,k);
end


%solve the numsolution
for n=deltatau:deltatau:endtau
    for k=1:numberx-1
       Unext(:,k)=Ucurrent(:,k)+Mtau\R(:,k)*deltatau;
    end
     if var(Ucurrent(1,:)-Unext(1,:))<tol*V0
         break
     end
    Ucurrent=Unext;   
    Rd=zeros(2,numberx-1);
    Rb=zeros(2,numberx-1);
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
    R(:,k)=Rd(:,k)+Rb(:,k);
end
       
end

Unumsolution(1,:)=Ucurrent(1,:);Unumsolution(2,:)=Ucurrent(2,:)/deltax;


end