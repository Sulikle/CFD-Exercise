function [Unumsolution,Grid,N]=subDGP0P1plusDGP0(Unit,CFL,endtau)
deltax=1/Unit;tol=0.0001;
nu=1;Lr=1/(2*pi);Tr=Lr^2/nu;
abslambda=sqrt(nu/Tr);
endx=1;
numberx=Unit+1;
%记录内点位置
Grid=zeros(1,numberx);
for i=2:numberx-1
    Grid(1,i)=(i-1)*deltax+(0.1*rand(1)-0.05)*deltax;
end
Grid(1,numberx)=endx;

deltatau=zeros(1,numberx-1);
for i=1:numberx-1
deltatau(i)=CFL*(Grid(1,i+1)-Grid(1,i))/abslambda;%伪时间变量
end

Ucurrent=zeros(2,numberx-1);
Unext=zeros(2,numberx-1);
%B1=1;
%C=[B1,0;0,B1/deltax];
%Mtau=[deltax,0;0,deltax/12+1/deltax];%此为推导出的U=CV中的C
A=[abslambda,0;0,abslambda];
R=zeros(2,numberx-1);
Rd=zeros(2,numberx-1);
Rb=zeros(2,numberx-1);
Fn=zeros(2,numberx);
V0=zeros(1,numberx-1);
%% solve the question
%initial  condition set up

for k=1:numberx-1
        x=Grid(k);
    Ucurrent(1,k)=(x+(Grid(k+1)-Grid(k))/2)^2-(x+(Grid(k+1)-Grid(k))/2);
    Ucurrent(2,k)=(2*(x+(Grid(k+1)-Grid(k))/2)-1)*(Grid(k+1)-Grid(k));%Ucurrent(2,:)存储的是Ux*deltx
end
V0=Ucurrent(1,:);
V1=Ucurrent(2,:);



%Rdomain
for k=1:numberx-1
    Rd(1,k)=gauss1(Grid(k),Grid(k+1));
    %Rd(1,k)=pi*(cos(pi*x)-cos(pi*(x+deltax)));
    %Rd(2,k)=-nu*Ucurrent(2,k)/deltax+(-pi/deltax*(deltax/2*cos(pi*(x+deltax))-(-deltax/2)*cos(pi*x)-1/pi*(sin(pi*(x+deltax))-sin(pi*x))))-Ucurrent(2,k)/(Tr*deltax);
    Rd(2,k)=-nu*Ucurrent(2,k)/(Grid(k+1)-Grid(k))-Ucurrent(2,k)/(Tr*(Grid(k+1)-Grid(k)))+gauss2(Grid(k),Grid(k+1));
end

%Rboundary
for iface=2:numberx-1
    ieL=iface-1;
    ieR=iface;
    B2L=1/2;
    B2R=-1/2;
    Fn(:,iface)=0.5*([-nu*Ucurrent(2,ieL)/(Grid(ieL+1)-Grid(ieL));-(Ucurrent(1,ieL)+B2L*Ucurrent(2,ieL))/Tr]+[-nu*Ucurrent(2,ieR)/(Grid(ieR+1)-Grid(ieR));-(Ucurrent(1,ieR)+B2R*Ucurrent(2,ieR))/Tr])-0.5*A*([Ucurrent(1,ieR)+B2R*Ucurrent(2,ieR);Ucurrent(2,ieR)/(Grid(ieR+1)-Grid(ieR))]-[Ucurrent(1,ieL)+B2L*Ucurrent(2,ieL);Ucurrent(2,ieL)/(Grid(ieL+1)-Grid(ieL))]);
    Rb(:,ieL)=Rb(:,ieL)-[1,B2L;0,1/(Grid(ieL+1)-Grid(ieL))]'*Fn(:,iface);
    Rb(:,ieR)=Rb(:,ieR)+[1,B2R;0,1/(Grid(ieR+1)-Grid(ieR))]'*Fn(:,iface);
end
Fn(:,1)=0.5*([-nu*Ucurrent(2,1)/(Grid(1+1)-Grid(1));0]+[-nu*Ucurrent(2,1)/(Grid(1+1)-Grid(1));-(Ucurrent(1,1)+B2R*Ucurrent(2,1))/Tr])-0.5*A*([Ucurrent(1,1)+B2R*Ucurrent(2,1);Ucurrent(2,1)/(Grid(1+1)-Grid(1))]-[0;Ucurrent(2,1)/(Grid(1+1)-Grid(1))]);
Fn(:,numberx)=0.5*([-nu*Ucurrent(2,numberx-1)/(Grid(numberx)-Grid(numberx-1));-(Ucurrent(1,numberx-1)+B2L*Ucurrent(2,numberx-1))/Tr]+[-nu*Ucurrent(2,numberx-1)/(Grid(numberx)-Grid(numberx-1));0])-0.5*A*([0;Ucurrent(2,numberx-1)/(Grid(numberx)-Grid(numberx-1))]-[Ucurrent(1,numberx-1)+B2L*Ucurrent(2,numberx-1);Ucurrent(2,numberx-1)/(Grid(numberx)-Grid(numberx-1))]);
Rb(:,1)=Rb(:,1)+[1,B2R;0,1/(Grid(2)-Grid(1))]'*Fn(:,1);
Rb(:,numberx-1)=Rb(:,numberx-1)-[1,B2L;0,1/(Grid(numberx)-Grid(numberx-1))]'*Fn(:,numberx);

%R组装
for k=1:numberx-1
    R(:,k)=Rd(:,k)+Rb(:,k);
end


%solve the numsolution

for n=1:floor(endtau/max(deltatau))
    for k=1:numberx-1
        deltax=(Grid(k+1)-Grid(k));
       Unext(:,k)=Ucurrent(:,k)+[deltax,0;0,deltax/12+1/deltax]\R(:,k)*deltatau(k);
    end
     if abs(Ucurrent(1,:)-Unext(1,:))<tol*V0 & abs(Ucurrent(2,:)-Unext(2,:))<tol*V1
         N=max(n*deltatau);
         break
     end
    Ucurrent=Unext;   
    Rd=zeros(2,numberx-1);
    Rb=zeros(2,numberx-1);
%Rdomain
for k=1:numberx-1
    Rd(1,k)=gauss1(Grid(k),Grid(k+1));
    %Rd(1,k)=pi*(cos(pi*x)-cos(pi*(x+deltax)));
    %Rd(2,k)=-nu*Ucurrent(2,k)/deltax+(-pi/deltax*(deltax/2*cos(pi*(x+deltax))-(-deltax/2)*cos(pi*x)-1/pi*(sin(pi*(x+deltax))-sin(pi*x))))-Ucurrent(2,k)/(Tr*deltax);
    Rd(2,k)=-nu*Ucurrent(2,k)/(Grid(k+1)-Grid(k))-Ucurrent(2,k)/(Tr*(Grid(k+1)-Grid(k)))+gauss2(Grid(k),Grid(k+1));
end

%Rboundary
for iface=2:numberx-1
    ieL=iface-1;
    ieR=iface;
    B2L=1/2;
    B2R=-1/2;
    Fn(:,iface)=0.5*([-nu*Ucurrent(2,ieL)/(Grid(ieL+1)-Grid(ieL));-(Ucurrent(1,ieL)+B2L*Ucurrent(2,ieL))/Tr]+[-nu*Ucurrent(2,ieR)/(Grid(ieR+1)-Grid(ieR));-(Ucurrent(1,ieR)+B2R*Ucurrent(2,ieR))/Tr])-0.5*A*([Ucurrent(1,ieR)+B2R*Ucurrent(2,ieR);Ucurrent(2,ieR)/(Grid(ieR+1)-Grid(ieR))]-[Ucurrent(1,ieL)+B2L*Ucurrent(2,ieL);Ucurrent(2,ieL)/(Grid(ieL+1)-Grid(ieL))]);
    Rb(:,ieL)=Rb(:,ieL)-[1,B2L;0,1/(Grid(ieL+1)-Grid(ieL))]'*Fn(:,iface);
    Rb(:,ieR)=Rb(:,ieR)+[1,B2R;0,1/(Grid(ieR+1)-Grid(ieR))]'*Fn(:,iface);
end
Fn(:,1)=0.5*([-nu*Ucurrent(2,1)/(Grid(1+1)-Grid(1));0]+[-nu*Ucurrent(2,1)/(Grid(1+1)-Grid(1));-(Ucurrent(1,1)+B2R*Ucurrent(2,1))/Tr])-0.5*A*([Ucurrent(1,1)+B2R*Ucurrent(2,1);Ucurrent(2,1)/(Grid(1+1)-Grid(1))]-[0;Ucurrent(2,1)/(Grid(1+1)-Grid(1))]);
Fn(:,numberx)=0.5*([-nu*Ucurrent(2,numberx-1)/(Grid(numberx)-Grid(numberx-1));-(Ucurrent(1,numberx-1)+B2L*Ucurrent(2,numberx-1))/Tr]+[-nu*Ucurrent(2,numberx-1)/(Grid(numberx)-Grid(numberx-1));0])-0.5*A*([0;Ucurrent(2,numberx-1)/(Grid(numberx)-Grid(numberx-1))]-[Ucurrent(1,numberx-1)+B2L*Ucurrent(2,numberx-1);Ucurrent(2,numberx-1)/(Grid(numberx)-Grid(numberx-1))]);
Rb(:,1)=Rb(:,1)+[1,B2R;0,1/(Grid(2)-Grid(1))]'*Fn(:,1);
Rb(:,numberx-1)=Rb(:,numberx-1)-[1,B2L;0,1/(Grid(numberx)-Grid(numberx-1))]'*Fn(:,numberx);

%R组装
for k=1:numberx-1
    R(:,k)=Rd(:,k)+Rb(:,k);
end

end
if n==floor(endtau/max(deltatau))
    N=max(n*deltatau);
end

Unumsolution(1,:)=Ucurrent(1,:);
for k=1:numberx-1
Unumsolution(2,k)=Ucurrent(2,k)/(Grid(k+1)-Grid(k));

end

end