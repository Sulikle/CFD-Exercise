function [Ul2error,Uxl2error]=Accuracy(Unit)
%% Pre-processing
%Some basic paramater
nu=0.01;a=10^4;x0=0.5;
Lr=1/max(2*pi,abs(a)/nu);Tr=Lr^2/nu;abslambda=sqrt(nu/Tr);A1=[abslambda+abs(a),0;0,abslambda];A=[a,-nu;-1/Tr,0];B1=1;
deltat=10^(-9);
CFLtau=100;
endtau=10;%伪时间阈值
endt=10^(-6);
tol=10^(-5);%跳出循环条件
belta=0;%网格扰动系数
endx=2;deltax=endx/Unit;numberx=endx/deltax+1;
Vcurrent=zeros(2,numberx-1);
Vn=zeros(2,numberx-1);
Vm=zeros(2,numberx-1);
Vm1=zeros(2,numberx-1);
dimention=2;
%LHS
LHS1=zeros(dimention*Unit,dimention*Unit);
LHS2=zeros(dimention*Unit,dimention*Unit);
LHS3=zeros(dimention*Unit,dimention*Unit);
%RHS
R=zeros(dimention*Unit,1);
Rd=zeros(dimention*Unit,1);
Rb=zeros(dimention*Unit,1);
Fn=zeros(2,numberx);
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
deltatau(i)=CFLtau*(Grid(1,i+1)-Grid(1,i))/(abslambda+abs(a));%伪时间变量
end
%记录每个物理时间步上的伪时间终止时刻
n=zeros(1,floor(endt/deltat));
 %solve the exasolution
 Uexasolution=zeros(2,numberx); Vnumsolution=zeros(2,numberx);
for k=1:numberx
    Uexasolution(1,k)=1/sqrt(4*endt+1)*exp(-(Grid(k)-a*endt-x0)^2/(nu*(4*endt+1)));
    Uexasolution(2,k)=1/sqrt(4*endt+1)*exp(-(Grid(k)-a*endt-x0)^2/(nu*(4*endt+1)))*(-2*(Grid(k)-a*endt-x0)/(nu*(4*endt+1)));
end

%% Proceeding

%构建LHS
%Mtau/deltatau
for i=1:Unit
    LHS1(dimention*(i-1)+1,dimention*(i-1)+1)=Deltax(i)/deltatau(i);
    LHS1(dimention*(i-1)+2,dimention*(i-1)+2)=(Deltax(i)/12+1/Deltax(i))/deltatau(i);
end
%Rdomain
for i=1:Unit
    LHS2(dimention*(i-1)+2,dimention*(i-1)+2)=(1/Tr+nu)/Deltax(i);
    LHS2(dimention*(i-1)+2,dimention*(i-1)+1)=-a;
end

%Rboundary
for iface=2:numberx-1
    ieL=iface-1;
    ieR=iface;
    CL=[B1,1/2;0,B1/Deltax(ieL)];
    CR=[B1,-1/2;0,B1/Deltax(ieR)];
    %diag
    LHS3(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieL-1)+1:dimention*ieL)=LHS3(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieL-1)+1:dimention*ieL)+0.5*CL'*(A+A1)*CL;
    LHS3(dimention*(ieR-1)+1:dimention*ieR,dimention*(ieR-1)+1:dimention*ieR)=LHS3(dimention*(ieR-1)+1:dimention*ieR,dimention*(ieR-1)+1:dimention*ieR)-0.5*CR'*(A-A1)*CR;
    %upper
    LHS3(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieR-1)+1:dimention*ieR)=LHS3(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieR-1)+1:dimention*ieR)+0.5*CL'*(A-A1)*CR;
    %lower
    LHS3(dimention*(ieR-1)+1:dimention*ieR,dimention*(ieL-1)+1:dimention*ieL)=LHS3(dimention*(ieR-1)+1:dimention*ieR,dimention*(ieL-1)+1:dimention*ieL)-0.5*CR'*(A+A1)*CL;
end

%边界：左
iface=1;
ieR=iface;
CR=[B1,-1/2;0,B1/Deltax(ieR)];
LHS3(dimention*(ieR-1)+1:dimention*ieR,dimention*(ieR-1)+1:dimention*ieR)=LHS3(dimention*(ieR-1)+1:dimention*ieR,dimention*(ieR-1)+1:dimention*ieR)-CR'*(0.5*(A+A1)*[0,0;0,1]*CR+0.5*(A-A1)*CR);
%边界：右
iface=numberx;
ieL=iface-1;
CL=[B1,1/2;0,B1/Deltax(ieL)];
LHS3(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieL-1)+1:dimention*ieL)=LHS3(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieL-1)+1:dimention*ieL)+CL'*(0.5*(A+A1)+0.5*(A-A1)*[0,0;0,1])*CL;

%组装LHS
LHS=LHS1+LHS2+LHS3;

%Mt/deltat
for i=1:Unit
    LHS(dimention*(i-1)+1,dimention*(i-1)+1)=LHS(dimention*(i-1)+1,dimention*(i-1)+1)+Deltax(i)/deltat;
    LHS(dimention*(i-1)+2,dimention*(i-1)+2)=LHS(dimention*(i-1)+2,dimention*(i-1)+2)+(Deltax(i)/12)/deltat;
end


%initial  condition set up
for k=1:numberx-1
        %利用Gauss积分计算Vcurrent
        t=[-sqrt(15)/5,0,sqrt(15)/5];
        W=[5/9,8/9,5/9];
           xci=(Grid(k+1)+Grid(k))/2;
           for ig=1:3
               xig=Deltax(k)/2*t(ig)+xci;
               const=W(ig)*0.5*Deltax(k);
          Vcurrent(1,k)=Vcurrent(1,k)+const*exp(-(xig-x0)^2/nu);
          Vcurrent(2,k)=Vcurrent(2,k)+const*exp(-(xig-x0)^2/nu)*(-2*(xig-x0)/nu);
           end 
           Vcurrent(1,k)=Vcurrent(1,k)/Deltax(k);
       

end
%对t=0时刻赋值
Vn=Vcurrent;
i=1;
for itime=0:deltat:endt
    Vm=Vn;
for itau=0:min(deltatau):endtau
        %组装RHS
        %利用Gauss积分计算Rdomain
        t=[-sqrt(15)/5,0,sqrt(15)/5];
        W=[5/9,8/9,5/9];
       for ie=1:Unit
           xci=(Grid(ie+1)+Grid(ie))/2;
           for ig=1:3
               xig=Deltax(ie)/2*t(ig)+xci;
               B2g=(xig-xci)/Deltax(ie);
               varphig=Vm(1,ie)+B2g*Vm(2,ie);
               Vg=Vm(2,ie)/Deltax(ie);
               S1g=0;
               S2g=-Vg/Tr;
               F1g=a*varphig-nu*Vg;
               F2g=-varphig/Tr;
               const=W(ig)*0.5*Deltax(ie);
               Rd(2*(ie-1)+1)=Rd(2*(ie-1)+1)+S1g*const;
               Rd(2*(ie-1)+2)=Rd(2*(ie-1)+2)+(S1g*B2g+S2g/Deltax(ie)+F1g/Deltax(ie))*const;
%         Rd(dimention*(ie-1)+3)=Rd(dimention*(ie-1)+3)+(S1g*B3g+S2g*B2g/Deltax(ie)+F1g*B2g/Deltax(ie)+F2g/Deltax(ie)^2)*const;
              
           end
       end



        %Rboundary
           for iface=2:numberx-1
               ieL=iface-1;
               ieR=iface;
               B2L=1/2;
               B2R=-1/2;
               varphiL=Vm(1,ieL)+B2L*Vm(2,ieL);
               varphiR=Vm(1,ieR)+B2R*Vm(2,ieR);
               VL=Vm(2,ieL)/Deltax(ieL);
               VR=Vm(2,ieR)/Deltax(ieR);
    
               Fn(:,iface)=0.5*([a*varphiL-nu*VL;-varphiL/Tr]+[a*varphiR-nu*VR;-varphiR/Tr])-0.5*A1*([varphiR;VR]-[varphiL;VL]);
               Rb(2*(ieL-1)+1)=Rb(2*(ieL-1)+1)-Fn(1,iface);
               Rb(2*(ieL-1)+2)=Rb(2*(ieL-1)+2)-Fn(1,iface)*B2L-Fn(2,iface)/Deltax(ieL);
%     Rb(dimention*(ieL-1)+3)=Rb(dimention*(ieL-1)+3)-Fn(1,iface)*B3L-Fn(2,iface)*B2L/Deltax(ieL);
               Rb(2*(ieR-1)+1)=Rb(2*(ieR-1)+1)+Fn(1,iface);
               Rb(2*(ieR-1)+2)=Rb(2*(ieR-1)+2)+Fn(1,iface)*B2R+Fn(2,iface)/Deltax(ieR);
%     Rb(dimention*(ieR-1)+3)=Rb(dimention*(ieR-1)+3)+Fn(1,iface)*B3R+Fn(2,iface)*B2R/Deltax(ieR);
            end
%边界：左
iface=1;
ieR=iface;
B2R=-1/2;
varphiR=Vm(1,ieR)+B2R*Vm(2,ieR);
VR=Vm(2,ieR)/Deltax(ieR);
Fn(:,iface)=0.5*([-nu*VR;0]+[a*varphiR-nu*VR;-varphiR/Tr])-0.5*A1*([varphiR;VR]-[0;VR]);
Rb(2*(ieR-1)+1)=Rb(2*(ieR-1)+1)+Fn(1,iface);
Rb(2*(ieR-1)+2)=Rb(2*(ieR-1)+2)+Fn(1,iface)*B2R+Fn(2,iface)/Deltax(ieR);
% Rb(dimention*(ieR-1)+3)=Rb(dimention*(ieR-1)+3)+Fn(1,iface)*B3R+Fn(2,iface)*B2R/Deltax(ieR);

%边界：右
iface=numberx;
ieL=iface-1;
B2L=1/2;
varphiL=Vm(1,ieL)+B2L*Vm(2,ieL);
VL=Vm(2,ieL)/Deltax(ieL);
Fn(:,iface)=0.5*([a*varphiL-nu*VL;-varphiL/Tr]+[-nu*VL;0])-0.5*A1*([0;VL]-[varphiL;VL]);
Rb(2*(ieL-1)+1)=Rb(2*(ieL-1)+1)-Fn(1,iface);
Rb(2*(ieL-1)+2)=Rb(2*(ieL-1)+2)-Fn(1,iface)*B2L-Fn(2,iface)/Deltax(ieL);
% Rb(dimention*(ieL-1)+3)=Rb(dimention*(ieL-1)+3)-Fn(1,iface)*B3L-Fn(2,iface)*B2L/Deltax(ieL);

%R组装
R=Rd+Rb;

for k=1:Unit
    Mt=[Deltax(k),0;0,Deltax(k)/12];
    R(2*k-1:2*k,1)=R(2*k-1:2*k,1)-Mt*(Vm(:,k)-Vn(:,k))/deltat;
end

    X=LUSGS(LHS,R,Unit);
    if max(abs(X))<tol&&itau>=min(deltatau)
        break;
    end
for k=1:Unit
    Vm1(:,k)=Vm(:,k)+X(2*k-1:2*k,1);
end
        
        Vm=Vm1;
        Rd=zeros(dimention*Unit,1);
        Rb=zeros(dimention*Unit,1);
        
    end
    n(i)=itau;i=i+1;    
    Vn=Vm1;    
end
    
Vnumsolution=Vn;



%calculate the accuracy of spaceU
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
      %对phi
       fi=(sin(pi*xi)-(Vnumsolution(1,K)+Vnumsolution(2,K)*(xi-xci)/deltaxi))^2;
       I1=I1+W(i)*fi*0.5*(Grid(K+1)-Grid(K));        
   end
   k=k+1;
end
Ul2error=sqrt(I1);

%calculate the accuracy of spaceUx
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
      %对phi
       fi=(pi*cos(pi*xi)-(Vnumsolution(2,K)/deltaxi))^2;
       I1=I1+W(i)*fi*0.5*(Grid(K+1)-Grid(K));        
   end
   k=k+1;
end
Uxl2error=sqrt(I1);

end
