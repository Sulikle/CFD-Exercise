clc
clear all
close all
%% Preproceeding
%Some basic paramater
Unit=8;%��Ԫ����
option=1;
afa=0.06;varphi0=50;
Lr=1/(2*pi);Tr=Lr^2/afa;abslambda=sqrt(afa/Tr);A=[abslambda,0;0,abslambda];
CFLt=0.01;CFLtau=0.01;
endt=10;
endtau=10;%αʱ����ֵ
tol=10^(-5);%����ѭ������
belta=0.05;%�����Ŷ�ϵ��
endx=1;deltax=endx/Unit;numberx=endx/deltax+1;
deltat=CFLt*deltax/abslambda;
Vcurrent=zeros(3,numberx-1);
Vn=zeros(3,numberx-1);
Vm=zeros(3,numberx-1);
Vm1=zeros(3,numberx-1);
%RHS
R=zeros(3,Unit);
Rd=zeros(3,Unit);
Rb=zeros(3,Unit);
%��¼�ڵ�λ��,���¸���������belta
Grid=zeros(1,numberx);
Deltax=zeros(1,Unit);
for i=2:numberx-1
    Grid(1,i)=(i-1)*deltax+(2*rand(1)-1)*belta*deltax;
end
Grid(1,numberx)=endx;
%��¼ÿ����Ԫ�����䳤��
for i=2:numberx
        Deltax(i-1)=Grid(1,i)-Grid(1,i-1);
end

%αʱ���ϵ�local time stepping
deltatau=zeros(1,numberx-1);%αʱ�����
for i=1:numberx-1
deltatau(i)=CFLtau*(Grid(1,i+1)-Grid(1,i))/abslambda;%αʱ�����
end
%��¼ÿ������ʱ�䲽�ϵ�αʱ����ֹʱ��
n=zeros(1,floor(endt/deltat));

 %solve the exasolution
 Uexasolution=zeros(2,numberx); Vnumsolution=zeros(3,numberx);
for k=1:numberx
    Uexasolution(1,k)=varphi0*sin(pi*Grid(k))*exp(-afa*pi^2*endt);
    Uexasolution(2,k)=pi*varphi0*cos(pi*Grid(k))*exp(-afa*pi^2*endt);
end

%% Proceeding
%initial  condition set up
for k=1:numberx-1
    Vcurrent(1,k)=-varphi0*(cos(pi*Grid(k+1))-cos(pi*Grid(k)))/(pi*Deltax(k));
    Vcurrent(2,k)=varphi0*(sin(pi*Grid(k+1))-sin(pi*Grid(k)));
    Vcurrent(3,k)=-varphi0*pi^2*sin(pi*0.5*(Grid(k+1)+Grid(k)))*Deltax(k)^2;
end

%��t=0ʱ�̸�ֵ
Vn=Vcurrent;
i=1;
for itime=0:deltat:endt
    Vm=Vn;
    for itau=0:min(deltatau):endtau
        %��װRHS
        %����Gauss���ּ���Rdomain
        t=[-sqrt(15)/5,0,sqrt(15)/5];
        W=[5/9,8/9,5/9];
       for ie=1:Unit
           xci=(Grid(ie+1)+Grid(ie))/2;
           for ig=1:3
               xig=Deltax(ie)/2*t(ig)+xci;
               B2g=(xig-xci)/Deltax(ie);
               B3g=0.5*B2g^2-1/24;
               varphig=Vm(1,ie)+B2g*Vm(2,ie)+B3g*Vm(3,ie);
               Vg=Vm(2,ie)/Deltax(ie)+Vm(3,ie)*B2g/Deltax(ie);
               S1g=0;
               S2g=-Vg/Tr;
               F1g=-afa*Vg;
               F2g=-varphig/Tr;
               const=W(ig)*0.5*Deltax(ie);
               Rd(3*(ie-1)+1)=Rd(3*(ie-1)+1)+S1g*const;
               Rd(3*(ie-1)+2)=Rd(3*(ie-1)+2)+(S1g*B2g+S2g/Deltax(ie)+F1g/Deltax(ie))*const;
               Rd(3*(ie-1)+3)=Rd(3*(ie-1)+3)+(S1g*B3g+S2g*B2g/Deltax(ie)+F1g*B2g/Deltax(ie)+F2g/Deltax(ie)^2)*const;
              
           end
       end



        %Rboundary
           for iface=2:numberx-1
               ieL=iface-1;
               ieR=iface;
               B2L=1/2;B3L=0.5*B2L^2-1/24;
               B2R=-1/2;B3R=0.5*B2R^2-1/24;
               varphiL=Vm(1,ieL)+B2L*Vm(2,ieL)+B3L*Vm(3,ieL);
               varphiR=Vm(1,ieR)+B2R*Vm(2,ieR)+B3R*Vm(3,ieR);
               VL=Vm(2,ieL)/Deltax(ieL)+Vm(3,ieL)*B2L/Deltax(ieL);
               VR=Vm(2,ieR)/Deltax(ieR)+Vm(3,ieR)*B2R/Deltax(ieR);
    
               Fn(:,iface)=0.5*([-afa*VL;-varphiL/Tr]+[-afa*VR;-varphiR/Tr])-0.5*A*([varphiR;VR]-[varphiL;VL]);
               Rb(3*(ieL-1)+1)=Rb(3*(ieL-1)+1)-Fn(1,iface);
               Rb(3*(ieL-1)+2)=Rb(3*(ieL-1)+2)-Fn(1,iface)*B2L-Fn(2,iface)/Deltax(ieL);
               Rb(3*(ieL-1)+3)=Rb(3*(ieL-1)+3)-Fn(1,iface)*B3L-Fn(2,iface)*B2L/Deltax(ieL);
               Rb(3*(ieR-1)+1)=Rb(3*(ieR-1)+1)+Fn(1,iface);
               Rb(3*(ieR-1)+2)=Rb(3*(ieR-1)+2)+Fn(1,iface)*B2R+Fn(2,iface)/Deltax(ieR);
               Rb(3*(ieR-1)+3)=Rb(3*(ieR-1)+3)+Fn(1,iface)*B3R+Fn(2,iface)*B2R/Deltax(ieR);
            end
%�߽磺��
iface=1;
ieR=iface;
B2R=-1/2;B3R=0.5*B2R^2-1/24;
varphiR=Vm(1,ieR)+B2R*Vm(2,ieR)+B3R*Vm(3,ieR);
VR=Vm(2,ieR)/Deltax(ieR)+Vm(3,ieR)*B2R/Deltax(ieR);
Fn(:,iface)=0.5*([-afa*VR;0]+[-afa*VR;-varphiR/Tr])-0.5*A*([varphiR;VR]-[0;VR]);
Rb(3*(ieR-1)+1)=Rb(3*(ieR-1)+1)+Fn(1,iface);
Rb(3*(ieR-1)+2)=Rb(3*(ieR-1)+2)+Fn(1,iface)*B2R+Fn(2,iface)/Deltax(ieR);
Rb(3*(ieR-1)+3)=Rb(3*(ieR-1)+3)+Fn(1,iface)*B3R+Fn(2,iface)*B2R/Deltax(ieR);

%�߽磺��
iface=numberx;
ieL=iface-1;
B2L=1/2;B3L=0.5*B2L^2-1/24;
varphiL=Vm(1,ieL)+B2L*Vm(2,ieL)+B3L*Vm(3,ieL);
VL=Vm(2,ieL)/Deltax(ieL)+Vm(3,ieL)*B2L/Deltax(ieL);
Fn(:,iface)=0.5*([-afa*VL;-varphiL/Tr]+[-afa*VL;0])-0.5*A*([0;VL]-[varphiL;VL]);
Rb(3*(ieL-1)+1)=Rb(3*(ieL-1)+1)-Fn(1,iface);
Rb(3*(ieL-1)+2)=Rb(3*(ieL-1)+2)-Fn(1,iface)*B2L-Fn(2,iface)/Deltax(ieL);
Rb(3*(ieL-1)+3)=Rb(3*(ieL-1)+3)-Fn(1,iface)*B3L-Fn(2,iface)*B2L/Deltax(ieL);

%R��װ
R=Rd+Rb;

if option==1
    for k=1:numberx-1
        Mtau=[Deltax(k),0,0;0,Deltax(k)/12+1/Deltax(k),0;0,0,1/720*Deltax(k)+1/(12*Deltax(k))];
        Mt=[Deltax(k),0,0;0,Deltax(k)/12,0;0,0,1/720*Deltax(k)];
        Vm1(:,k)=Vm(:,k)+Mtau\(R(:,k)-Mt*(Vm(:,k)-Vn(:,k))/deltat)*deltatau(k);
    end
end
if option==2
    for k=1:numberx-1
        Mtau=[Deltax(k),0,0;0,Deltax(k)/12+1/Deltax(k),0;0,0,1/720*Deltax(k)+1/(12*Deltax(k))];
        Mt=[Deltax(k),0,0;0,Deltax(k)/12,0;0,0,1/720*Deltax(k)];
        Vm1(:,k)=(Mtau/deltatau(k)+Mt/deltat)\(R(:,k)+Mt*Vn(:,k)/deltat+Mtau*Vm(:,k)/deltatau(k));
    end
end


        if max(max(abs(Vm1-Vm)))<tol &&itau>=min(deltatau)
            break;
        end
        
        Vm=Vm1;
        Rd=zeros(3,Unit);
        Rb=zeros(3,Unit);
        
    end
    n(i)=itau;i=i+1;    
    Vn=Vm1;    
end
Vnumsolution=Vn;

%% Post-proceeding
%U
figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Vnumsolution(1,k)+Vnumsolution(2,k)*(x-xci)/Deltax(k)+Vnumsolution(3,k)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on

 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Vnumsolution(1,k)+Vnumsolution(2,k)*(x-xci)/Deltax(k)+Vnumsolution(3,k)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);
 end
 
%plot the exact
plot(Grid,Uexasolution(1,:),'-b*','linewidth',1.5)
H2=plot(Grid,Uexasolution(1,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DG(P0P2)+DG(P0P1)','������');
lgd.FontSize=12;
xlabel('λ��x','fontsize',14)
ylabel('��ֵU','fontsize',14)
 title('DG(P0P2)+DG(P0P1) Explicit Euler��ֵ���������(U)','fontsize',16)
hold off

% Ux
 figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Vnumsolution(2,k)/Deltax(k)+Vnumsolution(3,k)*(x-xci)/Deltax(k)^2;
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on
 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
      xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Vnumsolution(2,k)/Deltax(k)+Vnumsolution(3,k)*(x-xci)/Deltax(k)^2;
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);
 end

%exact
plot(Grid,Uexasolution(2,:),'-b*','linewidth',1.5)
H2=plot(Grid,Uexasolution(2,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DG(P0P2)+DG(P0P1)','������');
lgd.FontSize=12;
xlabel('λ��x','fontsize',14)
ylabel('��ֵU','fontsize',14)
title('DG(P0P2)+DG(P0P1) Explicit Euler��ֵ���������(Ux)','fontsize',16)
hold off






