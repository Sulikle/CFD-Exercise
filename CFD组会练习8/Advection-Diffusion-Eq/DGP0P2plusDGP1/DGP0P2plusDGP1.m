clc
clear all
close all
%% Preproceeding
%Some basic paramater
Unit=8;%��Ԫ����
nu=1;a=1;
Lr=1/max(2*pi,abs(a)/nu);Tr=Lr^2/nu;abslambda=sqrt(nu/Tr);A=[abslambda+abs(a),0;0,abslambda];
CFLtau=0.01;
endtau=10;%αʱ����ֵ
tol=10^(-8);%����ѭ������
belta=0.05;%�����Ŷ�ϵ��
endx=1;deltax=endx/Unit;numberx=endx/deltax+1;
Vcurrent=zeros(3,numberx-1);
Vnext=zeros(3,numberx-1);
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
deltatau(i)=CFLtau*(Grid(1,i+1)-Grid(1,i))/(abslambda+abs(a));%αʱ�����
end

%����order
Acc=zeros(3,4);a1=[1/(2*8),1/(2*16),1/(2*32),1/(2*64)];a2=[1/(2*8),1/(2*16)];

 %solve the exasolution
 Uexasolution=zeros(2,numberx); Vnumsolution=zeros(3,numberx);
for k=1:numberx
    Uexasolution(1,k)=sin(pi*Grid(k));
    Uexasolution(2,k)=pi*cos(pi*Grid(k));
end

%% Proceeding
%initial  condition set up
for k=1:numberx-1
    Vcurrent(1,k)=((Grid(k+1)^3-Grid(k)^3)/3-(Grid(k+1)^2-Grid(k)^2)/2)/Deltax(k);
    Vcurrent(2,k)=Grid(k+1)^2-Grid(k)^2-(Grid(k+1)-Grid(k));
    Vcurrent(3,k)=2*Deltax(k)^2;
end

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
               varphig=Vcurrent(1,ie)+B2g*Vcurrent(2,ie)+B3g*Vcurrent(3,ie);
               Vg=Vcurrent(2,ie)/Deltax(ie)+Vcurrent(3,ie)*B2g/Deltax(ie);
               S1g=nu*pi^2*sin(pi*xig)+a*pi*cos(pi*xig);
               S2g=-Vg/Tr;
               F1g=a*varphig-nu*Vg;
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
               varphiL=Vcurrent(1,ieL)+B2L*Vcurrent(2,ieL)+B3L*Vcurrent(3,ieL);
               varphiR=Vcurrent(1,ieR)+B2R*Vcurrent(2,ieR)+B3R*Vcurrent(3,ieR);
               VL=Vcurrent(2,ieL)/Deltax(ieL)+Vcurrent(3,ieL)*B2L/Deltax(ieL);
               VR=Vcurrent(2,ieR)/Deltax(ieR)+Vcurrent(3,ieR)*B2R/Deltax(ieR);
    
               Fn(:,iface)=0.5*([a*varphiL-nu*VL;-varphiL/Tr]+[a*varphiR-nu*VR;-varphiR/Tr])-0.5*A*([varphiR;VR]-[varphiL;VL]);
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
varphiR=Vcurrent(1,ieR)+B2R*Vcurrent(2,ieR)+B3R*Vcurrent(3,ieR);
VR=Vcurrent(2,ieR)/Deltax(ieR)+Vcurrent(3,ieR)*B2R/Deltax(ieR);
Fn(:,iface)=0.5*([-nu*VR;0]+[a*varphiR-nu*VR;-varphiR/Tr])-0.5*A*([varphiR;VR]-[0;VR]);
Rb(3*(ieR-1)+1)=Rb(3*(ieR-1)+1)+Fn(1,iface);
Rb(3*(ieR-1)+2)=Rb(3*(ieR-1)+2)+Fn(1,iface)*B2R+Fn(2,iface)/Deltax(ieR);
Rb(3*(ieR-1)+3)=Rb(3*(ieR-1)+3)+Fn(1,iface)*B3R+Fn(2,iface)*B2R/Deltax(ieR);

%�߽磺��
iface=numberx;
ieL=iface-1;
B2L=1/2;B3L=0.5*B2L^2-1/24;
varphiL=Vcurrent(1,ieL)+B2L*Vcurrent(2,ieL)+B3L*Vcurrent(3,ieL);
VL=Vcurrent(2,ieL)/Deltax(ieL)+Vcurrent(3,ieL)*B2L/Deltax(ieL);
Fn(:,iface)=0.5*([a*varphiL-nu*VL;-varphiL/Tr]+[-nu*VL;0])-0.5*A*([0;VL]-[varphiL;VL]);
Rb(3*(ieL-1)+1)=Rb(3*(ieL-1)+1)-Fn(1,iface);
Rb(3*(ieL-1)+2)=Rb(3*(ieL-1)+2)-Fn(1,iface)*B2L-Fn(2,iface)/Deltax(ieL);
Rb(3*(ieL-1)+3)=Rb(3*(ieL-1)+3)-Fn(1,iface)*B3L-Fn(2,iface)*B2L/Deltax(ieL);

%R��װ
R=Rd+Rb;

    for k=1:numberx-1
        Mtau=[Deltax(k),0,0;0,Deltax(k)/12+1/Deltax(k),0;0,0,1/720*Deltax(k)+1/(12*Deltax(k))];
        Vnext(:,k)=Vcurrent(:,k)+Mtau\R(:,k)*deltatau(k);
    end



        if max(max(abs(Vnext-Vcurrent)))<tol &&itau>=min(deltatau)
            break;
        end
        
        Vcurrent=Vnext;
        Rd=zeros(3,Unit);
        Rb=zeros(3,Unit);
        
end
    
Vnumsolution=Vcurrent;

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
lgd=legend([H1,H2],'DG(P0P2)+DG(P1)','������');
lgd.FontSize=12;
xlabel('λ��x','fontsize',14)
ylabel('��ֵU','fontsize',14)
 title('DG(P0P2)+DG(P1) Explicit Euler��ֵ���������(U)','fontsize',16)
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
lgd=legend([H1,H2],'DG(P0P2)+DG(P1)','������');
lgd.FontSize=12;
xlabel('λ��x','fontsize',14)
ylabel('��ֵU','fontsize',14)
title('DG(P0P2)+DG(P1) Explicit Euler��ֵ���������(Ux)','fontsize',16)
hold off

%���㾫��
[Acc(1,1),Acc(2,1)]=Accuracy(8);
[Acc(1,2),Acc(2,2)]=Accuracy(16);
[Acc(1,3),Acc(2,3)]=Accuracy(32);
[Acc(1,4),Acc(2,4)]=Accuracy(64);

for k=1:3
accuracyU(k)=(log10(Acc(1,k+1))-log10(Acc(1,k)))./(log10(a1(1,k+1))-log10(a1(1,k)));
end
for k=1:3
accuracyUx(k)=(log10(Acc(2,k+1))-log10(Acc(2,k)))./(log10(a1(1,k+1))-log10(a1(1,k)));
end

%U����
figure
hold on
plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5)
H1=plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5);

H2=plot(log10(a2),1*log10(a2),'--','linewidth',1.5);
 plot(log10(a2),2*log10(a2),'--','linewidth',1.5)
H3=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
lgd=legend([H1,H2,H3],'DG(P0P2)+DG(P1)','Slope=1','Slope=2');
lgd.FontSize=12;
xlabel('Log(1/DOF)','fontsize',14)
ylabel('Log(episilo)','fontsize',14)
title('DG(P0P2)+DG(P1)���ȷ���(U)','fontsize',16)

%Ux����
figure
hold on
plot(log10(a1),log10(Acc(2,:)),'-c*','linewidth',1.5)
H1=plot(log10(a1),log10(Acc(2,:)),'-c*','linewidth',1.5);

H2=plot(log10(a2),1*log10(a2),'--','linewidth',1.5);
 plot(log10(a2),2*log10(a2),'--','linewidth',1.5)
H3=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
lgd=legend([H1,H2,H3],'DG(P0P2)+DG(P1)','Slope=1','Slope=2');
lgd.FontSize=12;
xlabel('Log(1/DOF)','fontsize',14)
ylabel('Log(episilo)','fontsize',14)
title('DG(P0P2)+DG(P1)���ȷ���(Ux)','fontsize',16)




