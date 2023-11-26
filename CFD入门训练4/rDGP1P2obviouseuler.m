clc
clear all
close all
%% Pre-processing
deltx=0.01;CFL=1;deltt=CFL*deltx;
endx=1;endt=2;
numberx=endx/deltx+1;
Ucurrent=zeros(2,numberx-1);
Unext=zeros(2,numberx-1);
Uxx=zeros(1,numberx-1);
Unumsolution=zeros(2,numberx-1);
Uexasolution=zeros(1,numberx);
M=[deltx,0;0,deltx/12];
R=zeros(2,1);R1=zeros(2,1);
Unumsolution1=zeros(1,2);
A=[1/12;1/2;1/12;-1/2];
%% solve the question
%initial  condition set up(dimensionless) 
%Uc
k=1;
for x=0:deltx:endx-deltx
  Ucurrent(1,k)=(cos(2*pi*x)-cos(2*pi*(x+deltx)))/(2*pi*deltx);
    k=k+1;
end
%Uxc
k=1;
for x=0:deltx:endx-deltx
   Ucurrent(2,k)=2*pi*cos(2*pi*(x+deltx/2))*deltx;
    k=k+1;
end

%Uxxc
for k=2:numberx-2
    b=[Ucurrent(1,k+1)-0.5*Ucurrent(2,k+1)-Ucurrent(1,k)-0.5*Ucurrent(2,k);
        Ucurrent(2,k+1)-Ucurrent(2,k);
        Ucurrent(1,k-1)+0.5*Ucurrent(2,k-1)-Ucurrent(1,k)+0.5*Ucurrent(2,k);
        Ucurrent(2,k-1)-Ucurrent(2,k)];
    Uxx(k)=A\b;
end
b=[Ucurrent(1,2)-0.5*Ucurrent(2,2)-Ucurrent(1,1)-0.5*Ucurrent(2,1);
        Ucurrent(2,2)-Ucurrent(2,1);
        Ucurrent(1,numberx-1)+0.5*Ucurrent(2,numberx-1)-Ucurrent(1,1)+0.5*Ucurrent(2,1);
        Ucurrent(2,numberx-1)-Ucurrent(2,1)];
Uxx(1)=A\b;
b=[Ucurrent(1,1)-0.5*Ucurrent(2,1)-Ucurrent(1,numberx-1)-0.5*Ucurrent(2,numberx-1);
        Ucurrent(2,1)-Ucurrent(2,numberx-1);
        Ucurrent(1,numberx-2)+0.5*Ucurrent(2,numberx-2)-Ucurrent(1,numberx-1)+0.5*Ucurrent(2,numberx-1);
        Ucurrent(2,numberx-2)-Ucurrent(2,numberx-1)];
Uxx(numberx-1)=A\b;

%solve the numsolution
for n=deltt:deltt:endt
    for k=2:numberx-1
        f1=Ucurrent(1,k-1)+Ucurrent(2,k-1)/2+1/12*Uxx(k-1);
        f2=Ucurrent(1,k)+Ucurrent(2,k)/2+1/12*Uxx(k);
        R(1,1)=f1-f2;
        R(2,1)=-0.5*(f1+f2)+Ucurrent(1,k);
        Unext(:,k)=Ucurrent(:,k)+M\R*deltt;
    end
    R1(1,1)=Ucurrent(1,numberx-1)+Ucurrent(2,numberx-1)/2+1/12*Uxx(numberx-1)-Ucurrent(1,1)-Ucurrent(2,1)/2-1/12*Uxx(1);
    R1(2,1)=-0.5*(Ucurrent(1,numberx-1)+Ucurrent(2,numberx-1)/2+1/12*Uxx(numberx-1)+Ucurrent(1,1)+Ucurrent(2,1)/2+1/12*Uxx(1))+Ucurrent(1,1);
    Unext(:,1)=Ucurrent(:,1)+M\R1*deltt;
    Ucurrent=Unext;
    for k=2:numberx-2
    b=[Ucurrent(1,k+1)-0.5*Ucurrent(2,k+1)-Ucurrent(1,k)-0.5*Ucurrent(2,k);
        Ucurrent(2,k+1)-Ucurrent(2,k);
        Ucurrent(1,k-1)+0.5*Ucurrent(2,k-1)-Ucurrent(1,k)+0.5*Ucurrent(2,k);
        Ucurrent(2,k-1)-Ucurrent(2,k)];
    Uxx(k)=A\b;
    end
    b=[Ucurrent(1,2)-0.5*Ucurrent(2,2)-Ucurrent(1,1)-0.5*Ucurrent(2,1);
        Ucurrent(2,2)-Ucurrent(2,1);
        Ucurrent(1,numberx-1)+0.5*Ucurrent(2,numberx-1)-Ucurrent(1,1)+0.5*Ucurrent(2,1);
        Ucurrent(2,numberx-1)-Ucurrent(2,1)];
    Uxx(1)=A\b;
    b=[Ucurrent(1,1)-0.5*Ucurrent(2,1)-Ucurrent(1,numberx-1)-0.5*Ucurrent(2,numberx-1);
        Ucurrent(2,1)-Ucurrent(2,numberx-1);
        Ucurrent(1,numberx-2)+0.5*Ucurrent(2,numberx-2)-Ucurrent(1,numberx-1)+0.5*Ucurrent(2,numberx-1);
        Ucurrent(2,numberx-2)-Ucurrent(2,numberx-1)];
    Uxx(numberx-1)=A\b;    
end

for i=1:numberx-1
Unumsolution(1,i)=Ucurrent(1,i)+Ucurrent(2,i)*(-1/2);
Unumsolution(2,i)=Ucurrent(1,i)+Ucurrent(2,i)*(1/2);
end

%solve the exasolution
k=1;
for x=0:deltx:endx
    Uexasolution(1,k)=sin(2*pi*(x-endt));
    k=k+1;
end


%% post-processing
figure 
hold on
x=0*deltx:deltx:1*deltx;
Unumsolution1(1,1)=Unumsolution(1,1);Unumsolution1(1,2)=Unumsolution(2,1);
plot(x,Unumsolution1,'-ro');hold on
H1=plot(x,Unumsolution1,'-ro');hold on
for i=2:numberx-1
    x=(i-1)*deltx:deltx:i*deltx;
    Unumsolution1(1,1)=Unumsolution(1,i);Unumsolution1(1,2)=Unumsolution(2,i);
    plot(x,Unumsolution1,'-ro')
end
y=0:deltx:endx;
plot(y,Uexasolution(1,:),'-b*')
H2=plot(y,Uexasolution(1,:),'-b*');hold on
legend('数值解');hold on
lgd=legend([H1,H2],'数值解','解析解');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
title('rDG(P1P2)显式欧拉(CFL=1，ENDT=2)','fontsize',16)
hold off


%calculate the accuracy of space
I=0;t=[-1/sqrt(5),0,1/sqrt(5)];W=[5/9,8/9,5/9];
for x=0:deltx:endx-deltx
   for i=1:3
       xi=deltx/2*t(i)+0.5*(2*x+deltx);
       for m=1:numberx-1
           if xi>(m-1)*deltx&&xi<m*deltx
              fi=(sin(2*pi*(xi-endt))-(Unumsolution(1,m)+Unumsolution(2,m)/deltx*(xi-((m-1)*deltx+deltx/2))))^2;
           end
       end
       I=I+W(i)*fi;
   end
end
I=I*0.5*deltx;
I=sqrt(I)