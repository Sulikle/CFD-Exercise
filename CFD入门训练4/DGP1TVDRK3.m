clc
clear all
close all
%% Pre-processing
deltx=0.01;CFL=1;deltt=CFL*deltx;
endx=1;endt=2;
numberx=endx/deltx+1;
Ucurrent=zeros(2,numberx-1);
Unext=zeros(2,numberx-1);
Unumsolution=zeros(2,numberx-1);
Uexasolution=zeros(1,numberx);
M=[deltx,0;0,deltx/12];
R=zeros(2,1);R1=zeros(2,1);
Unumsolution1=zeros(1,2);
Unk=zeros(4,numberx-1);%store the intermediate quantity of TVDRK3
Uhold=zeros(2,numberx-1);
afa=[0,1/4,2/3];beta=[1,1/4,2/3];gama=[1,3/4,1/3];
%% solve the question
%initial  condition set up(dimensionless) 
%Uc
k=1;
for x=0:deltx:endx-deltx
  Ucurrent(1,k)=sin(2*pi*(x+deltx/2));
    k=k+1;
end
%Uxc
k=1;
for x=0:deltx:endx-deltx
   Ucurrent(2,k)=2*pi*cos(2*pi*(x+deltx/2))*deltx;
    k=k+1;
end
Unk([1,2],:)=Ucurrent;Uhold=Ucurrent;

%solve the numsolution TVDRK3
for n=deltt:deltt:endt
    for istage=1:3
        R1(1,1)=Ucurrent(1,numberx-1)+Ucurrent(2,numberx-1)/2-Ucurrent(1,1)-Ucurrent(2,1)/2;
        R1(2,1)=-0.5*(Ucurrent(1,numberx-1)+Ucurrent(2,numberx-1)/2+Ucurrent(1,1)+Ucurrent(2,1)/2)+Ucurrent(1,1);
        Unk([3,4],1)=gama(istage)*Unk([1,2],1)+afa(istage)*Uhold([1,2],1)+M\R1*beta(istage)*deltt;
    for k=2:numberx-1
        f2=Ucurrent(1,k)+Ucurrent(2,k)/2;
        f1=Ucurrent(1,k-1)+Ucurrent(2,k-1)/2;
        R(1,1)=f1-f2;
        R(2,1)=-0.5*(f1+f2)+Ucurrent(1,k);
       Unk([3,4],k)=gama(istage)*Unk([1,2],k)+afa(istage)*Uhold([1,2],k)+M\R*beta(istage)*deltt;
    end
    Uhold([1,2],:)=Unk([3,4],:);
    end
    Unext([1,2],:)=Uhold([1,2],:);Ucurrent=Unext;Unk([1,2],:)=Unext([1,2],:);Uhold([1,2],:)=Unext([1,2],:);
    
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
title('DG(P1)TVDRK3(CFL=1，ENDT=2)','fontsize',16)
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