clc
clear all
close all
%% Pre-processing
deltx=0.01;CFL=1;deltt=CFL*deltx;
endx=1;endt=0.35;
numberx=endx/deltx+1;
Ucurrent=zeros(1,numberx-1);
Unext=zeros(1,numberx-1);
Unumsolution=zeros(1,numberx-1);
Uexasolution=zeros(1,numberx);
Uxr=zeros(1,numberx-1);
Unumsolution1=zeros(1,2);
Unk=zeros(2,numberx-1);%store the intermediate quantity of TVDRK3
Uhold=zeros(1,numberx-1);
afa=[0,1/4,2/3];beta=[1,1/4,2/3];gama=[1,3/4,1/3];

%% solve the question
%initial  condition set up
k=1;
for x=0:deltx:endx-deltx
    Ucurrent(1,k)=0.5*(sin(2*pi*x)+sin(2*pi*(x+deltx)));
    k=k+1;
end
for k=2:numberx-2
    Uxr(1,k)=(Ucurrent(1,k+1)-Ucurrent(1,k-1))/(2*deltx);
end
Uxr(1,1)=(Ucurrent(1,2)-Ucurrent(1,numberx-1))/(2*deltx);Uxr(1,numberx-1)=(Ucurrent(1,1)-Ucurrent(1,numberx-2))/(2*deltx);
Unk(1,:)=Ucurrent(1,:);Uhold(1,:)=Ucurrent(1,:);

%solve the numsolution TVDRK3
for n=deltt:deltt:endt
    for istage=1:3
        Unk(2,1)=gama(istage)*Unk(1,1)+afa(istage)*Uhold(1,1)+beta(istage)*deltt*(Ucurrent(numberx-1)+0.5*Uxr(numberx-1)*deltx-Ucurrent(1)+0.5*Uxr(1)*deltx)/deltx;
    for k=2:numberx-1
        f2=Ucurrent(k)+0.5*Uxr(k)*deltx;
        f1=Ucurrent(k-1)+0.5*Uxr(k-1)*deltx;
        R=f1-f2;
       Unk(2,k)=gama(istage)*Unk(1,k)+afa(istage)*Uhold(1,k)+beta(istage)*deltt*R/deltx;
    end
    Uhold(1,:)=Unk(2,:);
    end
    Unext(1,:)=Uhold(1,:);Ucurrent=Unext;Unk(1,:)=Unext(1,:);Uhold(1,:)=Unext(1,:);
    for k=2:numberx-2
    Uxr(1,k)=(Ucurrent(1,k+1)-Ucurrent(1,k-1))/(2*deltx);
    end
    Uxr(1,1)=(Ucurrent(1,2)-Ucurrent(1,numberx-1))/(2*deltx);Uxr(1,numberx-1)=(Ucurrent(1,1)-Ucurrent(1,numberx-2))/(2*deltx);
end
Unumsolution=Ucurrent;

%solve the exasolution
k=1;
for x=0:deltx:endx
    Uexasolution(1,k)=sin(2*pi*(x-endt));
    k=k+1;
end

%% post-processing
%calculate the exact value
 figure
 hold on
 x=0*deltx:deltx:1*deltx;
 Unumsolution1(1,1)=Unumsolution(1,1);Unumsolution1(1,2)=Unumsolution(1,1);
 plot(x,Unumsolution1,'-ro');hold on
 H1=plot(x,Unumsolution1,'-ro');hold on
for i=2:numberx-1
    x=(i-1)*deltx:deltx:i*deltx;
    Unumsolution1(1,1)=Unumsolution(1,i);Unumsolution1(1,2)=Unumsolution(1,i);
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
title('rDG(P0P1)TVDRK3(CFL=1，ENDT=0.35)','fontsize',16)
hold off

%calculate the accuracy of space
I=0;t=[-1/sqrt(5),0,1/sqrt(5)];W=[5/9,8/9,5/9];
k=1;%determine the correctness of the program
for x=0:deltx:endx-deltx
   for i=1:3
       xi=deltx/2*t(i)+0.5*(2*x+deltx);
       for m=1:numberx-1
           if xi>(m-1)*deltx&&xi<m*deltx
              fi=(sin(2*pi*(xi-endt))-Unumsolution(1,m))^2;k=k+1;
           end
       end
       I=I+W(i)*fi;
   end
end
I=I*0.5*deltx;
I=sqrt(I)