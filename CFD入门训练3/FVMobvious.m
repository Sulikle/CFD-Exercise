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
UL=0;UR=0;
Unumsolution1=zeros(1,2);
B=zeros(2,numberx-1);
%% solve the question
%initial  condition set up
k=1;
for x=0:deltx:endx-deltx
    if x<0.2&&x+deltx<0.2
        U=0;Ucurrent(1,k)=U;
    elseif x<0.2&&x+deltx>0.2&&x+deltx<=0.3
         U=0+(x+deltx-0.2)/deltx;Ucurrent(1,k)=U;
    elseif x>=0.2&&x<=0.3&&x+deltx>=0.2&&x+deltx<=0.3
        U=1;Ucurrent(1,k)=U;
    elseif x>=0.2&&x<=0.3&&x+deltx>0.3&&x+deltx<=0.4
        U=(0.3-x)/deltx+(0.5*(x+deltx-0.3)^4-(x+deltx-0.3)^3+x+deltx-0.3)/deltx;Ucurrent(1,k)=U;
    elseif x>0.3&&x<=0.4&&x+deltx>0.3&&x+deltx<=0.4
        U=(0.5*((x+deltx-0.3)^4-(x-0.3)^4)-((x+deltx-0.3)^3-(x-0.3)^3)+deltx)/deltx;Ucurrent(1,k)=U;
    elseif x>0.3&&x<=0.4&&x+deltx>0.4
         U=(0.5*((0.4-0.3)^4-(x-0.3)^4)-((0.4-0.3)^3-(x-0.3)^3)+0.4-x)/deltx;Ucurrent(1,k)=U;
    elseif x>0.4
         U=0;Ucurrent(1,k)=U;
    end
    k=k+1;
end

%solve the numsolution
for n=deltt:deltt:endt
    for k=2:numberx-1
       Unext(1,k)=Ucurrent(1,k)+CFL*(Ucurrent(1,k-1)-Ucurrent(1,k));
    end
    Unext(1,1)=UL;Ucurrent=Unext;
end
Unumsolution=Ucurrent;


k=1;
for x=0:deltx:endx
    if x-endt<0.2
        U=0;Uexasolution(1,k)=U;
    elseif x-endt>0.2&&x-endt<=0.3
         U=1;Uexasolution(1,k)=U;
    elseif x-endt>0.3&&x-endt<=0.4
        U=2*(x-endt-0.3)^3-3*(x-endt-0.3)^2+1;Uexasolution(1,k)=U;
    elseif x-endt>0.4
         U=0;Uexasolution(1,k)=U;
    end
    k=k+1;
end
%% post-processing
%calculate the exact value
 figure
 hold on
 x=0*deltx:deltx:1*deltx;
 Unumsolution1(1,1)=Unumsolution(1,1);Unumsolution1(1,2)=Unumsolution(1,1);
 plot(x,Unumsolution1,'-r.');hold on
 H1=plot(x,Unumsolution1,'-r.');hold on
for i=2:numberx-1
    x=(i-1)*deltx:deltx:i*deltx;
    Unumsolution1(1,1)=Unumsolution(1,i);Unumsolution1(1,2)=Unumsolution(1,i);
    plot(x,Unumsolution1,'-r.')
end
y=0:deltx:endx;
plot(y,Uexasolution(1,:),'-b*')
H2=plot(y,Uexasolution(1,:),'-b*');hold on
legend('数值解');hold on
lgd=legend([H1,H2],'数值解','解析解');
lgd.FontSize=12;
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
title('采取FVM方法计算一维线性波动方程(CFL=1)','fontsize',16)
hold off
%calculate the variance
for i=1:numberx-1
    B(1,i)=Unumsolution(1,i)-Uexasolution(1,i);
    B(2,i)=Unumsolution(1,i)-Uexasolution(1,i+1);
end

Var=var(B(1,:))+var(B(2,:))

