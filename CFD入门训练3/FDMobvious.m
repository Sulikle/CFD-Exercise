clc
clear all
close all
%% Pre-processing
deltx=0.01;CFL=1;deltt=CFL*deltx;
endx=1;endt=0.35;
numberx=endx/deltx+1;
Ucurrent=zeros(1,numberx);
Unext=zeros(1,numberx);
Unumsolution=zeros(1,numberx);
Uexasolution=zeros(1,numberx);
UL=0;UR=0;
%% solve the question
%initial  condition set up
Ucurrent(1,1)=UL;Ucurrent(1,numberx)=UR;k=2;
for x=deltx:deltx:endx-deltx
    if x<0.2
        U=0;Ucurrent(1,k)=U;
    elseif x>=0.2&&x<=0.3
        U=1;Ucurrent(1,k)=U;
    elseif x>0.3&&x<=0.4
        U=2*(x-0.3)^3-3*(x-0.3)^2+1;Ucurrent(1,k)=U;
    elseif x>0.4
         U=0;Ucurrent(1,k)=U;
    end
    k=k+1;
end

%solve the numsolution
for n=deltt:deltt:endt
    for k=2:numberx-1
       Unext(1,k)=CFL*(Ucurrent(1,k-1)-Ucurrent(1,k))+Ucurrent(1,k);
    end
    Unext(1,1)=UL;Unext(1,numberx)=UR;Ucurrent=Unext;
end
Unumsolution=Ucurrent;

%solve the exasolution
Uexasolution(1,1)=UL;Uexasolution(1,numberx)=UR;k=2;
for x=deltx:deltx:endx-deltx
    if x-endt<0.2
        U=0;Uexasolution(1,k)=U;
    elseif x-endt>=0.2&&x-endt<=0.3
        U=1;Uexasolution(1,k)=U;
    elseif x-endt>0.3&&x-endt<=0.4
        U=2*(x-endt-0.3)^3-3*(x-endt-0.3)^2+1;Uexasolution(1,k)=U;
    elseif  x-endt>0.4
         U=0;Uexasolution(1,k)=U;
    end
    k=k+1;
end

%% post-processing
%calculate the exact value
x=0:deltx:endx;
figure
scatter(x,Unumsolution)
hold on
plot(x,Uexasolution,'-r*')
legend('数值解','解析解')
xlabel('位置x','fontsize',14)
ylabel('数值U','fontsize',14)
title('采取FDM方法计算一维线性波动方程(CFL=1)','fontsize',16)

hold off
%calculate the variance
B=Uexasolution-Unumsolution;
Var=var(B);


