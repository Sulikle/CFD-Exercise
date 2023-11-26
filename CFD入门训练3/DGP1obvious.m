clc
clear all
close all
%% Pre-processing
deltx=0.01;CFL=1;deltt=CFL*deltx;
endx=1;endt=0.35;
numberx=endx/deltx+1;
Ucurrent=zeros(2,numberx-1);
Unext=zeros(2,numberx-1);
Unumsolution=zeros(2,numberx-1);
Uexasolution=zeros(1,numberx);
UL=0;UR=0;
M=[deltx,0;0,deltx^3/12];
R=zeros(2,1);
B=zeros(2,numberx-1);
Unumsolution1=zeros(1,2);
%% solve the question
%initial  condition set up
%Uc
k=1;
for x=0:deltx:endx-deltx
    if x+deltx/2<0.2
        U=0;Ucurrent(1,k)=U;
    elseif x+deltx/2>0.2&&x+deltx/2<=0.3
         U=1;Ucurrent(1,k)=U;
    elseif x+deltx/2>0.3&&x+deltx/2<=0.4
        U=2*(x+deltx/2-0.3)^3-3*(x+deltx/2-0.3)^2+1;Ucurrent(1,k)=U;
    elseif x+deltx/2>0.4
         U=0;Ucurrent(1,k)=U;
    end
    k=k+1;
end
%Uxc
k=1;
for x=0:deltx:endx-deltx
    if x+deltx/2<0.2
        U=0;Ucurrent(2,k)=U;
    elseif x+deltx/2>0.2&&x+deltx/2<=0.3
         U=0;Ucurrent(2,k)=U;
    elseif x+deltx/2>0.3&&x+deltx/2<=0.4
        U=6*(x+deltx/2-0.3)^2-6*(x+deltx/2-0.3);Ucurrent(2,k)=U;
    elseif x+deltx/2>0.4
         U=0;Ucurrent(2,k)=U;
    end
    k=k+1;
end
%solve the numsolution
for n=deltt:deltt:endt
    for k=2:numberx-1
        f1=Ucurrent(1,k-1)+Ucurrent(2,k-1)*deltx/2;
        f2=Ucurrent(1,k)+Ucurrent(2,k)*deltx/2;
        R(1,1)=f1-f2;
        R(2,1)=-deltx/2*(f1+f2)+deltx*Ucurrent(1,k);
        Unext(:,k)=Ucurrent(:,k)+M\R*deltt;
    end
    Unext(1,1)=UL;Unext(2,1)=UL;Ucurrent=Unext;
end
for i=1:numberx-1
Unumsolution(1,i)=Ucurrent(1,i)+Ucurrent(2,i)*(-deltx/2);
Unumsolution(2,i)=Ucurrent(1,i)+Ucurrent(2,i)*(deltx/2);
end

%solve the exasolution
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
figure 
hold on
x=0*deltx:deltx:1*deltx;
Unumsolution1(1,1)=Unumsolution(1,1);Unumsolution1(1,2)=Unumsolution(2,1);
plot(x,Unumsolution1,'-r.');hold on
H1=plot(x,Unumsolution1,'-r.');hold on
for i=2:numberx-1
    x=(i-1)*deltx:deltx:i*deltx;
    Unumsolution1(1,1)=Unumsolution(1,i);Unumsolution1(1,2)=Unumsolution(2,i);
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
title('采取DGM(P1)方法计算一维线性波动方程(CFL=1)','fontsize',16)
hold off
% x1=0:deltx:endx-deltx;x2=deltx:deltx:endx;y=0:deltx:endx;
% figure
% scatter(x1,Unumsolution(1,:),50,'b','p','filled')
% hold on
% scatter(x2,Unumsolution(2,:))
% plot(y,Uexasolution(1,:),'-r*')
% lgd=legend('数值解左端点','数值解左端点','解析解');
% lgd.FontSize=12;
% xlabel('位置x','fontsize',14)
% ylabel('数值U','fontsize',14)
% title('采取DGM(P1)方法计算一维线性波动方程(CFL=1)','fontsize',16)
% hold off

for i=1:numberx-1
    B(1,i)=Uexasolution(1,i)-Unumsolution(1,i);
end
for i=1:numberx-1
    B(2,i)=Uexasolution(1,i+1)-Unumsolution(2,i);
end
Var=var(B(1,:))+var(B(2,:))