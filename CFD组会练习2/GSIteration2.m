clc
close all
clear all
%% Pre-proceeding
n=20;
omega=1;%松弛因子
%构建矩阵
D=3*eye(n);
L=sparse(2:n,1:n-1,-1/2,n,n)+sparse(3:n,1:n-2,-1/4,n,n);
U=L';
A=D+L+U;
b=[1:n]';b1=b;
X0=zeros(n,1);XF=zeros(n,1);Xsgs=zeros(n,1);Xold=zeros(n,1);
%终止条件等
tol=10^(-5);endtimes=100;epsilon=10^(-12);
ratioRF=[0:endtimes;zeros(1,endtimes+1)];
ratioRsgs=[0:endtimes;zeros(1,endtimes+1)];
ratioRSORFGS=[0:endtimes;zeros(1,endtimes+1)];
ratioRSORSGS=[0:endtimes;zeros(1,endtimes+1)];
XF=X0;Xsgs=X0;Xsorfgs=X0;Xsorsgs=X0;
resnorm0=sqrt(sum((A*X0-b).^2))+epsilon;
resnorm=sqrt(sum((A*XF-b).^2));
ratio=resnorm/resnorm0;ratioRF(2,1)=ratio;ratioRsgs(2,1)=ratio;ratioRSORFGS(2,1)=ratio;ratioRSORSGS(2,1)=ratio;

%% Proceeding
%G-S-Forward sweep
for times=2:endtimes
    %calculate X
    for i=1:n
        for j=1:i-1
            b(i)=b(i)-A(i,j)*XF(j);
        end
        for j=i+1:n
            b(i)=b(i)-A(i,j)*XF(j);
        end
        XF(i)=b(i)/A(i,i);
    end
     b=b1;
    resnorm=sqrt(sum((A*XF-b).^2));
    ratio=resnorm/resnorm0;
    if ratio<tol
        break
    end
    ratioRF(2,times)=ratio;
end
ind=find(ratioRF(2,:),1,'last');
ratioRF(:,ind+1:endtimes+1)=[];
subplot(1,2,1)
plot(ratioRF(1,:),log10(ratioRF(2,:)),'-*','linewidth',1.5)
xlabel('Iteration times','fontsize',14)
ylabel('Log(ratio)','fontsize',14)
title('G-S Forward sweep','fontsize',16)
hold on
str1=num2str(ratioRF(2,:)');text(ratioRF(1,:),log10(ratioRF(2,:)),str1,'linewidth',1.5);
line([ratioRF(1,ind),ratioRF(1,ind)], [-5,1], 'color', 'b','linewidth',1.5);
str2=num2str(ind-1);
text(ind-1,0,str2,'linewidth',1.5);
xlim([0,ind])

%SOR-FGS 
R=b1-A*X0;
for times=2:endtimes
    %calculate deltaX
    Xold=Xsorfgs;
    ie=1;
    Xsorfgs(ie)=R(ie)/(D(ie,ie)/omega);
    ie=2;
    R(ie)=R(ie)-Xsorfgs(ie-1)*(L(ie,ie-1)/omega);
    Xsorfgs(ie)=R(ie)/(D(ie,ie)/omega);
    for ie=3:n
        R(ie)=R(ie)-Xsorfgs(ie-1)*(L(ie,ie-1)/omega)-Xsorfgs(ie-2)*(L(ie,ie-2)/omega);
        Xsorfgs(ie)=R(ie)/(D(ie,ie)/omega);
    end
    Xsorfgs=Xsorfgs+Xold;
    resnorm=sqrt(sum((A*Xsorfgs-b1).^2));
    ratio=resnorm/resnorm0;
    if ratio<tol
        break
    end
    R=b1-A*Xsorfgs;
    ratioRSORFGS(2,times)=ratio;
end
ind=find(ratioRSORFGS(2,:),1,'last');
ratioRSORFGS(:,ind+1:endtimes+1)=[];
subplot(1,2,2)
plot(ratioRSORFGS(1,:),log10(ratioRSORFGS(2,:)),'-*','linewidth',1.5)
xlabel('Iteration times','fontsize',14)
ylabel('Log(ratio)','fontsize',14)
title('G-S Forward sweep(SOR)','fontsize',16)
hold on
str1=num2str(ratioRSORFGS(2,:)');text(ratioRSORFGS(1,:),log10(ratioRSORFGS(2,:)),str1,'linewidth',1.5);
line([ratioRSORFGS(1,ind),ratioRSORFGS(1,ind)], [-5,1], 'color', 'b','linewidth',1.5);
str2=num2str(ind-1);
text(ind-1,0,str2,'linewidth',1.5);
xlim([0,ind])

%SGS
for times=2:endtimes
    %calculate X Forward sweep
    for i=1:n
        for j=1:i-1
            b(i)=b(i)-A(i,j)*Xsgs(j);
        end
        for j=i+1:n
            b(i)=b(i)-A(i,j)*Xsgs(j);
        end
        Xsgs(i)=b(i)/A(i,i);
    end
     b=b1;
     %Backward sweep
    for i=n:-1:1
        for j=n:-1:i+1
            b(i)=b(i)-A(i,j)*Xsgs(j);
        end
        for j=i-1:-1:1
            b(i)=b(i)-A(i,j)*Xsgs(j);
        end
        Xsgs(i)=b(i)/A(i,i);
    end
    b=b1;
    resnorm=sqrt(sum((A*Xsgs-b).^2));
    ratio=resnorm/resnorm0;
    if ratio<tol
        break
    end
    ratioRsgs(2,times)=ratio;
end
ind=find(ratioRsgs(2,:),1,'last');
ratioRsgs(:,ind+1:endtimes+1)=[];
figure
subplot(1,2,1)
plot(ratioRsgs(1,:),log10(ratioRsgs(2,:)),'-*','linewidth',1.5)
xlabel('Iteration times','fontsize',14)
ylabel('Log(ratio)','fontsize',14)
title('S-G-S','fontsize',16)
hold on
str1=num2str(ratioRsgs(2,:)');text(ratioRsgs(1,:),log10(ratioRsgs(2,:)),str1,'linewidth',1.5);
line([ratioRsgs(1,ind),ratioRsgs(1,ind)], [-5,1], 'color', 'b','linewidth',1.5);
str2=num2str(ind-1);
text(ind-1,0,str2,'linewidth',1.5);
xlim([0,ind])
hold off


%SOR-SGS
R=omega*(2-omega)*(b1-A*X0);
for times=2:endtimes
    %calculate deltaX   
    %Forward sweep
    Xold=Xsorsgs;
    ie=1;
    Xsorsgs(ie)=R(ie)/D(ie,ie);
    ie=2;
    R(ie)=R(ie)-Xsorsgs(ie-1)*(L(ie,ie-1)*omega);
    Xsorsgs(ie)=R(ie)/D(ie,ie);
    for ie=3:n
        R(ie)=R(ie)-Xsorsgs(ie-1)*(L(ie,ie-1)*omega)-Xsorsgs(ie-2)*(L(ie,ie-2)*omega);
        Xsorsgs(ie)=R(ie)/D(ie,ie);
    end
    %Backward sweep
    R=D*Xsorsgs;
    ie=n;
    Xsorsgs(ie)=R(ie)/D(ie,ie);
    ie=n-1;
    R(ie)=R(ie)-Xsorsgs(ie+1)*(U(ie,ie+1)*omega);
    Xsorsgs(ie)=R(ie)/D(ie,ie);
    for ie=n-2:-1:1
        R(ie)=R(ie)-Xsorsgs(ie+1)*(U(ie,ie+1)*omega)-Xsorsgs(ie+2)*(U(ie,ie+2)*omega);
        Xsorsgs(ie)=R(ie)/D(ie,ie);
    end
    %此时Xsorsgs里面存储的是deltaX
    Xsorsgs=Xsorsgs+Xold;
    resnorm=sqrt(sum((A*Xsorsgs-b1).^2));
    ratio=resnorm/resnorm0;
    if ratio<tol
        break
    end
    R=omega*(2-omega)*(b1-A*Xsorsgs);
    ratioRSORSGS(2,times)=ratio;
end
ind=find(ratioRSORSGS(2,:),1,'last');
ratioRSORSGS(:,ind+1:endtimes+1)=[];
subplot(1,2,2)
plot(ratioRSORSGS(1,:),log10(ratioRSORSGS(2,:)),'-*','linewidth',1.5)
xlabel('Iteration times','fontsize',14)
ylabel('Log(ratio)','fontsize',14)
title('S-G-S(SOR)','fontsize',16)
hold on
str1=num2str(ratioRSORSGS(2,:)');text(ratioRSORSGS(1,:),log10(ratioRSORSGS(2,:)),str1,'linewidth',1.5);
line([ratioRSORSGS(1,ind),ratioRSORSGS(1,ind)], [-5,1], 'color', 'b','linewidth',1.5);
str2=num2str(ind-1);
text(ind-1,0,str2,'linewidth',1.5);
xlim([0,ind])


