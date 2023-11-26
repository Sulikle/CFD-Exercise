clc
close all
clear all
%% Pre-proceeding
n=20;
D=3*eye(n);
L=sparse(2:n,1:n-1,-1/2,n,n)+sparse(3:n,1:n-2,-1/4,n,n);
U=L';
A=D+L+U;
b=100*rand(n,1);b1=b;
X0=rand(n,1);XF=zeros(n,1);XB=zeros(n,1);
tol=10^(-5);endtimes=100;epsilon=10^(-12);
ratioRF=[0:endtimes;zeros(1,endtimes+1)];
ratioRB=[0:endtimes;zeros(1,endtimes+1)];
XF=X0;XB=X0;
resnorm0=sqrt(sum((A*X0-b).^2))+epsilon;
resnorm=sqrt(sum((A*XF-b).^2));
ratio=resnorm/resnorm0;ratioRF(2,1)=ratio;ratioRB(2,1)=ratio;

%% Proceeding
%Forward sweep
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


%Backward sweep
for times=2:endtimes
    %calculate X
    for i=n:-1:1
        for j=n:-1:i+1
            b(i)=b(i)-A(i,j)*XB(j);
        end
        for j=i-1:-1:1
            b(i)=b(i)-A(i,j)*XB(j);
        end
        XB(i)=b(i)/A(i,i);
    end
    b=b1;
    resnorm=sqrt(sum((A*XB-b).^2));
    ratio=resnorm/resnorm0;
    if ratio<tol
        break
    end
    ratioRB(2,times)=ratio;
end
ind=find(ratioRB(2,:),1,'last');
ratioRB(:,ind+1:endtimes+1)=[];
subplot(1,2,2)
plot(ratioRB(1,:),log10(ratioRB(2,:)),'-*','linewidth',1.5)
xlabel('Iteration times','fontsize',14)
ylabel('Log(ratio)','fontsize',14)
title('G-S Backward sweep','fontsize',16)
hold on
str1=num2str(ratioRB(2,:)');text(ratioRB(1,:),log10(ratioRB(2,:)),str1,'linewidth',1.5);
line([ratioRB(1,ind),ratioRB(1,ind)], [-5,1], 'color', 'b','linewidth',1.5);
str2=num2str(ind-1);
text(ind-1,0,str2,'linewidth',1.5);
xlim([0,ind])


