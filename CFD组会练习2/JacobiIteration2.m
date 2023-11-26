clc
close all
clear all
%% Pre-proceeding
n=20;
D=3*eye(n);
L=sparse(2:n,1:n-1,-1/2,n,n)+sparse(3:n,1:n-2,-1/4,n,n);
U=L';b=[1:n]';X0=zeros(n,1);Xold=zeros(n,1);Xnew=zeros(n,1);
tol=10^(-5);endtimes=100;epsilon=10^(-12);
ratioR=[0:endtimes;zeros(1,endtimes+1)];
%% Proceeding
Xold=X0;resnorm0=sqrt(sum(((D+L+U)*X0-b).^2))+epsilon;
resnorm=sqrt(sum(((D+L+U)*Xold-b).^2));
ratio=resnorm/resnorm0;ratioR(2,1)=ratio;
for times=2:endtimes
    Xnew=-D\(L+U)*Xold+D\b;
    resnorm=sqrt(sum(((D+L+U)*Xnew-b).^2));
    ratio=resnorm/resnorm0;
    if ratio<tol
        break
    end
    ratioR(2,times)=ratio;
    Xold=Xnew;
end
ind=find(ratioR(2,:),1,'last');
ratioR(:,ind+1:endtimes+1)=[];
plot(ratioR(1,:),log10(ratioR(2,:)),'-*','linewidth',1.5)
xlabel('Iteration times','fontsize',14)
ylabel('Log(ratio)','fontsize',14)
title('Element on the main diagonal is 3','fontsize',16)
hold on
str1=num2str(ratioR(2,:)');text(ratioR(1,:),log10(ratioR(2,:)),str1,'linewidth',1.5);
line([ratioR(1,ind),ratioR(1,ind)], [-5,1], 'color', 'b','linewidth',1.5);
str2=num2str(ind-1);
text(ind-1,0,str2,'linewidth',1.5);
xlim([0,ind])

