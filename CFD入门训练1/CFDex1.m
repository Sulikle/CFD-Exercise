clc
clear all
close all

%% parameter set up
afa=0.06;
deltx=0.02;
deltt=0.002;
T0=50;
endt=10;endx=1;
numberx=endx/deltx+1;numbert=endt/deltt+1;
A=zeros(numbert,numberx);


%% solve the question
%initial  condition set up
k=1;
for x=0:deltx:endx
            T=T0*sin(pi*x);A(1,k)=T;k=k+1;
end

if A(1,k-1)~=0
    A(1,k-1)=0;
end

%solve 
for n=2:1:numbert
  for i=2:1:numberx-1
      Tin=A(n-1,i)+afa*deltt/(deltx)^2*(A(n-1,i+1)-2*A(n-1,i)+A(n-1,i-1));A(n,i)=Tin;%calculate  inner value
  end
  A(n,1)=0;A(n,numberx)=0;%boundary condition set up
end

%% post-processing
%calculate the exact value
B1=A(numbert,:);
B2=zeros(1,numberx);
p=1;
for x=0:deltx:endx
    T=T0*sin(pi*x)*exp((-afa*(pi)^2)*endt);B2(1,p)=T;p=p+1;
end
%calculate the variance
B=B2-B1;
Var=var(B)

%figure
x=0:deltx:endx;
scatter(x,B1)
hold on
plot(x,B2,'-r*')
legend('数值解','解析解')
xlabel('位置x','fontsize',14)
ylabel('温度T','fontsize',14)
title('t=10时，T的数值解与解析解','fontsize',16)
    













