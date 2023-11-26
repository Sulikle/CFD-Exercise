clc
clear all
close all
%% Pre-processing
afa=0.006;
T0=50;TL=0;TR=0;N=4;% if you want to change the numeric value and number of deltt, you just need to change N and deltt
deltx=0.02;deltt=[0.002,0.004,0.008,0.02];
endx=1;endt=10;
numberx=endx/deltx+1;
Tcurrent=zeros(1,numberx-2);
Tnext=zeros(1,numberx-2)';
Tnumsolution=zeros(1,numberx);
Texasolution=zeros(1,numberx);
Var=zeros(1,N);
for time=1:1:N
numbert=endt/deltt(time)+1;
sita=afa*deltt(time)/deltx^2;
A1=sparse(1:numberx-2,1:numberx-2,1+2*sita,numberx-2,numberx-2);
A2=sparse(1:numberx-3,2:numberx-2,-sita,numberx-2,numberx-2);
A=A1+A2+A2';
%% solve the question
%initial  condition set up
k=1;
for x=deltx:deltx:endx-deltx
            T=T0*sin(pi*x);Tcurrent(1,k)=T;k=k+1;
end
k=k-1;

%solve(左除) 
for n=2:1:numbert
    f=Tcurrent';f(1,1)=Tcurrent(1,1)+sita*TL;f(numberx-2,1)=Tcurrent(1,numberx-2)+sita*TR;
    Tnext=A\f;Tcurrent=Tnext';
end

for j=2:numberx-1
    Tnumsolution(1,j)=Tcurrent(1,j-1);
end
Tnumsolution(1,1)=TL;Tnumsolution(1,numberx)=TR;

%% post-processing
%calculate the exact value
p=2;
for x=deltx:deltx:endx-deltx
    T=T0*sin(pi*x)*exp((-afa*(pi)^2)*endt);Texasolution(1,p)=T;p=p+1;
end
Texasolution(1,1)=TL;Texasolution(1,p)=TR;

%figure
x=0:deltx:endx;
figure
scatter(x,Tnumsolution)
hold on
plot(x,Texasolution,'-r*')
legend('数值解','解析解')
xlabel('位置x','fontsize',14)
ylabel('温度T','fontsize',14)
titleName=strcat('t=10,delt为第',num2str(time),'个值时，T的数值解与解析解');
title(titleName)
hold off
%calculate the variance
B=Texasolution-Tnumsolution;
Var(time)=var(B);
end
figure
time=1:N;
plot(time,Var,'-r*')
xlabel('第i个deltt','fontsize',14)
ylabel('解析解与数值解的方差','fontsize',14)
title('方差与deltt取值的关系','fontsize',16)