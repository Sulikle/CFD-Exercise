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
Y=zeros(1,numberx-2)';
f=zeros(1,numberx-2)';
for time=1:1:N
numbert=endt/deltt(time)+1;
sita=afa*deltt(time)/(2*deltx^2);
B=ones(1,numberx-2)*(-(1+2*sita));A=ones(1,numberx-3)*(sita);C=ones(1,numberx-3)*(sita);U=zeros(1,numberx-2);L=zeros(1,numberx-3);
U(1,1)=-(1+2*sita);
for i=1:numberx-3
    L(1,i)=A(1,i)/U(1,i);
    U(1,i+1)=B(1,i+1)-L(1,i)*C(1,i);
end
%% solve the question
%initial  condition set up
k=1;
for x=deltx:deltx:endx-deltx
            T=T0*sin(pi*x);Tcurrent(1,k)=T;k=k+1;
end
k=k-1;

%solve(Thomas)
for n=2:1:numbert
    for k=2:numberx-3
        f(k,1)=-Tcurrent(1,k)-sita*(Tcurrent(1,k+1)-2*Tcurrent(1,k)+Tcurrent(1,k-1));
    end
    f(1,1)=-Tcurrent(1,1)-sita*(Tcurrent(1,2)-2*Tcurrent(1,1)+TL)-sita*TL;
    f(numberx-2,1)=-Tcurrent(1,numberx-2)-sita*(TR-2*Tcurrent(1,numberx-2)+Tcurrent(1,numberx-3))-sita*TR;
    Y(1,1)=f(1,1);
    for i=2:numberx-2
        Y(i,1)=f(i,1)-L(1,i-1)*Y(i-1,1);
    end
    Tnext(numberx-2,1)=Y(numberx-2,1)/U(1,numberx-2);
    for i=numberx-3:(-1):1
        Tnext(i,1)=(Y(i,1)-C(1,i)*Tnext(i+1,1))/U(1,i);
    end
    Tcurrent=Tnext';
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
Discreteerror=Texasolution-Tnumsolution;
Var(time)=var(Discreteerror);
end
figure
time=1:N;
plot(time,Var,'-r*')
xlabel('第i个deltt','fontsize',14)
ylabel('解析解与数值解的方差','fontsize',14)
title('方差与deltt取值的关系','fontsize',16)