clc
clear all
close all

%% parameter set up
format long
afa=0.06;
deltx=0.02;
T0=50;
endt=10;endx=1;n=(0.003-0.0001)/0.0001+1;
A=zeros(n,endx/deltx+1);
deltt=zeros(1,n);
for i=1:n
deltt(1,i)=0.0001*i;
end

%% calculate
for p=1:n
B=change(afa,deltx,deltt(1,p),T0,endt,endx);A(p,:)=B;
end
for i=1:n
    for j=1:endx/deltx+1
        if A(i,j)>0
            A(i,j)=1;
        elseif A(i,j)<0
                A(i,j)=-1;
        end
    end
end


% plot(deltt,A,'-*');
% xlabel('间隔时间deltt','fontsize',14)
% ylabel('解析解与数值解的误差方差','fontsize',14)
% title('T的数值解与解析解的误差方差','fontsize',16)