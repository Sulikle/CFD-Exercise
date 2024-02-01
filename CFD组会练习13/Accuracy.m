function [Ul2error,Uxl2error,Uxl2error0]=Accuracy(Unit)
%% Preproceeding
endx=1;deltax=endx/Unit;npoint=Unit+1;
A=zeros(npoint,npoint);b=zeros(npoint,1);%matrix default setting
unumsolution=zeros(npoint,1);
uxnumsolution=zeros(npoint,1);
uxnumsolution0=zeros(Unit,1);
% Ue=@(x) 5*x.*(x-2);
% Uex=@(x)10*x-10;



%% Proceeding
%set A
for element=1:Unit
    ip1=element;
    ip2=element+1;
    A(ip1,ip1)=A(ip1,ip1)+1/deltax;
    A(ip1,ip2)=A(ip1,ip2)-1/deltax;
    A(ip2,ip1)=A(ip2,ip1)-1/deltax;
    A(ip2,ip2)=A(ip2,ip2)+1/deltax;
end
%Due to initial setting
for ip=2:npoint
    A(1,ip)=0;
    A(ip,1)=0;
end

%set b
b(1,1)=0;%initial setting
for ip=2:npoint-1
    b(ip,1)=-10*deltax;
end
b(npoint,1)=-5*deltax;

%solve the equation
[L,U]=lu(A);
y=L\b;
unumsolution=U\y;

%calculate Ux
for ip=2:npoint-1
    uxnumsolution(ip,1)=0.5*(unumsolution(ip+1,1)/deltax-unumsolution(ip-1,1)/deltax);
end
uxnumsolution(1,1)=unumsolution(2,1)/deltax-unumsolution(1,1)/deltax;

uxnumsolution(npoint,1)=0;
%Ux是常值
for element=1:Unit
    uxnumsolution0(element,1)=unumsolution(element+1,1)/deltax-unumsolution(element,1)/deltax;
end


%calculate the accuracy of spaceU
I1=0;
t=[-sqrt(15)/5,0,sqrt(15)/5];
% t=[-1/sqrt(5),0,1/sqrt(5)];
W=[5/9,8/9,5/9];
k=1;%determine the correctness of the program
for K=1:Unit
   for i=1:3
       xi=deltax/2*t(i)+0.5*(2*K-1)*deltax;
      %对phi
       fi=(5*xi*(xi-2)-unumsolution(K,1)-(unumsolution(K+1,1)-unumsolution(K,1))/deltax*(xi-(K-1)*deltax))^2;
       I1=I1+W(i)*fi*0.5*(deltax);        
   end
   k=k+1;
end
Ul2error=sqrt(I1);

%calculate the accuracy of spaceUx
I1=0;
t=[-sqrt(15)/5,0,sqrt(15)/5];
% t=[-1/sqrt(5),0,1/sqrt(5)];
W=[5/9,8/9,5/9];
k=1;%determine the correctness of the program
for K=1:Unit
   for i=1:3
       xi=deltax/2*t(i)+0.5*(2*K-1)*deltax;
      %对phi
       fi=(10*xi-10-uxnumsolution(K,1)-(uxnumsolution(K+1,1)-uxnumsolution(K,1))/deltax*(xi-(K-1)*deltax))^2;
       I1=I1+W(i)*fi*0.5*(deltax);        
   end
   k=k+1;
end
Uxl2error=sqrt(I1);

%calculate the accuracy of spaceUx0
I1=0;
t=[-sqrt(15)/5,0,sqrt(15)/5];
% t=[-1/sqrt(5),0,1/sqrt(5)];
W=[5/9,8/9,5/9];
k=1;%determine the correctness of the program
for K=1:Unit
   for i=1:3
       xi=deltax/2*t(i)+0.5*(2*K-1)*deltax;
      %对phi
       fi=(10*xi-10-uxnumsolution0(K,1))^2;
       I1=I1+W(i)*fi*0.5*(deltax);        
   end
   k=k+1;
end
Uxl2error0=sqrt(I1);

end
