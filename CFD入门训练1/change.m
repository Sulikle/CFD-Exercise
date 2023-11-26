function B=change(afa,deltx,deltt,T0,endt,endx)
format long
numberx=endx/deltx+1;numbert=ceil(endt/deltt)+1;
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
k=randi([1,numbert]);
B1=A(k,:);
B2=zeros(1,numberx);
p=1;
for x=0:deltx:endx
    T=T0*sin(pi*x)*exp((-afa*(pi)^2)*(k-1)*deltt);B2(1,p)=T;p=p+1;
end
B2(1,p-1)=0;
%calculate the variance
B=B2-B1;
% Var=var(B);Vart=Var
end