function A2=Accuracy(Unit)
%% Pre-processing
deltax=1/Unit;endx=1;
numberx=endx/deltax+1;
Grid=zeros(1,numberx);
for i=2:numberx-1
    Grid(1,i)=(i-1)*deltax+(0.1*rand(1)-0.05)*deltax;
end
Grid(1,numberx)=endx;
f=@(x)x.^2+x+1;F=@(x)2*x+1;
g=@(y)sin(pi*y);G=@(y)pi*cos(pi*y);
Unumsolution=zeros(2,Unit);
Ureconstruct=zeros(1,Unit);
Un=zeros(1,numberx);

%% Proceeding
%¶Ôf
for k=1:Unit
    Unumsolution(1,k)=f(0.5*(Grid(k)+Grid(k+1)));
    Unumsolution(2,k)=F(0.5*(Grid(k)+Grid(k+1)));
end

%Reconstruct Uxx
for iface=2:numberx-1
    ieL=iface-1;
    ieR=iface;
    Un(iface)=((Grid(ieL+1)-Grid(ieL))*Unumsolution(2,ieL)+(Grid(ieR+1)-Grid(ieR))*Unumsolution(2,ieR))/((Grid(ieL+1)-Grid(ieL))+(Grid(ieR+1)-Grid(ieR)));
    Ureconstruct(ieL)=Ureconstruct(ieL)+Un(iface)/deltax;
    Ureconstruct(ieR)=Ureconstruct(ieR)-Un(iface)/deltax;
end
    Un(1)=Unumsolution(2,1);
    Ureconstruct(1)=Ureconstruct(1)-Un(1)/deltax;
    Un(numberx)=0.5*(Unumsolution(2,numberx-1)+Unumsolution(2,numberx-1));
    Ureconstruct(numberx-1)=Ureconstruct(numberx-1)+Un(numberx)/deltax;
%calculate the accuracy of space
I1=0;t=[-1/sqrt(5),0,1/sqrt(5)];W=[5/9,8/9,5/9];
k=1;%determine the correctness of the program
for K=1:numberx-1
   for i=1:3
       xi=(Grid(K+1)-Grid(K))/2*t(i)+0.5*(Grid(K+1)+Grid(K));
 %      for m=1:numberx-1
          % if xi>Grid(m)&&xi<Grid(m+1)
              fi=(xi^2+xi+1-(Unumsolution(1,K)+Unumsolution(2,K)*(xi-(Grid(K)+Grid(K+1))/2)+0.5*Ureconstruct(1,K)*(xi-(Grid(K)+Grid(K+1))/2)^2))^2;
           %end
  %     end
       I1=I1+W(i)*fi*0.5*(Grid(K+1)-Grid(K));        
   end
   k=k+1;
end
A2=sqrt(I1);
end
