function A2=Accuracyg(Unit)
%% Pre-processing
deltax=1/Unit;endx=1;
numberx=endx/deltax+1;
Grid=zeros(1,numberx);
for i=2:numberx-1
    %Grid(1,i)=(i-1)*deltax;
    Grid(1,i)=(i-1)*deltax+(0.1*rand(1)-0.05)*deltax;
end
Grid(1,numberx)=endx;
g=@(y)sin(pi*y);G=@(y)pi*cos(pi*y);
Unumsolution=zeros(2,Unit);
Ureconstruct=zeros(1,Unit);

%% Proceeding
%¶Ôg
for k=1:Unit
    Unumsolution(1,k)=(cos(pi*Grid(k))-cos(pi*Grid(k+1)))/(pi*(Grid(k+1)-Grid(k)));
    Unumsolution(2,k)=G(0.5*(Grid(k)+Grid(k+1)))*(Grid(k+1)-Grid(k));%store Uxc*deltax
end

%Reconstruct Uxx
for k=2:numberx-2
    xci=0.5*(Grid(k)+Grid(k+1));%xci
    xci1=0.5*(Grid(k+1)+Grid(k+2));%xci+1
    xci2=0.5*(Grid(k-1)+Grid(k));%xci-1
    A=[(xci1-xci)^2/(2*(Grid(k+1)-Grid(k))^2)+1/24*((Grid(k+2)-Grid(k+1))^2/(Grid(k+1)-Grid(k))^2-1);
        (xci1-xci)/(Grid(k+1)-Grid(k))^2;
        (xci2-xci)^2/(2*(Grid(k+1)-Grid(k))^2)+1/24*((Grid(k)-Grid(k-1))^2/(Grid(k+1)-Grid(k))^2-1);
        (xci2-xci)/(Grid(k+1)-Grid(k))^2];
    b=[Unumsolution(1,k+1)-Unumsolution(1,k)-Unumsolution(2,k)*(xci1-xci)/(Grid(k+1)-Grid(k));
        Unumsolution(2,k+1)/(Grid(k+2)-Grid(k+1))-Unumsolution(2,k)/(Grid(k+1)-Grid(k));
        Unumsolution(1,k-1)-Unumsolution(1,k)-Unumsolution(2,k)*(xci2-xci)/(Grid(k+1)-Grid(k));
        Unumsolution(2,k-1)/(Grid(k)-Grid(k-1))-Unumsolution(2,k)/(Grid(k+1)-Grid(k))];
    Ureconstruct(k)=A\b;
end
%boundary
%Left
    xci=0.5*(Grid(1)+Grid(2));%xci
    xci1=0.5*(Grid(2)+Grid(3));%xci+1
A=[(xci1-xci)^2/(2*(Grid(2)-Grid(1))^2)+1/24*((Grid(3)-Grid(2))^2/(Grid(2)-Grid(1))^2-1);
        (xci1-xci)/(Grid(1+1)-Grid(1))^2];
b=[Unumsolution(1,2)-Unumsolution(1,1)-Unumsolution(2,1)*(xci1-xci)/(Grid(2)-Grid(1));
        Unumsolution(2,1+1)/(Grid(1+2)-Grid(1+1))-Unumsolution(2,1)/(Grid(1+1)-Grid(1))];    
    Ureconstruct(1)=A\b;
%Right    
    xci=0.5*(Grid(numberx-1)+Grid(numberx-1+1));%xci
    xci2=0.5*(Grid(numberx-1-1)+Grid(numberx-1));%xci-1    
 A=[(xci2-xci)^2/(2*(Grid(numberx-1+1)-Grid(numberx-1))^2)+1/24*((Grid(numberx-1)-Grid(numberx-1-1))^2/(Grid(numberx-1+1)-Grid(numberx-1))^2-1);
        (xci2-xci)/(Grid(numberx-1+1)-Grid(numberx-1))^2];
 b=[Unumsolution(1,numberx-1-1)-Unumsolution(1,numberx-1)-Unumsolution(2,numberx-1)*(xci2-xci)/(Grid(numberx-1+1)-Grid(numberx-1));
        Unumsolution(2,numberx-1-1)/(Grid(numberx-1)-Grid(numberx-1-1))-Unumsolution(2,numberx-1)/(Grid(numberx-1+1)-Grid(numberx-1))];   
       Ureconstruct(numberx-1)=A\b;      
%calculate the accuracy of space
I1=0;t=[-1/sqrt(5),0,1/sqrt(5)];W=[5/9,8/9,5/9];
k=1;%determine the correctness of the program
for K=1:numberx-1
   for i=1:3
       xi=(Grid(K+1)-Grid(K))/2*t(i)+0.5*(Grid(K+1)+Grid(K));

              fi=(sin(pi*xi)-(Unumsolution(1,K)+Unumsolution(2,K)*(xi-(Grid(K)+Grid(K+1))/2)/(Grid(K+1)-Grid(K))+Ureconstruct(1,K)*(0.5*((xi-(Grid(K)+Grid(K+1))/2)/(Grid(K+1)-Grid(K)))^2-1/24)))^2;

       I1=I1+W(i)*fi*0.5*(Grid(K+1)-Grid(K));        
   end
   k=k+1;
end
A2=sqrt(I1);
end
