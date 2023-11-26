function A2=Accuracyh(Unit)
%% Pre-processing
deltax=1/Unit;endx=1;
numberx=endx/deltax+1;
Grid=zeros(1,numberx);
for i=2:numberx-1
    %Grid(1,i)=(i-1)*deltax;
    Grid(1,i)=(i-1)*deltax+(0.1*rand(1)-0.05)*deltax;
end
Grid(1,numberx)=endx;
h=@(y)sin(pi*y);H=@(y)pi*cos(pi*y);
Unumsolution=zeros(2,Unit);
Ureconstruct=zeros(1,Unit);

%% Proceeding
%对h
%对h
for k=1:Unit
    Unumsolution(1,k)=(cos(pi*Grid(k))-cos(pi*Grid(k+1)))/(pi*(Grid(k+1)-Grid(k)));
    Unumsolution(2,k)=sin(pi*Grid(k+1))-sin(pi*Grid(k));%store Uxc*deltax
    Unumsolution(3,k)=-pi^2*sin(pi*(Grid(k+1)+Grid(k))/2)*(Grid(k+1)-Grid(k))^2;
end

%Reconstruct Uxx*deltaxi^2
for k=2:numberx-2
    xci=0.5*(Grid(k)+Grid(k+1));%xci
    xci1=0.5*(Grid(k+1)+Grid(k+2));%xci+1
    xci2=0.5*(Grid(k-1)+Grid(k));%xci-1
    deltaxi=Grid(k+1)-Grid(k);deltaxi1=Grid(k+2)-Grid(k+1);deltaxi2=Grid(k)-Grid(k-1);
    B4i1=((deltaxi/2+deltaxi1)^4-(deltaxi/2)^4)/(24*deltaxi^3)-(xci1-xci)*deltaxi1/(24*deltaxi);
    B4i2=((deltaxi/2+deltaxi2)^4-(deltaxi/2)^4)/(24*deltaxi^3)-(xci2-xci)*deltaxi2/(24*deltaxi);
    B3i1=(xci1-xci)^2*deltaxi1/(2*deltaxi^2)+(deltaxi1^3/deltaxi^2-deltaxi1)/24;
    B3i2=(xci2-xci)^2*deltaxi2/(2*deltaxi^2)+(deltaxi2^3/deltaxi^2-deltaxi2)/24;
    A=[
         B4i1/deltaxi1;
        B3i1/(deltaxi1*deltaxi);
%        deltaxi1*(3*(xci1-xci)^2-deltaxi^2/4)/deltaxi^3;
%         (xci1-xci)^2*deltaxi1^2/deltaxi^3;
         (xci1-xci)/deltaxi^3;%额外加的
         B4i2/deltaxi2;
        B3i2/(deltaxi2*deltaxi);
%          deltaxi2*(3*(xci2-xci)^2-deltaxi^2/4)/deltaxi^3;
%          (xci2-xci)^2*deltaxi2^2/deltaxi^3];
         (xci2-xci)/deltaxi^3];%额外加的
        
        
        
        
    b=[
         Unumsolution(1,k+1)-Unumsolution(1,k)-Unumsolution(2,k)*(xci1-xci)/deltaxi-Unumsolution(3,k)*B3i1/deltaxi1;
        Unumsolution(2,k+1)/deltaxi1-Unumsolution(2,k)/deltaxi-Unumsolution(3,k)*(xci1-xci)/deltaxi^2;
%          Unumsolution(2,k+1)-Unumsolution(2,k)*deltaxi1/deltaxi-Unumsolution(3,k)*deltaxi1*(xci1-xci)/deltaxi^2;
%          Unumsolution(3,k+1)-Unumsolution(3,k)*deltaxi1^2/deltaxi^2;
        Unumsolution(3,k+1)/deltaxi1^2-Unumsolution(3,k)/deltaxi^2;%额外加的
        Unumsolution(1,k-1)-Unumsolution(1,k)-Unumsolution(2,k)*(xci2-xci)/deltaxi-Unumsolution(3,k)*B3i2/deltaxi2;
        Unumsolution(2,k-1)/deltaxi2-Unumsolution(2,k)/deltaxi-Unumsolution(3,k)*(xci2-xci)/deltaxi^2;
%          Unumsolution(2,k-1)-Unumsolution(2,k)*deltaxi2/deltaxi-Unumsolution(3,k)*deltaxi2*(xci2-xci)/deltaxi^2;
%          Unumsolution(3,k-1)-Unumsolution(3,k)*deltaxi2^2/deltaxi^2];
        Unumsolution(3,k-1)/deltaxi2^2-Unumsolution(3,k)/deltaxi^2];%额外加的
                
    Ureconstruct(k)=A\b;
end
%boundary
%Left
k=1;
    xci=0.5*(Grid(k)+Grid(k+1));%xci
    xci1=0.5*(Grid(k+1)+Grid(k+2));%xci+1
    deltaxi=Grid(k+1)-Grid(k);deltaxi1=Grid(k+2)-Grid(k+1); 
    B4i1=((deltaxi/2+deltaxi1)^4-(deltaxi/2)^4)/(24*deltaxi^3)-(xci1-xci)*deltaxi1/(24*deltaxi);
    B3i1=(xci1-xci)^2*deltaxi1/(2*deltaxi^2)+(deltaxi1^3/deltaxi^2-deltaxi1)/24;

    
A=[
     B4i1/deltaxi1;
   B3i1/(deltaxi1*deltaxi);
%     deltaxi1*(xci1-xci)^2/(2*deltaxi^3);
%     (xci1-xci)^2*deltaxi1^2/deltaxi^3];
  (xci1-xci)/deltaxi^3];%额外加的

b=[
     Unumsolution(1,k+1)-Unumsolution(1,k)-Unumsolution(2,k)*(xci1-xci)/deltaxi-Unumsolution(3,k)*B3i1/deltaxi1;
        Unumsolution(2,k+1)/deltaxi1-Unumsolution(2,k)/deltaxi-Unumsolution(3,k)*(xci1-xci)/deltaxi^2;
%         Unumsolution(2,k+1)-Unumsolution(2,k)*deltaxi1/deltaxi-Unumsolution(3,k)*deltaxi1*(xci1-xci)/deltaxi^2;
%          Unumsolution(3,k+1)-Unumsolution(3,k)*deltaxi1^2/deltaxi^2]; 
        Unumsolution(3,k+1)/deltaxi1^2-Unumsolution(3,k)/deltaxi^2];%额外加的
    Ureconstruct(1)=A\b;
%Right
k=numberx-1;
    xci=0.5*(Grid(k)+Grid(k+1));%xci
    xci2=0.5*(Grid(k-1)+Grid(k));%xci-1
    deltaxi=Grid(k+1)-Grid(k);deltaxi2=Grid(k)-Grid(k-1);    
    B4i2=((deltaxi/2+deltaxi2)^4-(deltaxi/2)^4)/(24*deltaxi^3)-(xci2-xci)*deltaxi2/(24*deltaxi);
    B3i2=(xci2-xci)^2*deltaxi2/(2*deltaxi^2)+(deltaxi2^3/deltaxi^2-deltaxi2)/24;
    
    
 A=[
      B4i2/deltaxi2;
        B3i2/(deltaxi2*deltaxi);
%         deltaxi2*(xci2-xci)^2/(2*deltaxi^3);
%          (xci2-xci)^2*deltaxi2^2/deltaxi^3];
        (xci2-xci)/deltaxi^3];%额外加的
    
 b=[  
      Unumsolution(1,k-1)-Unumsolution(1,k)-Unumsolution(2,k)*(xci2-xci)/deltaxi-Unumsolution(3,k)*B3i2/deltaxi2;
        Unumsolution(2,k-1)/deltaxi2-Unumsolution(2,k)/deltaxi-Unumsolution(3,k)*(xci2-xci)/deltaxi^2;
%         Unumsolution(2,k-1)-Unumsolution(2,k)*deltaxi2/deltaxi-Unumsolution(3,k)*deltaxi2*(xci2-xci)/deltaxi^2;
%          Unumsolution(3,k-1)-Unumsolution(3,k)*deltaxi2^2/deltaxi^2]; 
              Unumsolution(3,k-1)/deltaxi2^2-Unumsolution(3,k)/deltaxi^2];%额外加的
 Ureconstruct(numberx-1)=A\b;
   
%calculate the accuracy of space
I1=0;t=[-1/sqrt(5),0,1/sqrt(5)];W=[5/9,8/9,5/9];
k=1;%determine the correctness of the program
for K=1:numberx-1
    xci=(Grid(K+1)+Grid(K))/2;
    deltaxi=Grid(K+1)-Grid(K);
   for i=1:3
       xi=(Grid(K+1)-Grid(K))/2*t(i)+0.5*(Grid(K+1)+Grid(K));
       %对v
%       fi=(pi*cos(pi*xi)-(Unumsolution(2,K)/deltaxi+Unumsolution(3,K)*(xi-xci)/deltaxi^2+Ureconstruct(K)*((xi-xci)^2/(2*deltaxi^2)-1/24)/deltaxi))^2;
      %对phi
      fi=(sin(pi*xi)-(Unumsolution(1,K)+Unumsolution(2,K)*(xi-xci)/deltaxi+Unumsolution(3,K)*((xi-xci)^2/(2*deltaxi^2)-1/24)+Ureconstruct(K)*(xi-xci)*((xi-xci)^2-deltaxi^2/4)/(6*deltaxi^3)))^2;

       I1=I1+W(i)*fi*0.5*(Grid(K+1)-Grid(K));        
   end
   k=k+1;
end
A2=sqrt(I1);
end
