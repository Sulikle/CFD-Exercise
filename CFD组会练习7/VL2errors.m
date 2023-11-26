function l2errors=VL2errors(Unumsolution,Deltax,Grid)
numberx=size(Unumsolution,2)+1;

if size(Unumsolution,1)==2
    %calculate the accuracy of space
I1=0;
t=[-sqrt(15)/5,0,sqrt(15)/5];
W=[5/9,8/9,5/9];
%determine the correctness of the program
for K=1:numberx-1
    xci=0.5*(Grid(K+1)+Grid(K));
   for i=1:3
       xi=Deltax(K)/2*t(i)+xci;
       fi=(pi*cos(pi*xi)-(Unumsolution(2,K)/Deltax(K)))^2;
       I1=I1+W(i)*fi*0.5*Deltax(K);        
   end
end
l2errors=sqrt(I1);

elseif size(Unumsolution,1)==3    
%calculate the accuracy of space
I1=0;
t=[-sqrt(15)/5,0,sqrt(15)/5];
W=[5/9,8/9,5/9];
%determine the correctness of the program
for K=1:numberx-1
    xci=0.5*(Grid(K+1)+Grid(K));
   for i=1:3
       xi=Deltax(K)/2*t(i)+xci;
       fi=(pi*cos(pi*xi)-(Unumsolution(2,K)/Deltax(K)+Unumsolution(3,K)*(xi-xci)/Deltax(K)^2))^2;
       I1=I1+W(i)*fi*0.5*Deltax(K);        
   end
end
l2errors=sqrt(I1);

elseif size(Unumsolution,1)==4 
I1=0;
t=[0.8611363,0.3399810,-0.8611363,-0.3399810];
W=[0.3478548,0.6521452,0.3478548,0.6521452];
%determine the correctness of the program
for K=1:numberx-1
    xci=0.5*(Grid(K+1)+Grid(K));
   for i=1:4
       xi=Deltax(K)/2*t(i)+xci;
       fi=(pi*cos(pi*xi)-(Unumsolution(2,K)/Deltax(K)+Unumsolution(3,K)*(xi-xci)/Deltax(K)^2+Unumsolution(4,K)*(0.5*((xi-xci)/Deltax(K))^2-1/24)/Deltax(K)))^2;
       I1=I1+W(i)*fi*0.5*Deltax(K);        
   end
end
l2errors=sqrt(I1);        
end
      
end