function A1=AccuracyU(Unit,Unumsolution,Grid)
%% Pre-processing
deltx=1/Unit;endx=1;
numberx=endx/deltx+1;
%calculate the accuracy of space
I1=0;t=[-1/sqrt(5),0,1/sqrt(5)];W=[5/9,8/9,5/9];
k=1;%determine the correctness of the program
for K=1:numberx-1
   for i=1:3
       xi=(Grid(K+1)-Grid(K))/2*t(i)+0.5*(Grid(K+1)+Grid(K));
       %for m=1:numberx-1
           %if xi>Grid(m)&&xi<Grid(m+1)
              fi=(sin(pi*xi)-(Unumsolution(1,k)+Unumsolution(2,k)*(xi-(Grid(k)+Grid(k+1))/2)))^2;
           %end
       %end
       I1=I1+W(i)*fi*0.5*(Grid(K+1)-Grid(K));
   end
   k=k+1;
end

A1=sqrt(I1);
end