function A=accuracy(Unit,Unumsolution)
%% Pre-processing
deltx=1/Unit;endx=1;
numberx=endx/deltx+1;
%calculate the accuracy of space DGp0+DGP0
I1=0;t=[-1/sqrt(5),0,1/sqrt(5)];W=[5/9,8/9,5/9];
k=1;%determine the correctness of the program
for x=0:deltx:endx-deltx
   for i=1:3
       xi=deltx/2*t(i)+0.5*(2*x+deltx);
       for m=1:numberx-1
           if xi>(m-1)*deltx&&xi<m*deltx
              fi=(sin(pi*xi)-Unumsolution(1,m))^2;k=k+1;
           end
       end
       I1=I1+W(i)*fi;
   end
end
I1=I1*0.5*deltx;
A=sqrt(I1);
end