function I=gauss2(a,b)
I=0;t=[-1/sqrt(5),0,1/sqrt(5)];W=[5/9,8/9,5/9];deltax=b-a;
xc=a+deltax/2;
   for i=1:3
       xi=(b-a)/2*t(i)+0.5*(a+b);
       fi=pi^2*sin(pi*xi)*(xi-xc)/deltax;
       I=I+W(i)*fi;
   end
I=I*0.5*(b-a);
end