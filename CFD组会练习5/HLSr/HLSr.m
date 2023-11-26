clc
clear all
close all
%% Pre-proceeding
Unit=64;endx=1;deltax=endx/Unit;numberx=Unit+1;
%记录内点位置,上下浮动不超过百分之5
Grid=zeros(1,numberx);
for i=2:numberx-1
    Grid(1,i)=(i-1)*deltax+(0.1*rand(1)-0.05)*deltax;
end
Grid(1,numberx)=endx;
f=@(x)x.^3+x.^2+x+1;F=@(x)3*x.^2+2*x+1;
h=@(x)sin(pi*x);H=@(x)pi*cos(pi*x);
Unumsolution=zeros(3,Unit);
Ureconstruct=zeros(1,Unit);
Unumsolution1=zeros(1,2);
Unumsolution2=zeros(2,numberx-1);
Acc=zeros(3,4);a1=[1/(2*64),1/(2*128),1/(2*256),1/(2*512)];a2=[1/(2*64),1/(2*128)];
%% Proceeding
%对f
for k=1:Unit
    %Un1存储的是 mean phi
    Unumsolution(1,k)=(Grid(k+1)-Grid(k)+0.5*(Grid(k+1)^2-Grid(k)^2)+(Grid(k+1)^3-Grid(k)^3)/3+(Grid(k+1)^4-Grid(k)^4)/4)/(Grid(k+1)-Grid(k));
    %Unumsolution(2,k)=F(0.5*(Grid(k)+Grid(k+1)))*(Grid(k+1)-Grid(k));%store Uxc*deltax
    %Un2存储的是 mean phix *deltaxi
    Unumsolution(2,k)=Grid(k+1)^3-Grid(k)^3+Grid(k+1)^2-Grid(k)^2+Grid(k+1)-Grid(k);
    %Un3存储的是phixxci * deltaxi^2
    Unumsolution(3,k)=(6*(Grid(k)+Grid(k+1))/2+2)*(Grid(k+1)-Grid(k))^2;
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
 
%Post-proceeding
figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 deltaxi=Grid(k+1)-Grid(k);
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/deltaxi+Unumsolution(3,k)*((x-xci).^2/(2*deltaxi^2)-1/24)+Ureconstruct(k)*(x-xci).*((x-xci).^2-deltaxi^2/4)/(6*deltaxi^3);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on

 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 deltaxi=Grid(k+1)-Grid(k);
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/deltaxi+Unumsolution(3,k)*((x-xci).^2/(2*deltaxi^2)-1/24)+Ureconstruct(k)*(x-xci).*((x-xci).^2-deltaxi^2/4)/(6*deltaxi^3);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);
 end
 
hold on
x=Grid(1):0.01*(Grid(numberx)-Grid(1)):Grid(numberx);
plot(x,f(x),'-b','linewidth',1.5);
H2=plot(x,f(x),'-b','linewidth',1.5);
lgd=legend([H1,H2],'Reconstruct','Exact solution');
lgd.FontSize=12;
xlabel('Position x','fontsize',14)
ylabel('Numerical value','fontsize',14)
title('Hyperbolic HLSr(NUG) f=x^3+x^2+x+1','fontsize',16)
hold off

%计算精度
Acc(1,1)=Accuracy(64);
Acc(1,2)=Accuracy(128);
Acc(1,3)=Accuracy(256);
Acc(1,4)=Accuracy(512);

for k=1:3
accuracyf(k)=(log10(Acc(1,k+1))-log10(Acc(1,k)))./(log10(a1(1,k+1))-log10(a1(1,k)));
end

% figure
% hold on
% plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5)
% H1=plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5);
% 
% H2=plot(log10(a2),1*log10(a2),'--','linewidth',1.5);
% plot(log10(a2),2*log10(a2),'--','linewidth',1.5)
% H3=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
% lgd=legend([H1,H2,H3],'v reconstruct','Slope=1','Slope=2');
% lgd.FontSize=12;
% xlabel('Log(1/DOF)','fontsize',14)
% ylabel('Log(episilo)','fontsize',14)
% title('v 精度分析 DG(P0P2)+rDG(P0P1)','fontsize',16)

figure
hold on
plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5)
H1=plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5);

H2=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
plot(log10(a2),3*log10(a2),'--','linewidth',1.5)
H3=plot(log10(a2),3*log10(a2),'--','linewidth',1.5);
lgd=legend([H1,H2,H3],'f reconstruct','Slope=2','Slope=3');
lgd.FontSize=12;
xlabel('Log(1/DOF)','fontsize',14)
ylabel('Log(episilo)','fontsize',14)
title('f(x)精度分析 DG(P0P3)+rDG(P1P2)','fontsize',16)



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
 


%Post-proceeding
figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 deltaxi=Grid(k+1)-Grid(k);
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/deltaxi+Unumsolution(3,k)*((x-xci).^2/(2*deltaxi^2)-1/24)+Ureconstruct(k)*(x-xci).*((x-xci).^2-deltaxi^2/4)/(6*deltaxi^3);
 y=p(x);
 plot(x,y,'-rh','linewidth',1.5);hold on
 H1=plot(x,y,'-rh','linewidth',1.5);hold on

 for k=2:numberx-1
    x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 deltaxi=Grid(k+1)-Grid(k);
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/deltaxi+Unumsolution(3,k)*((x-xci).^2/(2*deltaxi^2)-1/24)+Ureconstruct(k)*(x-xci).*((x-xci).^2-deltaxi^2/4)/(6*deltaxi^3);
 y=p(x);
 plot(x,y,'-rh','linewidth',1.5);
 end
 
hold on
x=Grid(1):0.01*(Grid(numberx)-Grid(1)):Grid(numberx);
plot(x,h(x),'-b','linewidth',1.5);
H2=plot(x,h(x),'-b','linewidth',1.5);
lgd=legend([H1,H2],'Reconstruct','Exact solution');
lgd.FontSize=12;
xlabel('Position x','fontsize',14)
ylabel('Numerical value','fontsize',14)
title('Hyperbolic HLSr(NUG) g(x)=sin(pi*x)','fontsize',16)
hold off

%计算精度
Acc(1,1)=Accuracyh(8);
Acc(1,2)=Accuracyh(16);
Acc(1,3)=Accuracyh(32);
Acc(1,4)=Accuracyh(64);
for k=1:3
accuracyh(k)=(log10(Acc(1,k+1))-log10(Acc(1,k)))./(log10(a1(1,k+1))-log10(a1(1,k)));
end

figure
hold on
plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5)
H1=plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5);

H2=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
plot(log10(a2),3*log10(a2),'--','linewidth',1.5)
H3=plot(log10(a2),3*log10(a2),'--','linewidth',1.5);
lgd=legend([H1,H2,H3],'h reconstruct','Slope=2','Slope=3');
lgd.FontSize=12;
xlabel('Log(1/DOF)','fontsize',14)
ylabel('Log(episilo)','fontsize',14)
title('g(x)精度分析 DG(P0P3)+rDG(P1P2)','fontsize',16)






