clc 
clear
close all
%% Preproceeding
Acc=zeros(1,4);a1=zeros(1,4);a2=zeros(1,2);
for number=1:4
[Acc(1,number),a1(1,number)]=Accuracy(number);
end

for number=1:4
a1(1,number)=1/(3*a1(1,number));
end

a2=a1(1,3:4);
%º∆À„order
Accuracyvarphi=zeros(1,3);
for k=1:3
Accuracyvarphi(k)=(log10(Acc(1,k+1))-log10(Acc(1,k)))./(log10(a1(1,k+1))-log10(a1(1,k)));
end

figure
hold on
plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5)
H1=plot(log10(a1),log10(Acc(1,:)),'-c*','linewidth',1.5);

H2=plot(log10(a2),1*log10(a2),'--','linewidth',1.5);
 plot(log10(a2),2*log10(a2),'--','linewidth',1.5)
H3=plot(log10(a2),2*log10(a2),'--','linewidth',1.5);
lgd=legend([H1,H2,H3],'FEM','Slope=1','Slope=2');
lgd.FontSize=12;
xlabel('Log(1/DOF)','fontsize',14)
ylabel('Log(episilo)','fontsize',14)
title('2D FEM Velocity potential Accuracy analysis','fontsize',16)
grid on
hold off