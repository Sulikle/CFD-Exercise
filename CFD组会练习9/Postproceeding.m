clc
clear all
close all
%% Nelem=8
Re=[0,1,2,10,20];
%DGP0P1plusDGP0
U_DGP0P1plusDGP0=[0.0655,0.054,0.049,0.0098,0.0081];
Ux_DGP0P1plusDGP0=[0.2540,0.2625,0.2684,0.2781,0.2616];
%DGP0P2plusDGP1
U_DGP0P2plusDGP1=[0.0123,0.0188,0.0178,0.0045,0.0025];
Ux_DGP0P2plusDGP1=[0.0399,0.0602,0.0589,0.0738,0.0207];
%DGP0P2plusrDGP0P1
U_DGP0P2plusrDGP0P1=[0.0423,0.0414,0.0388,0.0185,0.0086];
Ux_DGP0P2plusrDGP0P1=[0.0678,0.0747,0.0764,0.0785,0.0732];

figure;
plot(Re,U_DGP0P1plusDGP0,'-r^','linewidth',1.5);
H1=plot(Re,U_DGP0P1plusDGP0,'-r^','linewidth',1.5);
hold on
plot(Re,U_DGP0P2plusDGP1,'-g^','linewidth',1.5);
H2=plot(Re,U_DGP0P2plusDGP1,'-g^','linewidth',1.5);
plot(Re,U_DGP0P2plusrDGP0P1,'-b^','linewidth',1.5);
H3=plot(Re,U_DGP0P2plusrDGP0P1,'-b^','linewidth',1.5);
lgd=legend([H1,H2,H3],'DG(P0P1)+DG(P0)','DG(P0P2)+DG(P1)','DG(P0P2)+rDG(P0P1)');
lgd.FontSize=12;
xlabel('雷诺数Re','fontsize',14)
ylabel('L2-errors','fontsize',14)
 title('L2-errors与Re的关系(U,Nelem=8)','fontsize',16)
hold off

figure;
plot(Re,Ux_DGP0P1plusDGP0,'-r^','linewidth',1.5);
H1=plot(Re,Ux_DGP0P1plusDGP0,'-r^','linewidth',1.5);
hold on
plot(Re,Ux_DGP0P2plusDGP1,'-g^','linewidth',1.5);
H2=plot(Re,Ux_DGP0P2plusDGP1,'-g^','linewidth',1.5);
plot(Re,Ux_DGP0P2plusrDGP0P1,'-b^','linewidth',1.5);
H3=plot(Re,Ux_DGP0P2plusrDGP0P1,'-b^','linewidth',1.5);
lgd=legend([H1,H2,H3],'DG(P0P1)+DG(P0)','DG(P0P2)+DG(P1)','DG(P0P2)+rDG(P0P1)');
lgd.FontSize=12;
xlabel('雷诺数Re','fontsize',14)
ylabel('L2-errors','fontsize',14)
 title('L2-errors与Re的关系(Ux,Nelem=8)','fontsize',16)
hold off

%% Nelem=64
Re=[0,1,2,10,20];
%DGP0P1plusDGP0
U_DGP0P1plusDGP0=[0.0087,0.0084,0.0080,0.0023,5.2662*10^(-4)];
Ux_DGP0P1plusDGP0=[0.0315,0.0317,0.0323,0.0340,0.0324];
%DGP0P2plusDGP1
U_DGP0P2plusDGP1=[1.45*10^(-4),1.43*10^(-4),1.4299*10^(-4),6.74*10^(-5),3.67*10^(-5)];
Ux_DGP0P2plusDGP1=[4.88*10^(-4),4.85*10^(-4),4.95*10^(-4),3.52*10^(-4),2.84*10^(-4)];
%DGP0P2plusrDGP0P1
U_DGP0P2plusrDGP0P1=[6.4574*10^(-4),6.1610*10^(-4),5.67*10^(-4),3.21*10^(-4),1.75*10^(-4)];
Ux_DGP0P2plusrDGP0P1=[6.0662*10^(-4),6.1027*10^(-4),6.12*10^(-4),6.27*10^(-4),6.01*10^(-4)];

figure;
plot(Re,U_DGP0P1plusDGP0,'-r^','linewidth',1.5);
H1=plot(Re,U_DGP0P1plusDGP0,'-r^','linewidth',1.5);
hold on
plot(Re,U_DGP0P2plusDGP1,'-g^','linewidth',1.5);
H2=plot(Re,U_DGP0P2plusDGP1,'-g^','linewidth',1.5);
plot(Re,U_DGP0P2plusrDGP0P1,'-b^','linewidth',1.5);
H3=plot(Re,U_DGP0P2plusrDGP0P1,'-b^','linewidth',1.5);
lgd=legend([H1,H2,H3],'DG(P0P1)+DG(P0)','DG(P0P2)+DG(P1)','DG(P0P2)+rDG(P0P1)');
lgd.FontSize=12;
xlabel('雷诺数Re','fontsize',14)
ylabel('L2-errors','fontsize',14)
 title('L2-errors与Re的关系(U,Nelem=64)','fontsize',16)
hold off

figure;
plot(Re,Ux_DGP0P1plusDGP0,'-r^','linewidth',1.5);
H1=plot(Re,Ux_DGP0P1plusDGP0,'-r^','linewidth',1.5);
hold on
plot(Re,Ux_DGP0P2plusDGP1,'-g^','linewidth',1.5);
H2=plot(Re,Ux_DGP0P2plusDGP1,'-g^','linewidth',1.5);
plot(Re,Ux_DGP0P2plusrDGP0P1,'-b^','linewidth',1.5);
H3=plot(Re,Ux_DGP0P2plusrDGP0P1,'-b^','linewidth',1.5);
lgd=legend([H1,H2,H3],'DG(P0P1)+DG(P0)','DG(P0P2)+DG(P1)','DG(P0P2)+rDG(P0P1)');
lgd.FontSize=12;
xlabel('雷诺数Re','fontsize',14)
ylabel('L2-errors','fontsize',14)
 title('L2-errors与Re的关系(Ux,Nelem=64)','fontsize',16)
hold off


