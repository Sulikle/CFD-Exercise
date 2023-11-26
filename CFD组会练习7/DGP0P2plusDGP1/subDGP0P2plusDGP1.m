function [n,Ul2errors,Vl2errors]=subDGP0P2plusDGP1(Unit,CFL,endtau,tol,belta,nsdv,nexplicit)
%% Preproceeding
%Some basic paramater
endx=1;deltax=endx/Unit;numberx=endx/deltax+1;
%��¼�ڵ�λ��,���¸���������belta
Grid=zeros(1,numberx);
Deltax=zeros(1,Unit);
for i=2:numberx-1
    Grid(1,i)=(i-1)*deltax+(2*rand(1)-1)*belta*deltax;
end
Grid(1,numberx)=endx;
for i=2:numberx
        Deltax(i-1)=Grid(1,i)-Grid(1,i-1);%��¼ÿ����Ԫ�����䳤��
end
Uexasolution=zeros(2,numberx);

 %solve the exasolution
for k=1:numberx
    Uexasolution(1,k)=sin(pi*Grid(k));
    Uexasolution(2,k)=pi*cos(pi*Grid(k));
end

%% ������ⷽ���ж�
%% Explicit
if nexplicit==1
        [Unumsolution,n]=P0P2plusP1Explicit(Unit,CFL,endtau,tol,Grid,Deltax,nsdv);
        Ul2errors=UL2errors(Unumsolution,Deltax,Grid);
        Vl2errors=VL2errors(Unumsolution,Deltax,Grid);
%% Explicit Euler    
    if nsdv==1

% Post-proceeding
%U
figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/Deltax(k)+Unumsolution(3,k)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on

 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/Deltax(k)+Unumsolution(3,k)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);
 end
 
%plot the exact

plot(Grid,Uexasolution(1,:),'-b*','linewidth',1.5)
H2=plot(Grid,Uexasolution(1,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DG(P0P2)+DG(P1)','������');
lgd.FontSize=12;
xlabel('λ��x','fontsize',14)
ylabel('��ֵU','fontsize',14)
 title('DG(P0P2)+DG(P1) Explicit Euler��ֵ���������(U)','fontsize',16)
hold off

% Ux
 figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(2,k)/Deltax(k)+Unumsolution(3,k)*(x-xci)/Deltax(k)^2;
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on
 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(2,k)/Deltax(k)+Unumsolution(3,k)*(x-xci)/Deltax(k)^2;
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);
 end

%exact
plot(Grid,Uexasolution(2,:),'-b*','linewidth',1.5)
H2=plot(Grid,Uexasolution(2,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DG(P0P2)+DG(P1)','������');
lgd.FontSize=12;
xlabel('λ��x','fontsize',14)
ylabel('��ֵU','fontsize',14)
title('DG(P0P2)+DG(P1) Explicit Euler��ֵ���������(Ux)','fontsize',16)
hold off

%% TVDRK3        
elseif nsdv==2
% Post-proceeding
%U
figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/Deltax(k)+Unumsolution(3,k)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on

 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/Deltax(k)+Unumsolution(3,k)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);
 end
 
%plot the exact

plot(Grid,Uexasolution(1,:),'-b*','linewidth',1.5)
H2=plot(Grid,Uexasolution(1,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DG(P0P2)+DG(P1)','������');
lgd.FontSize=12;
xlabel('λ��x','fontsize',14)
ylabel('��ֵU','fontsize',14)
 title('DG(P0P2)+DG(P1) TVDRK3 ��ֵ���������(U)','fontsize',16)
hold off       

% fprintf('L2���Ϊ %f',l2errors)
% Ux
 figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(2,k)/Deltax(k)+Unumsolution(3,k)*(x-xci)/Deltax(k)^2;
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on
 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(2,k)/Deltax(k)+Unumsolution(3,k)*(x-xci)/Deltax(k)^2;
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);
 end

%exact
plot(Grid,Uexasolution(2,:),'-b*','linewidth',1.5)
H2=plot(Grid,Uexasolution(2,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DG(P0P2)+DG(P1)','������');
lgd.FontSize=12;
xlabel('λ��x','fontsize',14)
ylabel('��ֵU','fontsize',14)
title('DG(P0P2)+DG(P1) TVDRK3 ��ֵ���������(Ux)','fontsize',16)
hold off



 end




%% Implicit
elseif nexplicit==0

[Unumsolution,n]=P0P2plusP1BDF1(Unit,CFL,endtau,tol,Grid,Deltax,nsdv);
Ul2errors=UL2errors(Unumsolution,Deltax,Grid);
Vl2errors=VL2errors(Unumsolution,Deltax,Grid);
%% Jacobi
if nsdv==1
% Post-proceeding
%U
figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/Deltax(k)+Unumsolution(3,k)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on

 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/Deltax(k)+Unumsolution(3,k)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);
 end
 
%plot the exact

plot(Grid,Uexasolution(1,:),'-b*','linewidth',1.5)
H2=plot(Grid,Uexasolution(1,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DG(P0P2)+DG(P1)','������');
lgd.FontSize=12;
xlabel('λ��x','fontsize',14)
ylabel('��ֵU','fontsize',14)
 title('DG(P0P2)+DG(P1)Jacobi��ֵ���������(U)','fontsize',16)
hold off
% Ux
 figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(2,k)/Deltax(k)+Unumsolution(3,k)*(x-xci)/Deltax(k)^2;
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on
 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(2,k)/Deltax(k)+Unumsolution(3,k)*(x-xci)/Deltax(k)^2;
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);
 end

%exact
plot(Grid,Uexasolution(2,:),'-b*','linewidth',1.5)
H2=plot(Grid,Uexasolution(2,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DG(P0P2)+DG(P1)','������');
lgd.FontSize=12;
xlabel('λ��x','fontsize',14)
ylabel('��ֵU','fontsize',14)
title('DG(P0P2)+DG(P1) Jacobi ��ֵ���������(Ux)','fontsize',16)
hold off

%% LUSGS
elseif nsdv==2

% Post-proceeding
%U
figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/Deltax(k)+Unumsolution(3,k)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on

 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(1,k)+Unumsolution(2,k)*(x-xci)/Deltax(k)+Unumsolution(3,k)*(0.5*((x-xci)/Deltax(k)).^2-1/24);
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);
 end
 
%plot the exact

plot(Grid,Uexasolution(1,:),'-b*','linewidth',1.5)
H2=plot(Grid,Uexasolution(1,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DG(P0P2)+DG(P1)','������');
lgd.FontSize=12;
xlabel('λ��x','fontsize',14)
ylabel('��ֵU','fontsize',14)
 title('DG(P0P2)+DG(P1)LUSGS��ֵ���������(U)','fontsize',16)
hold off


% Ux
 figure
 k=1;
 x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(2,k)/Deltax(k)+Unumsolution(3,k)*(x-xci)/Deltax(k)^2;
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);hold on
 H1=plot(x,y,'-r^','linewidth',1.5);hold on
 for k=2:numberx-1
     x=Grid(k):1*(Grid(k+1)-Grid(k)):Grid(k+1);
 xci=(Grid(k+1)+Grid(k))/2;
 p=@(x)Unumsolution(2,k)/Deltax(k)+Unumsolution(3,k)*(x-xci)/Deltax(k)^2;
 y=p(x);
 plot(x,y,'-r^','linewidth',1.5);
 end

%exact
plot(Grid,Uexasolution(2,:),'-b*','linewidth',1.5)
H2=plot(Grid,Uexasolution(2,:),'-b*','linewidth',1.5);hold on
lgd=legend([H1,H2],'DG(P0P2)+DG(P1)','������');
lgd.FontSize=12;
xlabel('λ��x','fontsize',14)
ylabel('��ֵU','fontsize',14)
title('DG(P0P2)+DG(P1) LUSGS ��ֵ���������(Ux)','fontsize',16)
hold off

end

end
fprintf('U_L2���Ϊ %f, V_L2���Ϊ %f',Ul2errors,Vl2errors)
end
