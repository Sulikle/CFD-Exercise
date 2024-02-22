function [l2errors,Nelement]=Accuracy(number)
%% Preproceeding
%Input data
addpath(genpath('/Users/mac/Documents/CFD-Exercise-main/CFDProject1/Proj1_Mesh'))
% store the nodal points corresponding to each element in Edata————INPOEL
if number==1
feflo_bump = fopen('Proj1_Mesh/feflo.domn.cylinder.coarse', 'rb');%open the file
elseif number==2
    feflo_bump = fopen('Proj1_Mesh/feflo.domn.cylinder.medium', 'rb');%open the file
elseif number==3
    feflo_bump = fopen('Proj1_Mesh/feflo.domn.cylinder.fine', 'rb');%open the file
elseif number==4
    feflo_bump = fopen('Proj1_Mesh/feflo.domn.cylinder.vfine', 'rb');%open the file
end
skiprows=fscanf(feflo_bump, '%lf',1);% skip some rows
for i = 1:skiprows+2
    line=fgetl(feflo_bump);
end
Ndata = fscanf(feflo_bump, '%lf',[2,1]);
%read dimention&Nnode
Ndim=Ndata(1,1);Nnode=Ndata(2,1);
for i = 1:2
    line=fgetl(feflo_bump);
end
Ndata = fscanf(feflo_bump, '%lf',[3,1]);
Nelement=Ndata(1,1);Npoint=Ndata(2,1);Nface=Ndata(3,1);%Information of elements&points&boundary

for i = 1:3
    line=fgetl(feflo_bump);
end
%read data per row
C = textscan(line, '%f');
numbers = length(str2num(line));
INPOEL = fscanf(feflo_bump, '%lf',[numbers,Nelement]);
INPOEL=[C{1,1},INPOEL];
INPOEL=INPOEL(2:4,:);%每一列存储的是这个单元的第1，2，3个点

%store the coordinates of the points————COORD
for i = 1:1% skip some rows
    line=fgetl(feflo_bump);
end
COORD = fscanf(feflo_bump, '%lf',[Ndim+1,Npoint]);
COORD=COORD(2:3,:);%每一列存储的是这个点的横纵坐标

%store Boundary conditions
for i = 1:2+Npoint+2
    line=fgetl(feflo_bump);
end
%由于部分数据缺失，因此BOCND采取这种读取方式
BCOND=zeros(3,Nface);
for i=1:Nface
    C = textscan(line, '%f');
    BCOND(:,i)=C{1,1}(2:4,1);
    line=fgetl(feflo_bump);
end

% C = textscan(line, '%f');
% numbers = length(str2num(line));
% BCOND=fscanf(feflo_bump, '%lf',[numbers,Nface]);
% BCOND=[C{1,1},BCOND];
% BCOND=BCOND(2:4,:);%每一列存储的是该边界处的左右点以及flag，2代表固定面，4代表无穷远处
fclose(feflo_bump);

%elements surrounding points
ESUP2=zeros(1,Npoint+1);
for ie=1:Nelement
    for in=1:Nnode
        ESUP2(INPOEL(in,ie)+1)=ESUP2(INPOEL(in,ie)+1)+1;%斜对面存储该point涉及到的单元数
    end
end

for IPOIN=2:Npoint+1
    ESUP2(IPOIN)=ESUP2(IPOIN)+ESUP2(IPOIN-1);%将这些单元标号
end
MESUP=ESUP2(Npoint+1);
ESUP1=zeros(1,MESUP);
for ie=1:Nelement
    for in=1:Nnode
        IPOIN=INPOEL(in,ie);%点的标号
        ISTOR=ESUP2(IPOIN)+1;%改点所用单元的初地址
        ESUP2(IPOIN)=ISTOR;%更新地址
        ESUP1(ISTOR)=ie;%记录单元
    end
end
for IPOIN=Npoint+1:-1:2
    ESUP2(IPOIN)=ESUP2(IPOIN-1);
end
ESUP2(1)=0;

%stiffness matrix(LHS)&load vector(RHS)
LHS=zeros(Npoint,Npoint);RHS=zeros(Npoint,1);
HP=zeros(Nnode,Nnode,Nelement);%help matrix
Bv=[1,0];%boundary speed
varphi_exact=zeros(Npoint,1);

% v_element=zeros(Ndim,Nelement);%v in element
% v_point=zeros(Ndim,Npoint);%v in point
% v_point_scalar=zeros(1,Npoint);%store the size of v
% v_point_scalar_face=zeros(1,Npoint);%store the size of v in boundary



%% Proceeding
%construct help matrix
HPcycle=zeros(1,Nnode+2);c_ie=zeros(1,3);%c per element
for ie=1:Nelement
    ipoin1=INPOEL(1,ie);ipoin2=INPOEL(2,ie);ipoin3=INPOEL(3,ie);
    HPcycle=[ipoin1,ipoin2,ipoin3,ipoin1,ipoin2];
    for k=1:Nnode%先计算c
           c_ie(k)=COORD(1,HPcycle(k+1))*COORD(2,HPcycle(k+2))-COORD(1,HPcycle(k+2))*COORD(2,HPcycle(k+1));
    end
    for i=1:Nnode
        a_i=COORD(2,HPcycle(i+1))-COORD(2,HPcycle(i+2));b_i=-(COORD(1,HPcycle(i+1))-COORD(1,HPcycle(i+2)));
        for j=1:Nnode
        a_j=COORD(2,HPcycle(j+1))-COORD(2,HPcycle(j+2));b_j=-(COORD(1,HPcycle(j+1))-COORD(1,HPcycle(j+2)));
        HP(i,j,ie)=(a_i*a_j+b_i*b_j)/(2*sum(c_ie));
        end
    end
end

%construct stiffness matrix(LHS)
for ie=1:Nelement
    for i=1:Nnode
        ipoin=INPOEL(i,ie);
        for j=1:Nnode
            jpoin=INPOEL(j,ie);
            LHS(ipoin,jpoin)=LHS(ipoin,jpoin)+HP(i,j,ie);
        end
    end
end


%construct load vector
for iface=1:Nface
    ip1=BCOND(1,iface);ip2=BCOND(2,iface);flag=BCOND(3,iface);ip3=-1;%信息提取
    gama_e=sqrt((COORD(1,ip1)-COORD(1,ip2))^2+(COORD(2,ip1)-COORD(2,ip2))^2);%length of face
    NESP=ESUP2(1,ip1+1)-ESUP2(1,ip1);%number of elements surrounding ip1
    je=-1;%element of this boundary
    for ie=ESUP2(1,ip1)+1:ESUP2(1,ip1+1)
        for in=1:Nnode
            if INPOEL(in,ESUP1(ie))==ip2
                je=ESUP1(ie);
                break;
            end
        end
        if je~=-1
            break;
        end
    end
    for in=1:Nnode
        if INPOEL(in,je)~=ip1&&INPOEL(in,je)~=ip2
            ip3=INPOEL(in,je);%remark the third point of this element
        end
    end

    Bvector1=[COORD(1,ip2)-COORD(1,ip1),COORD(2,ip2)-COORD(2,ip1)];
    if Bv(1)==0&&Bv(2)==0
        fprintf('Please input a vector which is not zero' );
    end
    
    dot_product=dot(Bv,Bvector1);abs_dot_product=abs(dot_product);%向量点积
    length_v1=norm(Bv);length_v2=norm(Bvector1);length_product=length_v1*length_v2;%向量长度


    Bvector2=[COORD(1,ip3)-COORD(1,ip1),COORD(2,ip3)-COORD(2,ip1)];%该单元内部点与边界点构成的向量
    V_a=COORD(1,ip2)-COORD(1,ip1);V_b=COORD(2,ip2)-COORD(2,ip1);
    if V_b==0
        V_product1=[0,1];V_product2=[0,-1];
    else
    V_product1=[-V_b,V_a]/norm([-V_b,V_a]);%法线向量
    V_product2=[V_b,-V_a]/norm([V_b,-V_a]);
    end

    if abs_dot_product-length_product==0%判断速度向量与切线是否平行
        g=0;
    elseif abs_dot_product-length_product~=0&&dot(Bvector2,V_product1)<0&&flag==4%判断外法线向量并判断是否是无穷远处
        g=dot(Bv,V_product1);
    elseif abs_dot_product-length_product~=0&&dot(Bvector2,V_product2)<0&&flag==4%判断外法线向量并判断是否是无穷远处
         g=dot(Bv,V_product2);
    elseif abs_dot_product-length_product~=0&&flag==2%判断是否是固定面
        g=0;
%     elseif abs_dot_product-length_product~=0&&dot(Bvector2,V_product2)<0&&flag==2
%         g=0;
    end
    RHS(ip1)=RHS(ip1)+0.5*g*gama_e;
    RHS(ip2)=RHS(ip2)+0.5*g*gama_e;
end
%固定边界面
a=0.5;V_x=Bv(1,1);C_big=1e10;
for iface=1:Nface
    ip1=BCOND(1,iface);ip2=BCOND(2,iface);
    varphi_exact_1=V_x*COORD(1,ip1)*(1+a^2/(COORD(1,ip1)^2+COORD(2,ip1)^2));
    varphi_exact_2=V_x*COORD(1,ip2)*(1+a^2/(COORD(1,ip2)^2+COORD(2,ip2)^2));
    LHS(ip1,ip1)=C_big*LHS(ip1,ip1);LHS(ip2,ip2)=C_big*LHS(ip2,ip2);
    RHS(ip1,1)=LHS(ip1,ip1)*varphi_exact_1;RHS(ip2,1)=LHS(ip2,ip2)*varphi_exact_2;
end

% [l,u]=lu(LHS);
% y=l\RHS;
% varphi=u\y;
%solve the equation 
endtimes=10000;tol=10^(-16);
A=LHS;b=RHS;b1=b;varphi=zeros(Npoint,1);%initial value
varphi_start=0;
for times=2:endtimes
    %calculate X Forward sweep
    for i=1:Npoint
        for j=1:i-1
            b(i)=b(i)-A(i,j)*varphi(j);
        end
        for j=i+1:Npoint
            b(i)=b(i)-A(i,j)*varphi(j);
        end
        varphi(i)=b(i)/A(i,i);
    end
     b=b1;
     %Backward sweep
    for i=Npoint:-1:1
        for j=Npoint:-1:i+1
            b(i)=b(i)-A(i,j)*varphi(j);
        end
        for j=i-1:-1:1
            b(i)=b(i)-A(i,j)*varphi(j);
        end
        varphi(i)=b(i)/A(i,i);
    end
    b=b1;
    ratio=max(abs(varphi-varphi_start));
    varphi_start=varphi;
    if ratio<tol
        break;
    end
end
t=[-sqrt(15)/5,0,sqrt(15)/5];
W=[5/9,8/9,5/9];
l2errors=0;
 for ie=1:Nelement
     ipoin1=INPOEL(1,ie);ipoin2=INPOEL(2,ie);ipoin3=INPOEL(3,ie);c_ie=zeros(1,Nnode);a_ie=zeros(1,Nnode);b_ie=zeros(1,Nnode);
     x_point=[COORD(1,ipoin1),COORD(1,ipoin2),COORD(1,ipoin3)];y_point=[COORD(2,ipoin1),COORD(2,ipoin2),COORD(2,ipoin3)];
    HPcycle=[ipoin1,ipoin2,ipoin3,ipoin1,ipoin2];
    for k=1:Nnode%先计算c
           c_ie(k)=COORD(1,HPcycle(k+1))*COORD(2,HPcycle(k+2))-COORD(1,HPcycle(k+2))*COORD(2,HPcycle(k+1));
           a_ie(k)=COORD(2,HPcycle(k+1))-COORD(2,HPcycle(k+2));b_ie(k)=-(COORD(1,HPcycle(k+1))-COORD(1,HPcycle(k+2)));
    end

    D=sum(c_ie);
    Jacobi=(x_point(1,2)-x_point(1,1))*(y_point(1,3)-y_point(1,1))-(x_point(1,3)-x_point(1,1))*(y_point(1,2)-y_point(1,1));
    for i=1:3
        for j=1:3
            x=x_point(1)*(1-(1-t(i))*(1+t(j))/4-(1+t(i))/2)+x_point(2)*(1-t(i))*(1+t(j))/4+x_point(3)*(1+t(i))/2;
            y=y_point(1)*(1-(1-t(i))*(1+t(j))/4-(1+t(i))/2)+y_point(2)*(1-t(i))*(1+t(j))/4+y_point(3)*(1+t(i))/2;
            f_ij=(1-t(i))/8*(V_x*x*(1+a^2/(x^2+y^2))-1/D*(varphi(ipoin1,1)*(a_ie(1)*x+b_ie(1)*y+c_ie(1))+varphi(ipoin2,1)*(a_ie(2)*x+b_ie(2)*y+c_ie(2))+varphi(ipoin3,1)*(a_ie(3)*x+b_ie(3)*y+c_ie(3))))^2;
            l2errors=l2errors+W(j)*W(i)*f_ij*Jacobi;
        end
    end
 end
 l2errors=sqrt(l2errors);








end
