clc
clear
close all
%% Preproceeding
%Input data
addpath(genpath('/Users/mac/Documents/CFD-Exercise-main/CFDProject1/Proj1_Mesh'))
% store the nodal points corresponding to each element in Edata————INPOEL
feflo_bump = fopen('Proj1_Mesh/feflo.domn.bump', 'rb');%open the file
skiprows=fscanf(feflo_bump, '%lf',1);% skip some rows
for i = 1:skiprows+2
    line=fgetl(feflo_bump);
end
Ndata = fscanf(feflo_bump, '%lf',[2,1]);
Ndim=Ndata(1,1);Nnode=Ndata(2,1);
for i = 1:2
    line=fgetl(feflo_bump);
end
Ndata = fscanf(feflo_bump, '%lf',[3,1]);
Nelement=Ndata(1,1);Npoint=Ndata(2,1);Nface=Ndata(3,1);%Information of elements&points&boundary
for i = 1:3
    line=fgetl(feflo_bump);
end

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
C = textscan(line, '%f');
numbers = length(str2num(line));
BCOND=fscanf(feflo_bump, '%lf',[numbers,Nface]);
BCOND=[C{1,1},BCOND];
BCOND=BCOND(2:4,:);%每一列存储的是该边界处的左右点以及flag，2代表固定面，4代表无穷远处
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
v_element=zeros(Ndim,Nelement);%v in element
v_point=zeros(Ndim,Npoint);%v in point
v_point_scalar=zeros(1,Npoint);%store the size of v
v_point_scalar_face=zeros(1,Npoint);%store the size of v in boundary



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
% %solve the equation 
% %SGS速度比较慢，故不采用
% endtimes=10000;tol=10^(-9);epsilon=10^(-12);
% A=LHS;b=RHS;b1=b;varphi=zeros(Npoint,1);varphi_start=varphi;%initial value
% resnorm_start=norm(A*varphi_start-b)+epsilon;
% for times=2:endtimes
%     %calculate X Forward sweep
%     for i=1:Npoint
%         for j=1:i-1
%             b(i)=b(i)-A(i,j)*varphi(j);
%         end
%         for j=i+1:Npoint
%             b(i)=b(i)-A(i,j)*varphi(j);
%         end
%         varphi(i)=b(i)/A(i,i);
%     end
%      b=b1;
%      %Backward sweep
%     for i=Npoint:-1:1
%         for j=Npoint:-1:i+1
%             b(i)=b(i)-A(i,j)*varphi(j);
%         end
%         for j=i-1:-1:1
%             b(i)=b(i)-A(i,j)*varphi(j);
%         end
%         varphi(i)=b(i)/A(i,i);
%     end
%     b=b1;
%     resnorm=norm(A*varphi-b);
%     ratio=resnorm/resnorm_start;
%     if ratio<tol
%         break;
%     end
% end

[l,u]=lu(LHS);
y=l\RHS;
varphi=u\y;

%solve v in element
for ie=1:Nelement
    ipoin1=INPOEL(1,ie);ipoin2=INPOEL(2,ie);ipoin3=INPOEL(3,ie);
    HPcycle=[ipoin1,ipoin2,ipoin3,ipoin1,ipoin2];c_ie=zeros(1,3);b_i=0;a_i=0;
    for k=1:Nnode%先计算c
           c_ie(k)=COORD(1,HPcycle(k+1))*COORD(2,HPcycle(k+2))-COORD(1,HPcycle(k+2))*COORD(2,HPcycle(k+1));
    end
    for i=1:Nnode
        a_i=COORD(2,HPcycle(i+1))-COORD(2,HPcycle(i+2));
        b_i=-(COORD(1,HPcycle(i+1))-COORD(1,HPcycle(i+2)));
        v_element(1,ie)=v_element(1,ie)+varphi(INPOEL(i,ie))*a_i/sum(c_ie);
        v_element(2,ie)=v_element(2,ie)+varphi(INPOEL(i,ie))*b_i/sum(c_ie);
    end
end

%solve v in point
for IPOIN=1:Npoint %loop over the points
    sum_weight_x=0;sum_weight_y=0;
    NESP=ESUP2(1,IPOIN+1)-ESUP2(1,IPOIN);%number of elements surrounding IPOIN
    WEIGHT=zeros(1,NESP);
    for ie=ESUP2(IPOIN)+1:ESUP2(IPOIN+1)%loop over the elements surrounding IPOIN
    ipoin1=INPOEL(1,ESUP1(ie));ipoin2=INPOEL(2,ESUP1(ie));ipoin3=INPOEL(3,ESUP1(ie));
    HPcycle=[ipoin1,ipoin2,ipoin3,ipoin1,ipoin2];
    for k=1:Nnode%计算第ie-ESUP2(IPOIN)单元的D
           WEIGHT(ie-ESUP2(IPOIN))=WEIGHT(ie-ESUP2(IPOIN))+COORD(1,HPcycle(k+1))*COORD(2,HPcycle(k+2))-COORD(1,HPcycle(k+2))*COORD(2,HPcycle(k+1));
    end
    end
    for ie=ESUP2(IPOIN)+1:ESUP2(IPOIN+1)
        sum_weight_x=sum_weight_x+v_element(1,ESUP1(ie))*WEIGHT(ie-ESUP2(IPOIN));
        sum_weight_y=sum_weight_y+v_element(2,ESUP1(ie))*WEIGHT(ie-ESUP2(IPOIN));        
    end
    v_point(1,IPOIN)=sum_weight_x/sum(WEIGHT);v_point(2,IPOIN)=sum_weight_y/sum(WEIGHT);
end
%calculate the size of v in points
for IPOIN=1:Npoint
    v_point_scalar(1,IPOIN)=sqrt(v_point(1,IPOIN)^2+v_point(2,IPOIN)^2);
end

%calculate the size of v in boundarys
for IFACE=1:Nface  
    ipoin1=BCOND(1,IFACE);ipoin2=BCOND(2,IFACE);flag=BCOND(3,IFACE);
    if flag==2
    v_point_scalar_face(1,ipoin1)=sqrt(v_point(1,ipoin1)^2+v_point(2,ipoin1)^2);
    v_point_scalar_face(1,ipoin2)=sqrt(v_point(1,ipoin2)^2+v_point(2,ipoin2)^2);
    end
end

%% Postproceeding
%导出网格
% 打开文件以进行写入
INPOEL=INPOEL';
DATA_mix=[COORD',v_point',varphi,v_point_scalar',v_point_scalar_face'];
DATA = fopen('data.dat', 'w');
% 使用 fprintf 将数据写入文件
fprintf(DATA,'Variables=x,y,u,v,Potential,Velocity,Velocity_B\n');
fprintf(DATA,'Zone n=');fprintf(DATA,num2str(Npoint));fprintf(DATA,',e=');fprintf(DATA,num2str(Nelement));fprintf(DATA,',f=fepoint,et=triangle\n');
fprintf(DATA,'%f\t%f\t%f\t%f\t%f\t%f\t%f\n',DATA_mix');
fprintf(DATA,'%d\t%d\t%d\n',INPOEL');
% 关闭文件
fclose(DATA);

%The second method constructing load vector(有瑕疵)
%construct load vector
% for iface=1:Nface
%     ip1=BCOND(1,iface);ip2=BCOND(2,iface);flag=BCOND(3,iface);%信息提取
%     gama_e=sqrt((COORD(1,ip1)-COORD(1,ip2))^2+(COORD(2,ip1)-COORD(2,ip2))^2);%length of face
%     if flag==4&&(abs(COORD(1,ip1))<1e-5&&abs(COORD(1,ip2))<1e-5)%left face
%         g=-1;
%     elseif flag==4&&(abs(COORD(1,ip1)-3)<1e-5&&abs(COORD(1,ip2)-3)<1e-5)%right face
%         g=1;
%     end
%     if flag==2&&(abs(COORD(2,ip1)-1)<1e-5&&abs(COORD(2,ip2)-1)<1e-5)%upper face
%         g=0;
%     elseif flag==2&&(abs(COORD(2,ip1))<1e-5&&abs(COORD(2,ip2))<1e-5)%lowwer face
%         g=0;
%     elseif flag==2&&(abs(COORD(2,ip1))>1e-5&&abs(COORD(2,ip2))>1e-5)&&(abs(COORD(2,ip1)-1)>1e-5&&abs(COORD(2,ip2)-1)>1e-5)%圆弧处
%         if COORD(1,ip2)>COORD(1,ip1)%判断点的前后
%             V_a=COORD(1,ip2)-COORD(1,ip1);V_b=COORD(2,ip2)-COORD(2,ip1);
%             if V_b>0
%                 V_x=sqrt(V_b^2/(V_a^2+V_b^2));V_y=-V_a/V_b*V_x;
%             elseif V_b<0
%                 V_x=-sqrt(V_b^2/(V_a^2+V_b^2));V_y=-V_a/V_b*V_x;
%             elseif V_b==0
%                 V_x=0;
%             end
%         elseif COORD(1,ip2)<COORD(1,ip1)
%             V_a=COORD(1,ip1)-COORD(1,ip2);V_b=COORD(2,ip1)-COORD(2,ip2);
%             if V_b>0
%                 V_x=sqrt(V_b^2/(V_a^2+V_b^2));V_y=-V_a/V_b*V_x;
%             elseif V_b<0
%                 V_x=-sqrt(V_b^2/(V_a^2+V_b^2));V_y=-V_a/V_b*V_x;
%             elseif V_b==0
%                 V_x=0;
%             end
%         end
%         g=V_x;
%     end
%     RHS(ip1)=RHS(ip1)+0.5*g*gama_e;
%     RHS(ip2)=RHS(ip2)+0.5*g*gama_e;
% end





