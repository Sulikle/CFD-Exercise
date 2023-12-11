function X=LUSGS(LHS,R,Unit)
dimention=size(LHS,1)/Unit;
numberx=Unit+1;
%LU-SGS 解三对角矩阵
%取出我们所需要的D
D=zeros(dimention*Unit,dimention*Unit);
for iface=2:numberx
    ieL=iface-1;
    D(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieL-1)+1:dimention*ieL)=LHS(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieL-1)+1:dimention*ieL);
end
%取出我们所需要的L
L=zeros(dimention*Unit,dimention*Unit);
for iface=2:numberx-1
    ieR=iface;
    ieL=iface-1;
    L(dimention*(ieR-1)+1:dimention*ieR,dimention*(ieL-1)+1:dimention*ieL)=LHS(dimention*(ieR-1)+1:dimention*ieR,dimention*(ieL-1)+1:dimention*ieL);
end

%取出我们所需要的U
U=zeros(dimention*Unit,dimention*Unit);
for iface=2:numberx-1
    ieR=iface;
    ieL=iface-1;
    U(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieR-1)+1:dimention*ieR)= LHS(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieR-1)+1:dimention*ieR);
end


X=0;
b=R;%用来存储最初的rhs
    %Forward sweep
    ie=1;
    X(ie:ie+dimention-1,1)=D(ie:ie+dimention-1,ie:ie+dimention-1)\R(ie:ie+dimention-1,1);
    for ie=2:numberx-1
        R(dimention*(ie-1)+1:dimention*ie,1)=R(dimention*(ie-1)+1:dimention*ie,1)-L(dimention*(ie-1)+1:dimention*ie,dimention*(ie-2)+1:dimention*(ie-1))*X(dimention*(ie-2)+1:dimention*(ie-1),1);
        X(dimention*(ie-1)+1:dimention*ie,1)=D(dimention*(ie-1)+1:dimention*ie,dimention*(ie-1)+1:dimention*ie)\R(dimention*(ie-1)+1:dimention*ie,1);
    end
    %Backward sweep
    for ie=1:numberx-1
        R(dimention*(ie-1)+1:dimention*ie,1)= D(dimention*(ie-1)+1:dimention*ie,dimention*(ie-1)+1:dimention*ie)*X(dimention*(ie-1)+1:dimention*ie,1);
    end
    
    ie=numberx-1;
    X(dimention*(ie-1)+1:dimention*ie)=D(dimention*(ie-1)+1:dimention*ie,dimention*(ie-1)+1:dimention*ie)\R(dimention*(ie-1)+1:dimention*ie,1);
    for ie=numberx-2:-1:1
        R(dimention*(ie-1)+1:dimention*ie,1)=R(dimention*(ie-1)+1:dimention*ie,1)-U(dimention*(ie-1)+1:dimention*ie,dimention*(ie)+1:dimention*(ie+1))*X(dimention*(ie)+1:dimention*(ie+1),1);
        X(dimention*(ie-1)+1:dimention*ie,1)=D(dimention*(ie-1)+1:dimention*ie,dimention*(ie-1)+1:dimention*ie)\R(dimention*(ie-1)+1:dimention*ie,1);
    end

end