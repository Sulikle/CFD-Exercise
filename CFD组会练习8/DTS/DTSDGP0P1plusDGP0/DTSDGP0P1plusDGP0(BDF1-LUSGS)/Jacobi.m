function X=Jacobi(LHS,R,Unit)
dimention=size(LHS,1)/Unit;
numberx=Unit+1;
%取出我们所需要的D
D=zeros(dimention*Unit,dimention*Unit);
for iface=2:numberx
    ieL=iface-1;
    D(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieL-1)+1:dimention*ieL)=LHS(dimention*(ieL-1)+1:dimention*ieL,dimention*(ieL-1)+1:dimention*ieL);
end

X=D\R;

end