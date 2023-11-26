function Unext=ExplicitTVDRK3(Ucurrent,Uhold,Deltax,Unit,dimention,R,deltatau,afa,beta,gama)
numberx=Unit+1;
Unext=zeros(dimention,Unit);

if dimention==2
     for k=1:numberx-1
        Mtau=[Deltax(k),0;0,Deltax(k)/12+1/Deltax(k)];
       Unext(:,k)=gama*Ucurrent(:,k)+afa*Uhold(:,k)+Mtau\R(:,k)*deltatau(k)*beta;
    end       
elseif dimention==3

    for k=1:numberx-1
        Mtau=[Deltax(k),0,0;0,Deltax(k)/12+1/Deltax(k),0;0,0,1/720*Deltax(k)+1/(12*Deltax(k))];
       Unext(:,k)=gama*Ucurrent(:,k)+afa*Uhold(:,k)+Mtau\R(:,k)*deltatau(k)*beta;
    end
end  


end

