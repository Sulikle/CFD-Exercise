function [n,UL2errors,Vl2errors]=Judges(Unit,CFL,endtau,tol,belta,norder,nreconstruct,nexplicit,nsdv,omega0,omega1,omega2,omega3,omegab)

if norder==0&&nreconstruct==0%DG(P0P1)+DG(P0)
   [n,UL2errors,Vl2errors]=subDGP0P1plusDGP0(Unit,CFL,endtau,tol,belta,nsdv,nexplicit);

elseif norder==1&&nreconstruct==0%DG(P0P2)+DG(P1)
   [n,UL2errors,Vl2errors]=subDGP0P2plusDGP1(Unit,CFL,endtau,tol,belta,nsdv,nexplicit);
   
elseif norder==0&&nreconstruct==1%DG(P0P2)+rDG(P0P1)
    [n,UL2errors,Vl2errors]=subDGP0P2plusrDGP0P1(Unit,CFL,endtau,tol,belta,nsdv,nexplicit,omega0,omega1,omega2,omega3,omegab);
    
    elseif norder==0&&nreconstruct==2%DG(P0P3)+rDG(P0P2)
    [n,UL2errors,Vl2errors]=subDGP0P3plusrDGP0P2(Unit,CFL,endtau,tol,belta,nsdv,nexplicit,omega0,omega1,omega2,omega3,omegab);
    
          elseif norder==1&&nreconstruct==1%DG(P0P3)+rDG(P1P2)
    [n,UL2errors,Vl2errors]=subDGP0P3plusrDGP1P2(Unit,CFL,endtau,tol,belta,nsdv,nexplicit,omega0,omega1,omega2,omega3,omegab);
    
else
    
    if nreconstruct==0
        fprintf('DG(P0P%d)+DG(P%d)这种情况暂时还没考虑呢！前面的区域请以后再来探索吧!\n',norder+nreconstruct+1,norder);
    n=nan;UL2errors=nan;Vl2errors=nan;
    else
         fprintf('DG(P0P%d)+rDG(P%dP%d)这种情况暂时还没考虑呢！前面的区域请以后再来探索吧!\n',norder+nreconstruct+1,norder,norder+nreconstruct);
    n=nan;UL2errors=nan;Vl2errors=nan;
    end
    
end


end
