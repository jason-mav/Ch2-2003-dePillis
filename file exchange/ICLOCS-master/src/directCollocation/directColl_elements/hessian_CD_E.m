function [ Ezz ] = hessian_CD_E( Ezz, E, x0, xf, u0, uf, p, t0, tf, e, e2, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if data.FD.FcnTypes.Etype
    idx=data.FD.index.Ey;nfd=size(idx,2);                               
    et0=data.FD.vector.Ey.et0;etf=data.FD.vector.Ey.etf;ep=data.FD.vector.Ey.ep;
    ex0=data.FD.vector.Ey.ex0;eu0=data.FD.vector.Ey.eu0;
    exf=data.FD.vector.Ey.exf;euf=data.FD.vector.Ey.euf;ez=e*data.FD.vector.Ey.ez;

    for i=1:nfd
       for j=1:i
        if j==i
         Ep1=E(x0+ex0(:,i)*ez(i),xf+exf(:,i)*ez(i),u0+eu0(:,i)*ez(i),uf+euf(:,i)*ez(i),p+ep(:,i)*ez(i),t0+et0(:,i),tf+etf(:,i),vdat);
         Eo=E(x0,xf,u0,uf,p,t0,tf,vdat);
         Ep2=E(x0-ex0(:,i)*ez(i),xf-exf(:,i)*ez(i),u0-eu0(:,i)*ez(i),uf-euf(:,i)*ez(i),p-ep(:,i)*ez(i),t0-et0(:,i),tf-etf(:,i),vdat);
         Ezz(idx(i),idx(j))=(Ep1-2*Eo+Ep2)/e2; 
        else
         Epp=E(x0+ex0(:,i)*ez(i)+ex0(:,j)*ez(i),xf+exf(:,i)*ez(i)+exf(:,j)*ez(i),u0+eu0(:,i)*ez(i)+eu0(:,j)*ez(i),...
              uf+euf(:,i)*ez(i)+euf(:,j)*ez(i),p+ep(:,i)*ez(i)+ep(:,j)*ez(i),t0+et0(:,i)+et0(:,j),tf+etf(:,i)+etf(:,j),vdat);
         Epm=E(x0+ex0(:,i)*ez(i)-ex0(:,j)*ez(i),xf+exf(:,i)*ez(i)-exf(:,j)*ez(i),u0+eu0(:,i)*ez(i)-eu0(:,j)*ez(i),...
              uf+euf(:,i)*ez(i)-euf(:,j)*ez(i),p+ep(:,i)*ez(i)-ep(:,j)*ez(i),t0+et0(:,i)-et0(:,j),tf+etf(:,i)-etf(:,j),vdat);

         Emp=E(x0-ex0(:,i)*ez(i)+ex0(:,j)*ez(i),xf-exf(:,i)*ez(i)+exf(:,j)*ez(i),u0-eu0(:,i)*ez(i)+eu0(:,j)*ez(i),...
              uf-euf(:,i)*ez(i)+euf(:,j)*ez(i),p-ep(:,i)*ez(i)+ep(:,j)*ez(i),t0-et0(:,i)+et0(:,j),tf-etf(:,i)+etf(:,j),vdat);
         Emm=E(x0-ex0(:,i)*ez(i)-ex0(:,j)*ez(i),xf-exf(:,i)*ez(i)-exf(:,j)*ez(i),u0-eu0(:,i)*ez(i)-eu0(:,j)*ez(i),...
              uf-euf(:,i)*ez(i)-euf(:,j)*ez(i),p-ep(:,i)*ez(i)-ep(:,j)*ez(i),t0-et0(:,i)-et0(:,j),tf-etf(:,i)-etf(:,j),vdat);

        Ezz(idx(i),idx(j))=(Epp+Emm-Epm-Emp)/e2/4;  
        end
      end
    end
end

end

