function [rcnp] = AccGradDescent(impars,datapars,regpars,initim,data,lambda,niter,fovmask)
 
told=1;
rcnqold=fovmask.*initim;
if regpars.pos == 1
   rcnqold(rcnqold<0)=0;
end
rcnp=rcnqold;

for n=1:niter
    residual = joseph(rcnp,impars,datapars)-data;
    iminc = josephT(residual,impars,datapars);
    iminc = iminc+regpars.beta*pengrad(impars,regpars.mode,regpars.delta,rcnp);
    rcnqnew=fovmask.*(rcnp-lambda*iminc);
    tnew=0.5*(1+sqrt(1+4*told^2));
    rcnp=rcnqnew+(told-1)/tnew*(rcnqnew-rcnqold);
    if regpars.pos == 1
        rcnp(rcnp<0)=0;
    end
    rcnqold=rcnqnew;
    told=tnew;
    imagesc(flipud(rcnp),[0.85,1.15]),axis image,title(num2str(n)),
    colorbar,colormap gray,pause(1);
end

end



