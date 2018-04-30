function [rcn] = GradDescent(impars,datapars,regpars,initim,data,lambda,niter,fovmask)
 
rcn=fovmask.*initim;

for n=1:niter
    residual = joseph(rcn,impars,datapars)-data;
    iminc = josephT(residual,impars,datapars);
    iminc = iminc+regpars.beta*pengrad(impars,regpars.mode,regpars.delta,rcn);
    rcn=fovmask.*(rcn-lambda*iminc);
    if regpars.pos == 1
        rcn(rcn<0)=0;
    end
    imagesc(flipud(rcn),[0.85,1.15]),axis image,title(num2str(n)),
    colorbar,colormap gray,pause(1);
end

end



