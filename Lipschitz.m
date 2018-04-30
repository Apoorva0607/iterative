function [mu] = Lipschitz(impars,datapars,niter,fovmask)

u=ones(impars.ny,impars.nx).*fovmask;
u=u/norm(u(:));
for k=1:niter
    subz=joseph(u,impars,datapars);
    z=fovmask.*josephT(subz,impars,datapars);
    mu=u(:)'*z(:);
    u=z/norm(z(:));
end






