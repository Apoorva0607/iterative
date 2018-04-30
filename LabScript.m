%%
%data simulation
datapars.nth=60;
datapars.ns=100;
datapars.ds=1.01;
datapars.soff=50.25;

smallph;
nsub=10;
soff=datapars.soff;
sino=zeros(datapars.nth,datapars.ns);
for k=1:nsub 
    shift=0.5-1/(2*nsub)-(k-1)/nsub;
    datapars.soff=soff+shift;
    sino=sino+sinogram(datapars,E,1);
end;
sino=sino/nsub;
datapars.soff=soff;

%with noise
noise=1;           %INVESTIGATE noise=0 and noise=1
if noise
   muw=0.0183;
   Nin=6e4;
   Nout=Nin*exp(-muw*sino);
   sino=(muw*sino+randn(size(sino))./sqrt(Nout))/muw;
end

%%
%image reconstruction set-up
impars.nx=100;
impars.dx=1;
impars.xoff=50;
impars.ny=100;
impars.dy=1;
impars.yoff=50;
radius=48;

fovmask=phantom(impars,[0 0 radius radius 0 1]);
fovmask(1,:)=0;fovmask(end,:)=0;
fovmask(:,1)=0;fovmask(:,end)=0;

truth=phantom(impars,E);
figure(1),imagesc(truth,[0.85,1.15]),axis image,colormap gray

%%
%Lipschitz constant evaluation
niter=5;
[L] = Lipschitz(impars,datapars,niter,fovmask) %CREATE THIS FUNCTION

%%
%iterative reconstruction
lambda= 0.95*1/L;  %ADJUST
niter=  100;       %ADJUST

initim=zeros(impars.ny,impars.nx);

%reconstruction with quadratic penalty
regpars.pos=1;    %1: with positivy constraint, else without
regpars.mode=0;   %1: edge preserving, else quadratic INVESTIGATE BOTH OPTIONS
regpars.beta=2.0; %INVESTIGATE USING 0.0 and 2.0
regpars.delta=0.005;

%CREATE THE TWO FUNCTIONS BELOW AS ALTERNATIVE OPTIONS
%DISPLAY THE RESULT IN figure(1) AT EACH ITERATION WITHIN THE FUNCTION

figure(1),
%[rcn1] = GradDescent(impars,datapars,regpars,initim,sino,lambda,niter,fovmask);
[rcn2] = AccGradDescent(impars,datapars,regpars,initim,sino,lambda,niter,fovmask);







