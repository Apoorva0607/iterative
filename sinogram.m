%function [sino] = sinogram(soff,ds,ns,nth,E,rot,mu)
%simulation of parallel-beam data 
%rot = 1 for data over 180 degrees
%    = 2 for data over 360 degrees
%Frederic Noo (noo@ucair.med.utah.edu), December 4, 2017.

function [sino] = sinogram(datapars,E,rot)

soff=datapars.soff;
ds=datapars.ds;
ns=datapars.ns;
nth=datapars.nth;

scoord=((0:ns-1)-soff)*ds;

th=(0:nth-1)*rot*pi/nth;
coth=cos(th);
sith=sin(th);

sino=zeros(nth,ns);

for k=1:length(E(:,1)) 

    x0=E(k,1);
    y0=E(k,2);
    a=E(k,3);
    b=E(k,4);
    phi=E(k,5)*pi/180;
    f=E(k,6);

    alpha=1/a;
    beta=1/b; 

    sino_ell=zeros(nth,ns);

    for ith=1:nth; %create each projection successively
        angle=th(ith)-phi;
        co=cos(angle);
        si=sin(angle);
        s1=scoord-x0*coth(ith)-y0*sith(ith);        
        s2=s1*alpha*beta/sqrt((si*alpha)^2+(co*beta)^2);        
        c=1/sqrt((si*alpha)^2+(co*beta)^2);
        
        i=find(abs(s2) <= 1);    
        sino_ell(ith,i)=2*f*c*sqrt(1-s2(i).*s2(i)); 
    end;
    sino=sino+sino_ell;

end;

