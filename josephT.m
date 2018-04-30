%Transpose operation for Joseph's program for estimation of sinogram from a 2D image.
%The output image is im[ny][nx], with y increasing from top to bottom, and x
%increasing from left to right.
%function[im]= josephT(sino,xoff,yoff,dx,dy,nx,ny,soff,ds,ns,nth)

function[im]= josephT(sino,impars,datapars)

xoff=impars.xoff;
yoff=impars.yoff;
dx=impars.dx;
dy=impars.dy;
nx=impars.nx;
ny=impars.ny;

soff=datapars.soff;
ds=datapars.ds;
ns=datapars.ns;
nth=datapars.nth;

theta=(0:nth-1)*pi/nth;
costheta=cos(theta);
sintheta=sin(theta);
tantheta=tan(theta);
cotantheta=1./tantheta;

indexcos=find(abs(costheta) > 1/sqrt(2));
tmp=ones(size(theta));
tmp(indexcos)=0;
indexsin=find(tmp>0);

sval=((0:ns-1)-soff)*ds;
xval=((0:nx-1)-xoff)*dx;
yval=((0:ny-1)-yoff)*dy;

im=zeros(ny,nx);

%sinogram part corresponding to indexcos
th=theta(indexcos);
costh=costheta(indexcos);
tanth=tantheta(indexcos);
factor=dy./abs(costh);factor=factor'*ones(1,ns);
tmp=sino(indexcos,:).*factor;
for j=1:ny

   scoord=ones(length(th),1)*sval;
   thcoord=th'*ones(1,ns);
   costhcoord=costh'*ones(1,ns);
   tanthcoord=tanth'*ones(1,ns);

   y=yval(j);
   x=scoord./costhcoord-tanthcoord*y;
   x=x/dx+xoff;
   ix=floor(x);
   wint=x-ix;
   ix=ix+1;

   i=find( (ix>0) & (ix<nx) );

   p=zeros(nx,1);
   for k=1:length(i)
       ik=i(k);
       p(ix(ik))=p(ix(ik))+tmp(ik)*(1-wint(ik));
       p(ix(ik)+1)=p(ix(ik)+1)+tmp(ik)*wint(ik);
   end
   im(j,:)=p';
   
end

%second part of the sinogram, corresponding to indexsin
th=theta(indexsin);
sinth=sintheta(indexsin);
cotanth=cotantheta(indexsin);
factor=dx./abs(sinth);factor=factor'*ones(1,ns);
tmp=sino(indexsin,:).*factor;
for j=1:nx

   scoord=ones(length(th),1)*sval;
   thcoord=th'*ones(1,ns);
   sinthcoord=sinth'*ones(1,ns);
   cotanthcoord=cotanth'*ones(1,ns);

   x=xval(j);
   y=scoord./sinthcoord-cotanthcoord*x;
   y=y/dy+yoff;
   iy=floor(y);
   wint=y-iy;
   iy=iy+1;

   i=find(iy>0 & iy<ny);

   p=zeros(ny,1);
   for k=1:length(i)
       ik=i(k);
       p(iy(ik))=p(iy(ik))+tmp(ik)*(1-wint(ik));
       p(iy(ik)+1)=p(iy(ik)+1)+tmp(ik)*wint(ik);
   end
   im(:,j)=im(:,j)+p;
end



