%function [image] = phantom(xoff,yoff,dx,dy,nx,ny,E)
%Frederic Noo (noo@ucair.med.utah.edu), December 4, 2017

function [image] = phantom(impars,E)

xoff=impars.xoff;
yoff=impars.yoff;
dx=impars.dx;
dy=impars.dy;
nx=impars.nx;
ny=impars.ny;

ox=ones(nx,1);
oy=ones(ny,1);

%create a 2D matrix with elements equal to the x coord. of pixels

xcoord=-xoff*dx+(0:nx-1)*dx;
xcoord=oy*xcoord;

%create a 2D matrix with elements equal to the y coord. of pixels

ycoord=-yoff*dy+(0:ny-1)*dy;
ycoord=((ox*ycoord)');

%initialize the image to zeroes
image=zeros(ny,nx);

for k=1:length(E(:,1)) 

    %create the image of each ellipse
    x0=E(k,1);
    y0=E(k,2);
    a=E(k,3);
    b=E(k,4);
    phi=E(k,5)*pi/180;
    f=E(k,6);

    alpha=1/a;
    beta=1/b; 
    p1=(cos(phi)*alpha)^2+(sin(phi)*beta)^2;
    p2=(sin(phi)*alpha)^2+(cos(phi)*beta)^2;
    p3=sin(phi)*cos(phi)*(alpha^2-beta^2);
  
    equation=p1*(xcoord-x0).^2+p2*(ycoord-y0).^2;
    equation=equation+2*p3*(xcoord-x0).*(ycoord-y0);

    i=find(equation<1.0);
    ell=zeros(size(equation));
    ell(i)=1;
    image=image+f*ell;

end;

