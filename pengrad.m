function [g] = pengrad(impars,mode,delta,im)

g=zeros(impars.ny,impars.nx);
for i=2:impars.ny-1
    for j=2:impars.nx-1
        g1=psiprime(im(i,j+1)-im(i,j),mode,delta);
        g2=psiprime(im(i,j)-im(i,j-1),mode,delta);
        g3=psiprime(im(i+1,j)-im(i,j),mode,delta);
        g4=psiprime(im(i,j)-im(i-1,j),mode,delta);
        g(i,j)=g2-g1+g4-g3;
    end
end

end

function [psiprimeval] = psiprime(t,mode,delta)

switch mode
    case 1
        psiprimeval=t/sqrt(t^2+delta^2);
    otherwise
        psiprimeval=t/delta;
end

end





