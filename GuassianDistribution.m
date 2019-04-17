function [Radii,Weights,Distribution,NX]=GuassianDistribution(u,sigma,Resolution,NumofSTDEV,k)
checkLOW=u-NumofSTDEV*sigma;
checkUP=u+NumofSTDEV*sigma;
x=0:1/Resolution:3*u;
X=k*x;
DIST=(1/sqrt(2*pi*sigma^2))*exp((-(x-u).^2)/(2*sigma^2));
DISTNX=(1/sqrt(2*pi*sigma^2))*exp((-(X-u).^2)/(2*sigma^2));
count=1;
for n=1:length(DIST)
        if x(n) <= checkUP && x(n) >=checkLOW
            boundPDF(count)=DIST(n);
            boundNX(count)=DISTNX(n);
            %nPDF(count)=boundPDF(count)/max(DIST);
            Radii(count)=x(n);
            count=count+1;
        end
end
NX=boundNX;
Weights=boundPDF;
Distribution=DIST;
end

