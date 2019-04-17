function [Ssc,Sext,Cscat,Cext,Qscat,Qext,P1P2,P1,P2,theta,P,g,HGPF1,HGPF2,g2]=MieFunction(a,m,lambda,k,AngLim,DTH,Density)
x=k*a;
Nc=round(x+4*x^(1/3)+2); %covergence number
Nd=2*Nc;

%** Radial Declarations **%
theta=DTH:DTH:(AngLim*pi/180);
Ntheta=length(theta);
u=cos(theta)';

% 1 is Pi0 and 2 ps Pi1
PI(:,1)=zeros(Ntheta,1);
PI(:,2)=ones(Ntheta,1);
TAU(:,1)=zeros(Ntheta,1);
TAU(:,2)=u;

% starts at -1
Ex(1)=cos(x)-i*sin(x);
Ex(2)=sin(x)+i*cos(x);

% Calculation of PI and TAU
for n=3:Nc %starts at 3 but it is n=2 
    N=n-1;
    PI(:,n)=((2*N-1)/(N-1))*u.*PI(:,n-1)-(N/(N-1))*PI(:,n-2);
    TAU(:,n)=u.*(PI(:,n)-PI(:,n-2))-(2*N-1)*(1-u.^2).*PI(:,n-1)+TAU(:,n-2);
end

% Calculation of E(x) or Zeta Function in A and B
for n=3:2*Nc
    N=n-2;
    Ex(n)=((2*N-1)/x)*Ex(n-1)-Ex(n-2);
end
Zeta=Ex(2:Nc+1);

% Calculation of Log Derivative Backward Recursion
Dmx(Nd)=0;
for n=Nd:-1:2
    N=n-1;
    Dmx(n-1)=N/(m*x)-1/(Dmx(n)+N/(m*x));
end
Dmx=Dmx(1:Nc); 

S1=0; S2=0; A(1)=0; B(1)=0; Ssc=0; Sext=0; MeanCOSsum=0;
% Calculation of the Coefficients A and B
for n=2:Nc
    N=n-1;
    A(n)= ((Dmx(n)/m+N/x)*real(Zeta(n))-real(Zeta(n-1)))/ ...
        ((Dmx(n)/m+N/x)*Zeta(n)-Zeta(n-1));
    
    B(n)= ((m*Dmx(n)+N/x)*real(Zeta(n))-real(Zeta(n-1)))/ ...
        ((m*Dmx(n)+N/x)*Zeta(n)-Zeta(n-1));
    
    S1=S1+((2*N+1)/(N*(N+1)))*(A(n)*PI(:,n)+B(n)*TAU(:,n));
    S2=S2+((2*N+1)/(N*(N+1)))*(B(n)*PI(:,n)+A(n)*TAU(:,n));
    Ssc=Ssc+(2*N+1)*(abs(A(n))^2+abs(B(n))^2);
    Sext=Sext+(2*N+1)*(real(A(n)+B(n)));
end

%** Calculating Other Useful Information **%
Vol=(4/3)*pi*a^3;
Area=pi*a^2;
Cscat=(2*pi/(k^2))*Ssc;
Cext=(2*pi/(k^2))*Sext;
Qscat=Cscat/(Area); %Qscat=(2/x^2)*Sscat;
Qext=Cext/(Area); %Qext=(2/x^2)*Sext;

% %** Other things not used **%
% Qabs=Qext-Qscat;
% I180=((x/2)^2)*((real(m)-1)^2+imag(m)^2)/((real(m)+1)^2+imag(m)^2);
% QeffBest=2-(7.680/x)*sin(0.684*x)+1.84*x^(-2/3);


%** Phase functions **%
Normalization=1;
P1P2=((Normalization)*(abs(S1).^2+abs(S2).^2))/2; %UnPolarized 0.5*S+0.5*P
P1=((Normalization)*abs(S1).^2);               %Perpendicular S - Polarization
P2=((Normalization)*abs(S2).^2);             %Parallel      P - Polarization

%** Degree of Polarization **%
P=(abs(S1).^2-abs(S2).^2)./(abs(S1).^2+abs(S2).^2);

%** Asymmetry Factor (mean COSINE) or g **%
for n=2:Nc-1
    N=n;
    MeanCOSsum=MeanCOSsum+((N*(N+2))/(N+1))*real(A(n)*A(n+1)'+B(n)*B(n+1)')+...
        ((2*N+1)/(N*(N+1)))*real(A(n)*B(n)');
end
g=(4/(x^2*Qscat))*MeanCOSsum;

%** Henyey-Greenstein Scatteirng phase function
g1=g;
BB=20;
HGPF1=(1-g1^2)./(1+g1^2-2*g1*cos(theta)).^(3/2);
G=((10*double(eulergamma)^2+lambda*BB)/(8*pi))*(8/(10+(5*x)/(double(eulergamma)*pi)^2));
g2=cos(G)^2/(1+cos(G));
HGPF2=(1-g2^2)./(1+g2^2-2*g2*cos(theta)).^(3/2);

end
