close all;clear;clc
tic
% This code features a polar plotting function that was influenced by 
% Mazin Mustafa which can be found at: 
% Mazin Mustafa (2019). Polar_dB 
%(https://www.mathworks.com/matlabcentral/fileexchange/67721-polar_db), 
% MATLAB Central File Exchange.
%** Variable Inputs **%
AngLim = 180;             %Display Angle Limit in Degrees
N0 = 1.000;                   %Index of medium
N1 = 1.335;1.00029944-0i; 1.3290-2.93e-7i;      %Ind0ex of Particles
lambda = .850;            %wavelegth in um

D = 10;.64e-4; %Particle Diameter in um
a = D/2; 

m = N1/N0;    %relative index of refraction
k = (2*pi*N0)/(lambda); %wave number in medium
x = k*a; %Size Factor
Tres = 0.001; %Theta Resolution
%*************************************************************************%
Nc = round(x+4*x^(1/3)+2)+1; %covergence number
Nd = 2*Nc;
Density =  1;2.5455e25;
DePol = 0.035;

%** Radial Declarations **%
DTH=Tres; theta=(DTH:DTH:(AngLim*pi/180));
Ntheta=length(theta);
u=cos(theta)';

%** Computations **%
% 1 is Pi0 and 2 ps Pi1
PI(:,1) = zeros(Ntheta,1);
PI(:,2) = ones(Ntheta,1);
TAU(:,1) = zeros(Ntheta,1);
TAU(:,2) = u;

% starts at -1
Ex(1) = cos(x)-1i*sin(x);
Ex(2) = sin(x)+1i*cos(x);

% Calculation of PI and TAU
for n = 3:Nc %starts at 3 but it is n=2 
    N = n-1;
    PI(:,n) = ((2*N-1)/(N-1))*u.*PI(:,n-1)-(N/(N-1))*PI(:,n-2);
    TAU(:,n) = u.*(PI(:,n)-PI(:,n-2))-(2*N-1)*(1-u.^2).*PI(:,n-1)+TAU(:,n-2);
end

% Calculation of E(x) or Zeta Function in A and B
for n = 3:2*Nc %starts at 3 but n=1 1(-1) 2(0) 3(1)
    N = n-2;
    Ex(n) = ((2*N-1)/x)*Ex(n-1)-Ex(n-2);
end
Zeta = Ex(2:Nc+1);

% Calculation of Log Derivative Backward Recursion
Dmx(Nd) = 0;
for n = Nd:-1:2
    N = n-1;
    Dmx(n-1) = N/(m*x)-1/(Dmx(n)+N/(m*x));
end

Dmx = (Dmx(1:Nc)); S1=0; S2=0; A(1)=0; B(1)=0; 
Sscat = 0; Sext = 0; count=0; MeanCOSsum=0;

% Calculation of the Coefficients A and B
for n=2:Nc
    N = n-1;
    A(n) = ((Dmx(n)/m+N/x)*real(Zeta(n))-real(Zeta(n-1)))/ ...
        ((Dmx(n)/m+N/x)*Zeta(n)-Zeta(n-1));
    
    B(n) = ((m*Dmx(n)+N/x)*real(Zeta(n))-real(Zeta(n-1)))/ ...
        ((m*Dmx(n)+N/x)*Zeta(n)-Zeta(n-1));
    
    S1 = S1+((2*N+1)/(N*(N+1)))*(A(n)*PI(:,n)+B(n)*TAU(:,n));
    S2 = S2+((2*N+1)/(N*(N+1)))*(B(n)*PI(:,n)+A(n)*TAU(:,n));
    Sscat = Sscat+(2*N+1)*(abs(A(n))^2+abs(B(n))^2);
    Sext = Sext+(2*N+1)*(real(A(n))+real(B(n)));
%     if (abs(A(n))^2+abs(B(n))^2) < 1e-14
%         break;
%     end
%     count=count+1;
end

%** Calculating Other Useful Information **%
Vol = (4/3)*pi*a^3;
Area = pi*a^2;
Cscat = (2*pi/(k^2))*Sscat;
Cext = (2*pi/(k^2))*Sext;
Qscat = Cscat/(Area); %Qscat=(2/x^2)*Sscat;
Qext = Cext/(Area); %Qext2=(2/x^2)*Sext
Qabs = Qext-Qscat;
ExtCoeR1 = Density*Cext*1e-8*1e2;
ExtCoeR2 = Density*Qext*Area*1e-8*1e2; %Another way to calculate \sigma 
ExtCoe = Density*(pi/k^2)*Qext*x^2*1e-8*1e2;

%** Degree of Polarization **%
P=(abs(S1).^2-abs(S2).^2)./(abs(S1).^2+abs(S2).^2);
DegPol=abs(P);

%** Phase functions **%
Normalization=(4*pi)/(k^2*ExtCoeR1*1e8*1e-2); % or 2*x^2 2*pi*Sscat
P1P2=.5*((Normalization)*(abs(S1).^2+abs(S2).^2)); %UnPolarized 0.5*S+0.5*P
P1=((Normalization)*abs(S1).^2);               %Perpendicular S - Polarization
P2=((Normalization)*abs(S2).^2);             %Parallel      P - Polarization

%** Asymmetry Factor (mean COSINE) or g **%
for n=2:Nc-1
    N=n;
    MeanCOSsum=MeanCOSsum+((N*(N+2))/(N+1))*real(A(n)*A(n+1)+B(n)*B(n+1))+...
        ((2*N+1)/(N*(N+1)))*real(A(n)*B(n));
end
MeanCOS=(4/(x^2*Qscat))*MeanCOSsum;
MeanCOS2 = simpson(theta,P1P2'.*cos(theta).*sin(theta));

%** Henyey-Greenstein Scatteirng phase function
g1=MeanCOS;
BB=20;
HGPF1=(1/(4*pi))*((1-g1^2)./(1+g1^2-2*g1*cos(theta)).^(3/2));
G=((10*double(eulergamma)^2+lambda*BB)/(8*pi))*(8/(10+(5*x)/(double(eulergamma)*pi)^2));
g2=cos(G)^2/(1+cos(G));
HGPF2=(1/(4*pi))*((1-g2^2)./(1+g2^2-2*g2*cos(theta)).^(3/2));

%** Plotting **%
figure(1);
semilogy( theta*180/pi , P1P2/(4*pi)'); hold on;
semilogy( theta*180/pi , P1/(4*pi)');
semilogy( theta*180/pi , P2/(4*pi)');

title(['Normalized Phase Function Through Air r = ' num2str(a) '\mum at ' num2str(lambda*1000) ' nm '])
xlabel('Angle [Degrees]')
ylabel('Relative Intensity []')
legend('Unpolarized','Perp','Parallel','Henyey-Greenstein','Iout')

figure(5); 
PolarPlotdB(theta,10*log10(P1P2/(4*pi)),[-30 25],5,'-b'); hold on;
%polarplot(theta,10*log10(P1P2/(4*pi)));
%Polar_dB(theta,10*log10(P1/(4*pi)),[-100 20],5,'-r'); hold on;
%Polar_dB(theta,10*log10(P2/(4*pi)),[-100 20],5,'-k'); hold on;
%polardb(theta',10*log10((P1P2/(4*pi))),-100,'-k'); hold on;
%polardb(theta',10*log10((P1/(4*pi))),-100,'-k'); hold on;
%polardb(theta',10*log10((P2/(4*pi))),-100,'-k'); hold on;
title(['Unpolarized Incident Light with Particle Radius = ' num2str(a) ' um '])


figure(2);plot(theta*180/pi,P);
%title('Degree of Polarization vs. Angle')
xlabel('Scattering Angle [Degrees]')
ylabel('Degree of Polarization')
axis([0 AngLim -1 1])

figure(3); plot(1:Nc,real(Dmx))
title('Plot of Logrithemic Derivative')

figure(4); semilogy(theta*180/pi,HGPF1,theta*180/pi,HGPF2)
title('Henyey-Greenstein Phase Function')
xlabel('Angle (\theta) (Degrees)')
ylabel('Relative Intensity')
legend('Actual','Estimated')

% figure(5); semilogy(theta*180/pi,IoutUP,theta*180/pi,IoutPep,theta*180/pi,IoutPara);
% title('Scattered Intensity')
% xlabel('Angle (\theta) (Degrees)')
% ylabel('Intensity (Watt)')
% legend('Unpolarized','Perp','Parallel')
time=toc;
%** Display Results **%
disp('        %** Results **%')
disp(['Relative index of refraction (real)= ' num2str(real(m))])
disp(['Relative index of refraction (Imaginary)= ' num2str(imag(m))])
disp(['Size Factor = ' num2str(x)])
disp(['Number of Terms = ' num2str(Nc)])
disp(['Derivative Number of Terms = ' num2str(Nd)]);
disp(['Wave Number (Inverse um) = ' num2str(k)])
disp(['Particle Cross Section (um^2) =' num2str(Area)])
disp(['Particle Volume (um^3) =' num2str(Vol)])
disp(['Scattering Sum = ' num2str(Sscat)])
disp(['Extinction Sum = ' num2str(Sext)])
disp(['Scattering Cross Section (um^2)= ' num2str(Cscat)])
disp(['Extinction Cross Section (um^2)= ' num2str(Cext)])
%disp(['Absorption Efficiency = ' num2str(Qabs)])
disp(['Scattering Efficiency = ' num2str(Qscat)])
disp(['Extinction Efficiency = ' num2str(Qext)])
disp(['Extinction Coefficient 1 (1/m) = ' num2str(ExtCoeR1)])
disp(['Extinction Coefficient 2 (dB/km) = ' num2str((ExtCoeR2*1e4)*(log10(exp(1))))])
%disp(['Extinction Coefficient 3 (1/km) = ' num2str(ExtCoeX)])
disp(['Asymmetry Parameter (Actual) = ' num2str(MeanCOS)])
disp(['Asymmetry Parameter (Estimated) = ' num2str(g2)])

disp(['Run Time = ' num2str(time)])