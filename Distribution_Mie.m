clear;close all;clc


%** Mie Theory Parameters **%
AngLim=360;             %Display Angle Limit in Degrees
N0=1;               %Index of medium
N1=1.335;1.315-0.0143i;       %Ind0ex of Particles
lambda=0.850;            %wavelegth in um

k=(2*pi*N0)/(lambda);   %wave number in medium
m=N1/N0;                %relative index of refraction
Tres=0.001;180/pi;              %theta resolution
%*************************************************************************%
% 
% % ** Guassian Parameters **%
% u=5; %Mean Radius 
% CV=10; %Constant of Variation Percent
% sigma=(CV/100)*u;
% NumofSTDEV=4;
% Res=100;
% Density=100;
% 
% [Radii,WeightR,WeightX,SizeFactPDF]=GuassianDistribution(u,sigma,Res,NumofSTDEV,k);
% WeightR=WeightR;
% WeightX=SizeFactPDF/Res;
% %figure(1);plot(SizeFactPDF,WeightX);
% figure(2);plot(Radii,WeightR);
% xlabel('Particle Radii [um]') 
% ylabel('Distribution Weight')

%** Gamma Distribution Parameters **%
%** (1/(gamma(k)*theta^k))*D.^(k-1).*exp(-D/theta);
Density=100;
Rc=8;
ALPHA=6;
Phi=1;
MaxR=25;       %in um
MinR=0.1;        %in um
Res=100;

%[Radii,Weight,NX]=RadiusGammaDistribution(Density,ALPHA,Rc,MaxR,MinR,Res,k);
[Radii,WeightX1,WeightR1]=SizeParamGammaDistribution(Density,ALPHA,Phi,Rc,MaxR,MinR,Res,k);
WeightX=WeightX1/Res;
WeightR=WeightR1/Res;
figure(1); plot(Radii*k,WeightX);
xlabel('Particle Size Factor []')
ylabel('Distribution Weight')
%title(['Gamma distribution where R_c = ' num2str(Rc) ' Normalized to ' num2str(Density) ' [1/cm^3]' ])
figure(2); plot(Radii,WeightR);
xlabel('Particle Radii [um]')
ylabel('Distribution Weight')
% %title(['Gamma distribution where R_c = ' num2str(Rc) ' Normalized to ' num2str(Density) ' [1/cm^3]' ])
% %** Log-Normal (Junge) Distribution **%
% u=.5; %Particle RWadius 
% sigma=.25;
% NumofSTDEV=4;
% Resolution=100;
% Len=10;
% Start=1/Resolution;
% 
% [Radii,Weight,NX]=LogNormalDistribution(u,sigma,Resolution,NumofSTDEV,Len,Start,k);
% figure(1);plot(Radii,Weight)

tic;
%** Declarations **%
SscatTOT=0; SextTOT=0; CscatTOT=0; CextTOT=0; QscatTOT=0; QextTOT=0; 
P1P2TOT=0; P1TOT=0; P2TOT=0; DegPol=0; Qextmat=0; QabsTOT=0;
HGPF1TOT=0; HGPF2TOT=0; g=0; IoutParaTOT=0; IoutPepTOT=0; IoutUPTOT=0;
SextMat=0; gMat=0; gTOT=0;
for n=1:length(Radii)
    a=Radii(n); %Particle Radius
    %a=SizeFactPDF(n)/k;
    x=k*a;                  %Size Factor
    
    [Sscat,Sext,Cscat,Cext,Qscat,Qext,P1P2,P1,P2,theta,P,g] = MieFunction(a,m,lambda,k,AngLim,Tres,Density);
    %Useful values toQs know: Divide by Res to match the single particle case in a small sigma limit   
   
    QextMat(n)=Qext;
    Area(n)=pi*a^2;
    CextMat(n)=Cext;
    SizeFactor(n)=x;
    QscatMat(n)=Qscat;
    CscatMat(n)=Cscat;
    SextMat(n)=Sext;
    gMat(n)=g;
    P1P2Mat(n,:)=P1P2;
    P1Mat(n,:)=P1;
    P2Mat(n,:)=P2;
%     %** Scaled Parameters by Weight in Distribution
    CscatTOT=CscatTOT+Cscat*WeightR(n);
    CextTOT=CextTOT+Cext*WeightR(n);
    QscatTOT=QscatTOT+Qscat*WeightR(n);
    QextTOT=QextTOT+Qext*WeightR(n);
    QabsTOT=Qext*WeightR(n)-Qscat*WeightR(n);
    gTOT=gTOT+g*WeightR(n);
    DegPol = DegPol +P*WeightR(n);

    P1P2TOT=P1P2TOT+P1P2*WeightR(n);
    P1TOT=P1TOT+P1*WeightR(n);
    P2TOT=P2TOT+P2*WeightR(n); 
end
%** NOTE FOR CHANGING VARIABLE **%
% If we are working with n(x) change the weights on the matrix elements
% P1... to be P1*WeightX(n) and change the Normalization to
% (4*pi)/(k^3*Sigma3) instead since all must be relative to n(x)
Sigma1=sum(Area.*QextMat.*WeightR);
SigmaT=trapz(Area.*QextMat.*WeightR)*1e-8*1e2*1e3;
Sigma2=sum(CextTOT);
Sigma3=(pi/k^3)*sum(QextMat.*WeightX.*SizeFactor.^2);

%** Finding Angular Scattering functions normalization **%
Normalize=(4*pi)/(k^2*Sigma1); %(4*pi)/(k^3*Sigma3)
P1P2TOT=Normalize*(P1P2TOT);
P1TOT=Normalize*P1TOT;
P2TOT=Normalize*P2TOT;

%** Calculating Sigma in terms of 1/m or 1/km
    %** If Density of 100 cm^-3
    Sigma1=Sigma1*1e-8*1e2*1e3;
    Sigma2=Sigma2*1e-8*1e2*1e3;
    Sigma3=Sigma3*1e-8*1e2*1e3;
%     %** if Density is 100 m^-3 = 0.0001 cm^-3
%     Sigma1=sum(Area.*QextMat.*WeightR)*1e-12;
%     Sigma2=sum(CextTOT)*1e-12;
%     Sigma3=(pi/k^3)*sum(QextMat.*WeightX.*SizeFactor.^2)*1e-12;

%** Plotting **%
figure(3);
semilogy(theta*180/pi,P1P2TOT'/(4*pi)); hold on;
semilogy(theta*180/pi,P1TOT'/(4*pi));
semilogy(theta*180/pi,P2TOT'/(4*pi));
%title('Normalized phase function')
xlabel('Angle [Degrees]')
ylabel('Relative Intensity []')
legend('Unpolarized','Perp','Parallel')%,'B(\theta)')

figure(4);plot(theta*180/pi,DegPol/Density);
%title('Degree of Polarization vs. Angle')
xlabel('Scattering Angle [Degrees]')
ylabel('Degree of Polarization [%]')
%axis([0 AngLim -20 80])
time=toc;

% 
% figure(5); semilogy(theta*180/pi,IoutUPTOT,theta*180/pi,IoutPepTOT,theta*180/pi,IoutParaTOT);
% title('Scattered Intensity')
% xlabel('Angle (\theta) (Degrees)')
% ylabel('Intensity (Watt)')
% legend('Unpolarized','Perp','Parallel')
%%
figure(5); 
Polar_dB(theta,10*log10(P1P2TOT/(4*pi)),[-30 30],5,'-b'); hold on;
%Polar_dB(theta,10*log10(P1/(4*pi)),[-100 20],5,'-r'); hold on;
%Polar_dB(theta,10*log10(P2/(4*pi)),[-100 20],5,'-k'); hold on;
%polardb(theta',10*log10((P1P2/(4*pi))),-100,'-k'); hold on;
%polardb(theta',10*log10((P1/(4*pi))),-100,'-k'); hold on;
%polardb(theta',10*log10((P2/(4*pi))),-100,'-k'); hold on;
title(['Polar Scattering Distribution Plot with Unpolarized Incident Light '])

disp('        %** Results **%')
disp(['Relative index of refraction (real)= ' num2str(real(m))])
disp(['Relative index of refraction (Imaginary)= ' num2str(imag(m))])
disp(['Size Factor = ' num2str(x)])
%disp(['Number of Terms = ' num2str(Nc)])
%disp(['Number of Terms for Dn =' num2str(x)])
disp(['Wave Number (Inverse um) = ' num2str(k)])
%disp(['Particle Cross Section (um^2) =' num2str(Area)])
%disp(['Particle Volume (um^3) =' num2str(Vol)])
disp(['Scattering Sum = ' num2str(SscatTOT)])
disp(['Extinction Sum = ' num2str(SextTOT)])
disp(['Scattering Cross Section (um^2)= ' num2str(CscatTOT)])
disp(['Extinction Cross Section (um^2)= ' num2str(CextTOT)])
disp(['Absorption Efficiency = ' num2str(QabsTOT)])
disp(['Scattering Efficiency = ' num2str(QscatTOT)])
disp(['Extinction Efficiency = ' num2str(QextTOT)])
disp(['Extinction Coefficient 1 (1/km)= ' num2str(Sigma1)])
disp(['Extinction Coefficient 2 (1/km)= ' num2str(Sigma2)])
disp(['Extinction Coefficient 3 (1/km)= ' num2str(Sigma3)])
disp(['Visibility 1&2 (km)= ' num2str(-log(0.05)/Sigma1)])
disp(['Visibility 3 (km)= ' num2str(-log(0.05)/Sigma3)])
disp(['Asymmetry Parameter = ' num2str(gTOT)])
%disp(['Extinction Bex = ' num2str(Bext)])
%disp(['Volume Scattering Coefficient = ' num2str()])
%disp(['Paper Extinction Coefficient (1/m) = ' num2str(SigmaTOT)])
disp(['Run Time = ' num2str(time)])