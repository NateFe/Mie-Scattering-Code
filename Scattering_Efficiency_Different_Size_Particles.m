clear;close all; clc;
%** Mie Theory Parameters **%
AngLim=180;             %Display Angle Limit in Degrees
N0=1;               %Index of medium
N1=1.335-0i;       %Ind0ex of Particles
%lambda=.850;            %wavelegth in um
%lambda=.45:.2:.85;
lambda = [0.45 0.65 0.85];
Density=1;

k=(2*pi*N0)./(lambda);   %wave number in medium
m=N1/N0;                %relative index of refraction
Tres=pi/180;              %theta resolution

Radii=0.01:0.01:10; QextMat=0; Qs=0;

for nn=1:length(lambda)
    for n=1:length(Radii)
        a=Radii(n);
        [Sscat,Sext,Cscat,Cext,Qscat,Qext,P1P2,P1,P2,theta,P,ExtCoe] = MieFunction(a,m,lambda(nn),k(nn),AngLim,Tres,Density);
        QextMat(nn,n)=Qext;
        if k*Radii(n) < 1
            Qs(nn,n) = (32/27)*(m-1)^2*(k(nn)*Radii(n))^4; % Rayleight Gans Approximation for x << 1
        end
    end
end

plot(Radii,QextMat); hold on;
plot(Radii(1:length(Qs)),Qs);
%title('Extinction Efficiency vs. Particle Radius')
xlabel('Particle Radius [um]')
ylabel('Extinction Efficiency []')
legend('450 \mum','650 \mum','850 \mum')
%legend('0.45','0.5','0.55','0.65','0.7','0.75','0.8','0.85','0.9','0.95')