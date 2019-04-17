function [R,WeightsNX,WeightsNR,bx]=SizeParamGammaDistribution(N,ALPHA,g,Rc,MaxR,MinR,Resolution,K)
    R=MinR:1/Resolution:MaxR;
    X=K*R; Xc=K*Rc;
    bx=ALPHA./(g*(Xc/K).^g);
    ax=g*K*N./(bx.^(-(ALPHA+1))*gamma((ALPHA+1)/g));
   
    b=ALPHA./Rc;
    a=N./(b.^(-ALPHA-1)*gamma(ALPHA+1));
    
    %** Allocation **%
    WeightsNX=zeros(length(Rc),length(R));
    WeightsNR=zeros(length(Rc),length(R));
    
    for mm=1:length(Rc)
% %         %** Cloud C model **%
%         WeightsNR(mm,:)=2.373*(R.^6).*exp(-1.5*R);
%         WeightsNX(mm,:)=2.373*K*(X/K).^6.*exp(-1.5*(X/K));
        %** Haze M Model **%
%         WeightsNR(mm,:) = 5.33e4.*R.*exp(-8.944*sqrt(R));
%         WeightsNX(mm,:) = 5.33e4*K*(X/K).*exp(-8.944*sqrt(X/K));
        %** Standard Modified Gamma Distributuion **%
        WeightsNX(mm,:)=ax(mm)*((X/K).^ALPHA).*exp(-bx(mm)*(X/K).^g);
        WeightsNR(mm,:)=a(mm)*R.^ALPHA.*exp(-b(mm)*R);
    end
end