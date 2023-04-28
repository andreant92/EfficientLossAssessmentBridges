function [FRAGpd, alpha, beta, PowerLawParameters]=cloudFRAGls(IM,EDP,DS,plotter)

% FRAGILITY CURVE ATTAINMENT with least square regression (Jalayer 2017)

%% Fitting Logarithmic Linear Model
lnIMtemp=log(IM);
lnEDPtemp=log(EDP);

[lnIM, lnEDP] = prepareCurveData( lnIMtemp, lnEDPtemp );

ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
%opts.Robust = 'Bisquare';
% opts.Lower = [1 -Inf];
[fitresult, gof] = fit( lnIM, lnEDP, ft, opts );

% 1st order linear regression: find ln(EDP)=b*ln(IM)+ln(a)
b = fitresult.p1; q = fitresult.p2; a=exp(q);
PowerLawParameters = [b a];

% 
% % Unlock to use polyfit with power law
% [xData, yData] = prepareCurveData( IM, EDP );
% % 
% % Set up fittype and options.
% ft = fittype( 'power1' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% % opts.Lower = [-Inf 1];
% % opts.Robust = 'Bisquare';
% % opts.StartPoint = [0.264793642240269 1.02473925712347];
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% a = fitresult.a; b = fitresult.b;
% PowerLawParameters = [b a];

% parameters of the probabilistic seismic demand model: EDP = aIM^b; ln(EDP) = ln(a) + b*ln(IM)
% figure
% hold on
% scatter(log(IM), log(EDP))
% plot(log(IM),  log(a) + b*log(IM),'g')
% plot(log(IM),  log(a1) + b1*log(IM))


%% PLOT
if strcmp(plotter,'plot')
    figure
    
    %logarithmic plane
    subplot(1,2,1)
    hold on
    scatter(log(IM), log(EDP),5,'filled','Displayname','DATA')
    plot(log(IM),b*log(IM)+log(a),'Displayname',['ln(EDP)=',num2str(b,2),'ln(IM)+ln(',num2str(a,2),')'])
    legend
    xlabel('IM'); ylabel('DCR')
    
    %power law
    subplot(1,2,2)
    hold on
    scatter(IM, EDP,'filled','Displayname','DATA')
    [IMord,ord]=sort(IM);
    plot(IMord, a*IMord.^b,'Displayname',['EDP=',num2str(a,2),'IM^{',num2str(b,2),'}'])
    legend
    xlabel('IM'); ylabel('DCR')
end

%% Fragility curve (standardized cumulative distribution)

%Initialize
alpha=zeros(1,length(DS)); beta=zeros(1,length(DS));
FRAGpd=cell(1,length(DS)); 

%Calculate lognormal Gaussian distribution
for ds=1:length(DS)
    alpha(ds)=exp((log(DS(ds)/a))/b); 
    
    IMrev=exp(lnIM); EDPrev=exp(lnEDP);
    res=(log(EDPrev)-log(a*IMrev.^b));
    beta(ds)=((sum(res.^2)/(length(IMrev)-2))^0.5)/b;
    
    FRAGpd{ds}=makedist('lognormal','mu',log(alpha(ds)),'sigma',beta(ds));
%     frag{ds}=[x; logncdf(x, log(alpha(ds)), beta(ds))];
end

% %% NEW method
% for ds=1:length(DS)
%     % calculate the error of the log(data) to the line
%     resid     = log(EDP) - (log(a) + b.*log(IM));
%     
%     FragMean   = log( DS ./ a ) ./ b;
%     beta(ds)  = std( resid ) / b;
%     alpha = exp(FragMean); % intersection of the power law with the EDP threshold line
%     FRAGpd{ds}=makedist('lognormal','mu',log(alpha(ds)),'sigma',beta(ds));
% end

%% PLOT
if strcmp(plotter,'plot')
        x=0:0.01:5;
colour=[.3 .6];

    figure
    hold on
    for ds=1:length(DS)
%         plot(x,frag{ds}(2,:),'color',ones(1,3)*colour(ds),'Linewidth',2,'Displayname',['DS=',num2str(DS(ds),2)])
        plot(x, cdf(FRAGpd{ds},x),'color',ones(1,3)*colour(ds),'Linewidth',2,'Displayname',['DS=',num2str(DS(ds),2)])
%         plot(x, normcdf(log(x/alpha(ds))/beta(ds)),'color',ones(1,3)*colour(ds),'Linewidth',3,'Displayname',['DS=',num2str(DS(ds),2)])
        scatter(alpha(ds),0.5,50,ones(1,3)*colour(ds),'filled','Displayname',['IM=',num2str(alpha(ds),2)])
        legend
    end
    axis([min(IM) max(IM) 0 1])
end

end