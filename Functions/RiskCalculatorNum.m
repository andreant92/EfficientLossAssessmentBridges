function lambda=RiskCalculatorNum(HazardCurve, FragilityCurve, varargin)
% Calculation of risk as MAF (Mean Annual Frequency) of exceeding a LS (Eads 2013)

% Optional input

% Maximum number of optional inputs
numvarargs = length(varargin);
if numvarargs > 3
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 3 optional inputs');
end

% set defaults for optional inputs
optargs = {'nolin','noplot'};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[ linearize, plotter ] = optargs{:};

%% Numerical calculation of MAF (Eads 2013)
%1- define dIM (differential operator) and prepare data
dIM=0.0001;
xH=HazardCurve(1,1):dIM:HazardCurve(end,1); %Abscissas for hazard curve [m/s^2 or g]

xF=xH(1:end-1)+dIM/2; %Abscissas for fragility curve [m/s^2 or g]
% for ds=1:size(FragilityParam,2)
%     FragilityCurve(ds,:)=normcdf(log(xF./FragilityParam(1,ds))./FragilityParam(2,ds));
% end

for ds=1:size(FragilityCurve,2)
    FragilityCurveResampl(ds,:)=interp1(FragilityCurve{ds}(:,1),FragilityCurve{ds}(:,2),xF,'linear','extrap');
end

%2- Interp Hazard and Fragility curves
%Find Hazard function (poly4 in log log plane) (Eads 2013)

if strcmp(linearize,'lin')
    [fitresult] = FitHazardCurve((HazardCurve(:,1)), (HazardCurve(:,2)), 'noplot');
    %Hfunct=exp(fitresult(log(xH)));
    Hfunct=fitresult(xH)';
else
    Hfunct=interp1(HazardCurve(:,1),HazardCurve(:,2),xH,'spline','extrap');
end


%%% Control Plot
% figure
% subplot(2,1,1)
% hold on
% scatter(HazardCurve(:,1),HazardCurve(:,2))
% plot(xH,Hfunct)
% 
% subplot(2,1,2)
% hold on
% scatter(FragilityCurve{ds}(:,1),FragilityCurve{ds}(:,2))
% plot(xF,FragilityCurveResampl)


for ds=1:size(FragilityCurve,2)
    Ffunct=FragilityCurveResampl(ds,:);
    
    %3- %Find derivative of Hazard Curve
    derivH=abs(diff(Hfunct)/dIM);
    
    %4- %Find Deaggregation for lambda (MAF) (Mean Annual Frequency of
    %exceeding a limit state given a value of IM)
    DeagLambda=derivH.*Ffunct;
    
    %5- %Integration for risk calculation lambda (MAF)
    lambda(ds)=sum(DeagLambda*dIM);
end

% PLOT?
if strcmpi(plotter,'plot')
    Xlim=max(Hfunct);
    
    figure
    subplot(3,1,1) %Fragility curve
    hold on
    plot(xF,Ffunct,'k','Linewidth',2)
    scatter(xF,Ffunct,1,[0 0 0],'o','filled')
    xlim([0 Xlim]) 
    ylim([0 1])
    
    subplot(3,1,2) %Hazard Curve and Derivative
    hold on
    plot(xH,Hfunct,'k','Linewidth',2)
    plot(xF,derivH,'r','Linewidth',2)
    scatter(xF,derivH,1,[0 0 0],'o','filled')
    
    subplot(3,1,3) %deaggregation of lambda given IM
    hold on    
    plot(xF,DeagLambda,'k','Linewidth',2)
    scatter(xF,DeagLambda,20,[0 0 0],'o','filled')
    
end

end

function [fitresult, gof] = FitHazardCurve(PGAhaz, Freq, plotter)

%% Fitting with LS
% PGAhaz=exp(logPGAHaz);
% Freq=exp(logFreq)
[xData, yData] = prepareCurveData( PGAhaz, Freq );

% Set up fittype and options.
ft = fittype( 'exp1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

if strcmp(plotter,'plot')
% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'logXFreq vs. logPGAHaz', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'logPGAHaz', 'Interpreter', 'none' );
ylabel( 'logXFreq', 'Interpreter', 'none' );
grid on
end
end


