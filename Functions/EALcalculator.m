function [EAL, lossExceedanceCurve]=EALcalculator(HazardCurve,LossCurve,varargin);
%This function calculate Expected Annual Losses by means of numerical
%integration between an hazard curve and a vulnerability (loss exceedance) curve

% Optional input

% Maximum number of optional inputs
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 1 optional inputs');
end

% set defaults for optional inputs
optargs = {'noplot'};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[ plotter ] = optargs{:};

%% RUN

%1- define dIM (differential operator) and prepare data
dIM=0.01;
xH=HazardCurve(1,1):dIM:HazardCurve(end,1); %Abscissas for hazard curve [m/s^2 or g]

xL=xH(1:end-1)+dIM/2; %Abscissas for fragility curve [m/s^2 or g]

%Resampling Loss Curve
LossCurveResampl=interp1(LossCurve(:,1),LossCurve(:,2),xL,'linear','extrap');

%Resampling Hazard Curve
HCurveResampl=interp1(HazardCurve(:,1),HazardCurve(:,2),xH,'spline','extrap');

%%- Calculate EAL via numerical integration

%2- %Find derivative of Hazard Curve
derivH=abs(diff(HCurveResampl)/dIM);

%4- %Find Deaggregation for lambda (MAF) (Mean Annual Frequency of
%exceeding a limit state given a value of IM)
tempProd=derivH.*LossCurveResampl;

%5- %Integration for risk calculation lambda (MAF)
EAL=sum(tempProd*dIM);

%6 - Calculate MAFE | TRC
LossMAFE=interp1(LossCurve(:,1),LossCurve(:,2),...
    xH,'linear','extrap');
lossExceedanceCurve=[LossMAFE',HCurveResampl'];

%% Final plot
if strcmpi(plotter,'plot')
    Xlim=max(xH);
    
    figure
    subplot(4,1,1)
    hold on
    plot(xH,HCurveResampl,'k','Linewidth',2)
    xlim([0 Xlim])

    subplot(4,1,2)
    hold on
    plot(xL,LossCurveResampl,'k','Linewidth',2)
    xlim([0 Xlim])

    subplot(4,1,3)
    loglog(NaN,NaN,'Parent',gca)
    hold on
    loglog(lossExceedanceCurve(:,1),lossExceedanceCurve(:,2),'k','Linewidth',2)

    subplot(4,1,4)
    %loglog(NaN,NaN,'Parent',gca)
    hold on
    plot(xL,tempProd,'k','Linewidth',2)
   
end

end
