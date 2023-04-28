function meanLoss=VulnCurveCalculator(DLRMatrix, FragilityCurve, IMvect, varargin)
% Calculation of Vulnerability Curves ( Expected Loss given IM)

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

%% Calculate probability for given DS

% tent1=[ones(1,length(FragilityCurve))-FragilityCurve(1,:);
% FragilityCurve(1,:)-FragilityCurve(2,:);
% FragilityCurve(2,:)-FragilityCurve(3,:);
% FragilityCurve(3,:)-zeros(1,length(FragilityCurve))];

Pds=[ones(1,length(FragilityCurve)); FragilityCurve]...
    -[FragilityCurve; zeros(1,length(FragilityCurve))];

meanLoss(:,1)=IMvect;
meanLoss(:,2)=sum(Pds.*DLRMatrix');


%% PLOT?
if strcmpi(plotter,'plot')
    
    colors={'r','g','b','y','c'};

    figure
    subplot(2,1,1) %Fragility curve
    hold on
    for ds=1:size(FragilityCurve,1)
    plot(IMvect,FragilityCurve(ds,:),'color',colors{ds},'Linewidth',2)
    xlim([0 2])
    end
    
    subplot(2,1,2) %Hazard Curve and Derivative
    hold on
    plot(meanLoss(:,1),meanLoss(:,2),'k','Linewidth',2);

    for ds=1:size(FragilityCurve,1)
    plot(IMvect,Pds(ds+1,:)*DLRMatrix(ds+1),'color',colors{ds},'Linewidth',2)
    end
    

end

end



