function [deltaSDOFbr, VbaseSDOFbr, massSDOFbr, kSDOFbr, tSDOFbr, ethaSDOFbr, csiSDOFbr, edpDSbr, CrMembBr, displSub] = DBAssL(FDpiers, FDbears, Mdeck, Mpiers, DSthresh, varargin )

%% Optional input

% maximum number of optional inputs
numvarargs = length(varargin);
if numvarargs > 5
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 5 optional inputs');
end

% set defaults for optional inputs
optargs = {'noplot', 0, 1, 0.001, 0.5};

% now put these defaults into the valuesToUse cell array, 
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
[ plotter,backfint, deltaMAXcontrol, passoANALISI, powerETHA ] = optargs{:};

%% FD for each subassembly

[deltaSDOF, VbaseSDOF, csiSDOF, ~, massSDOF, EDPds, crMemb, displSUBASS] = DBAssSeries(FDpiers, FDbears, Mdeck, Mpiers, DSthresh, backfint, 'noplot');

for pr=1:length(deltaSDOF)
    FD{pr}=[deltaSDOF{pr}',VbaseSDOF{pr}'];
    FD{pr}(1,:)=[0 0];
end

%% Inizio ciclo for x V punti della curva (senza controllo spostamento ultimo)        
if strcmp(plotter,'plot')
    f = figure(666);
    figure(f)
    hold on
end
% possible values for the displacement of the control node
dispCONTROL = 0.001 : passoANALISI : deltaMAXcontrol;

%% Pre-allocation of results

deltaSDOFbr  = zeros(numel(dispCONTROL), 1);
massSDOFbr   = zeros(numel(dispCONTROL), 1);
keff         = zeros(numel(dispCONTROL), 1);
csiSDOFbr    = zeros(numel(dispCONTROL), 1);
teff         = zeros(numel(dispCONTROL), 1);
ethaSDOFbr   = zeros(numel(dispCONTROL), 1);
VbaseSDOFbr  = zeros(numel(dispCONTROL), 1);

%% calculate curve

for j= 2:length(dispCONTROL)
    
    %% Calcolo vettore dei tagli alla base
    %interpolazione nelle leggi costitutive    
    deformata=dispCONTROL(j);

    Vi= zeros(1,length(FD));   
    for pr=1:(size(FD, 2))
        temp=linspace(1,1+10^-4,length(FD{pr}(:,2))); %Avoid unique points
        Vi(pr) = interp1(FD{pr}(:,1),FD{pr}(:,2).*temp',deformata,'linear','extrap');
        massSUB(pr) = interp1(FD{pr}(:,1),massSDOF{pr}'.*temp',deformata,'linear','extrap');
        csi(pr) = interp1(FD{pr}(:,1),csiSDOF{pr}'.*temp',deformata,'linear','extrap');
    end 
    
    % Calcolo taglio totale 
    VbaseSDOFbr(j)= sum(Vi);
    
    % Definizione DELTA eqSDOF
    deltaSDOFbr(j) = sum(massSUB.*(deformata.^2))/sum(massSUB.*deformata);
    
    % Definizione DELTA eqSDOF
    massSDOFbr(j) = sum(massSUB.*deformata)/deltaSDOFbr(j);       
    kSDOFbr(j) =  VbaseSDOFbr(j)/deltaSDOFbr(j);
    tSDOFbr(j) = 2*pi*sqrt(massSDOFbr(j)/kSDOFbr(j));  
    
    % Calcolo smorzamento       
    csiSDOFbr(j)   = sum(csi .* Vi .* deformata)./sum(Vi .* deformata);

    % reduction factor for the demand
    ethaSDOFbr(j) = (0.07/(0.02+csiSDOFbr(j)))^powerETHA;                     % eta riduzione spettro
    
    
     %% Checking plot
    % Show the deformed shape

    if strcmp(plotter,'plot')
        f=figure(666);
        hold on
        scatter(deformata,VbaseSDOFbr(j),20,'k')
        axis([0 deformata+0.2 0 max(VbaseSDOFbr)])
    end
    
end

%% arrange subassembly displacement profile
% interp subassembly FDlaw to obtain displacement profile of the
% subassembly step by step

for pr=1:length(deltaSDOF)   
    displSub{pr}(1,:)=dispCONTROL;
    displSub{pr}(2,:)=interp1(displSUBASS{pr}(1,:),displSUBASS{pr}(2,:),dispCONTROL,'linear','extrap');    
end

%% EDP calculation
%find damage state threshold of the bridge
for ds=1:size(EDPds,2)
    [edpDSbr(ds), crsub(ds)]=min(EDPds(:,ds));
    CrMembBr{1,ds}=crMemb(crsub(ds));
    CrMembBr{2,ds}=crsub(ds);
end

end
