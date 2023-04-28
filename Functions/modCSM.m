function [ PPgood, PP, numberSolutions, INspectrum ] = modCSM( ADRSspectrum, pushover, deltaY, plotter, varargin)
%version6b
%Capspectrum is a 2-col vector [#acc [m/s^2] #displ [m]]
% ADRSspectrum is a 3-col vector [#periods [s] #acc [g] #displ [m]]

%EVD is calculated in 3 ways according to strategy in varargin
% 1 Kowalsky 2002, varargin={deltay, 0.5, 1, 0.444}
% 2 Grant 2005 for Takeda (requires Teff) varargin={0, 0.5, 2, 1, 1, #vectorTeff}
% 3 Sadan 2012 (requires wighted average of EVD calculated with Kowalsky or
% Grant for each member step by step)  varargin={0, 0.5, 3, 1, 1, #vectorEVDsdof}
% vectorEVDsdof in [%]


%The performance of New Building Standerd is calculated and the CSM is
%performed even if the NBS<1

%% Example

% CASE 1: NBS>100%

% ADRSspectrum= [0,0.356939000000000;0.000188707833496794,0.474636000000000;0.000197305354375608,0.476990000000000;0.000207250202706988,0.479638000000000;0.000217487361514447,0.482287000000000;0.000228019146043329,0.484935000000000;0.000238849186280760,0.487583000000000;0.000252495891044776,0.490820000000000;0.000265297306034601,0.493762000000000;0.000279819960536491,0.496999000000000;0.000296195889566287,0.500530000000000;0.000313135550388310,0.504061000000000;0.000332131634261809,0.507886000000000;0.000351807838722964,0.511711000000000;0.000375369253572282,0.516125000000000;0.000399863067162004,0.520538000000000;0.000425305581112906,0.524952000000000;0.000455304662100784,0.529954000000000;0.000488439968950824,0.535250000000000;0.000524975350849041,0.540841000000000;0.000565195151404757,0.546726000000000;0.000611564580633147,0.553199000000000;0.000635599985137318,0.556436000000000;0.000662475598276850,0.559967000000000;0.000690042553940875,0.563498000000000;0.000718308433028745,0.567029000000000;0.000749726981842801,0.570854000000000;0.000784500343125123,0.574973000000000;0.000820259915965028,0.579093000000000;0.000857014957194100,0.583212000000000;0.000897516412943632,0.587626000000000;0.000942009263002863,0.592333000000000;0.000987856232472527,0.597041000000000;0.00103807206952885,0.602044000000000;0.00109294995354858,0.607340000000000;0.00115281244030097,0.612930000000000;0.00121469499481279,0.618521000000000;0.00128204746626547,0.624406000000000;0.00135524808896455,0.630585000000000;0.00143470344543515,0.637058000000000;0.00152466942151269,0.644120000000000;0.00161812347293649,0.651182000000000;0.00169871010507974,0.657067000000000;0.00178179622463277,0.662952000000000;0.00186741692835183,0.668837000000000;0.00200067562723999,0.677664000000000;0.00209278628501490,0.683549000000000;0.00223594429063482,0.692376000000000;0.00233476425222237,0.698261000000000;0.00248815058299698,0.707088000000000;0.00264783180177572,0.715915000000000;0.00281393023606115,0.724743000000000;0.00298655675398161,0.733570000000000;0.00322709360051803,0.745340000000000;0.00341542169501413,0.754167000000000;0.00367732420679258,0.765937000000000;0.00395181672691765,0.777706000000000;0.00423918990394306,0.789476000000000;0.00461693791863132,0.804188000000000;0.00493430819830476,0.815958000000000;0.00535044062235630,0.830670000000000;0.00587901378799521,0.848325000000000;0.00634447338124399,0.863037000000000;0.00695335985993069,0.873332000000000;0.00742731940969782,0.873332000000000;0.00800002053233310,0.873332000000000;0.00868057783456283,0.873332000000000;0.00903127317907916,0.873332000000000;0.00938891298586315,0.873332000000000;0.00984572839440702,0.873332000000000;0.0102189932412932,0.873332000000000;0.0106953399499649,0.873332000000000;0.0111825373809297,0.873332000000000;0.0117814972515145,0.873332000000000;0.0122925662715244,0.873332000000000;0.0129201720489633,0.873332000000000;0.0135634028665044,0.873332000000000;0.0142222587241477,0.873332000000000;0.0150106722059719,0.873332000000000;0.0158203531034908,0.873332000000000;0.0167717444341588,0.873332000000000;0.0177509136138975,0.873332000000000;0.0187578606427068,0.873332000000000;0.0199238792603345,0.873332000000000;0.0211250542181921,0.873332000000000;0.0226409001225526,0.873332000000000;0.0240645148874209,0.873332000000000;0.0258301444189710,0.873332000000000;0.0276582741109299,0.873332000000000;0.0297092776387913,0.873332000000000;0.0321669662382019,0.873332000000000;0.0347223113382513,0.873332000000000;0.0377364249768574,0.873332000000000;0.0410645585186288,0.873332000000000;0.0449274156550092,0.873332000000000;0.0491702650860977,0.873332000000000;0.0542536114660177,0.873332000000000;0.0600426888238876,0.873332000000000;0.0635270351154456,0.826989000000000;0.0671832848857745,0.781983000000000;0.0714108427183093,0.735690000000000;0.0762095747757297,0.689364000000000;0.0815797030755833,0.643986000000000;0.0878638269291779,0.597927000000000;0.0951763202548472,0.551988000000000;0.103859835458552,0.505837000000000;0.114257286226744,0.459806000000000;0.119970020080725,0.437910000000000;0.126825698041336,0.414240000000000;0.134823562561345,0.389666000000000;0.142821685436619,0.367845000000000;0.151962279914402,0.345719000000000;0.163388093421329,0.321543000000000;0.175956029451730,0.298575000000000;0.190809713646313,0.275333000000000;0.207948541329019,0.252641000000000;0.228514572453488,0.229903000000000;0.253651396082672,0.207120000000000;0.285642594341353,0.183922000000000;0.304411017817532,0.149768000000000;0.304411985612325,0.110475000000000;0.304410438139586,0.0765650000000000;0.304400498531471,0.0122500000000000]; 
% pushover    = [0, 0; 12.54/1000, 132.5; 90/1000, 132.5]; % capacity curve [m, kN]
% massa       = 0.8 * 3278.5 / 9.81; % effective mass [Ton]
% deltaY      = pushover(2,1);
% plot        = 'plot';
%
% [NBS, dispDEM, accDEM] = CSMnbs( ADRSspectrum, pushover, massa, deltaY );

%% Optional input

% Maximum number of optional inputs
numvarargs = length(varargin);
if numvarargs > 5
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 5 optional inputs');
end


% set defaults for optional inputs
optargs = { 0.5, 0.0005, 1, 0.444, 1};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
[ alpha, step, strategy, Cevd, EVDsdof] = optargs{:};

%% Initial stuff
% periods = ADRSspectrum(:,1);

%Vector elaboration in order to avoid unique points and prevent
%interpolation errors
pushover(1,:)=[0 0];

[~, indexes, ~] = unique(pushover(:,1),'rows');
pushover=pushover(indexes,:);
[~, indexes, ~] = unique(ADRSspectrum(:,1),'rows');
ADRSspectrum=ADRSspectrum(indexes,:);

% calculate the capacity spectrum

% interpolate pushover every 0.005m to reduce computational time (for the intersections)
% calculate the capacity spectrum, periods and evdsdof (if required)
Cspectrum=zeros(length([0 : step : pushover(end,1)]),2);
Cspectrum(:,1) = [ 0 : step : pushover(end,1) ]; 
Cspectrum(:,2) = [ interp1(pushover(1:end,1), pushover(1:end,2), Cspectrum(1:end,1))]/9.81;
Teff = 2*pi * ( Cspectrum(:,1)./(Cspectrum(:,2)*9.81) ).^0.5;

ADRSnew(:,1)=Teff;
ADRSnew(:,2)= interp1(ADRSspectrum(:,1),ADRSspectrum(:,2),Teff);
ADRSnew(:,3)=interp1(ADRSspectrum(:,1),ADRSspectrum(:,3),Teff);

if strategy==3
      temp = [ interp1(pushover(1:end,1), EVDsdof, Cspectrum(2:end-1,1)); EVDsdof(end)]; 
      EVDsdofnew=temp;
      EVDsdofnew(isnan(EVDsdofnew))=5;
end

clear pushover

%% check if CSM is needed 

% If the structure is not able to widstand the demand compatible with the
% last point of the capacity spectrum (max damping), no further action is required. 
%
% Oterwise, the Capacity Spectrum Method is needed to calculate the
% displacement demand related to the 100% spectrum

CSMdisp = Cspectrum(end,1);

% guess the displacement demand (Last point)
Ddemand = Cspectrum(end,1);

% ductility demand
MUdemand = Ddemand / deltaY;

% equivalent viscous damping
EVD = (0.05 + Cevd*(MUdemand - 1)/(pi * MUdemand)) * 100; % [%]

% reduction factor
etha = ( 7 / (2+EVD) )^alpha;

% reduce spectrum
INspectrum = etha * ADRSspectrum;

% displacement demand for secant-to-ultimate 

[NdispCSM] = intersections(Cspectrum(:,1), Cspectrum(:,2) , INspectrum(:,3), INspectrum(:,2));

% figure
% hold on
% plot(Cspectrum(:,1), Cspectrum(:,2))
% plot(INspectrum(:,3), INspectrum(:,2))

%% Run CSM or Inverse CSM

if ~isnan(NdispCSM) % CSM
    % guess displacement demand
    PPxOLD = Cspectrum(2:end,1);
    
    % acceleration on the capacity spectrum
    PPyOLD = Cspectrum(2:end,2);
    
    % ductility demand based on each trial performance point
    MUdemand = max(1, PPxOLD / deltaY);
    
    %Calculate equivalent viscous damping step-by-step
    %and inelastic demand
    
    % Are there elastic intersection'
    ELsteps=find(Cspectrum(:,1)<deltaY);
    % Elastic field intersection?
    eldem=intersections(Cspectrum(ELsteps,1),Cspectrum(ELsteps,2), ADRSspectrum(:,3),ADRSspectrum(:,2));
%     figure 
%     hold on
%     plot(Cspectrum(ELsteps,1),Cspectrum(ELsteps,2))
%     plot( ADRSspectrum(:,3),ADRSspectrum(:,2))
    if isempty(eldem)
        if strategy==1
            % equivalent viscous damping %Kowalsky 2002
            %----------------------------------------------------------------------
            EVD = (0.05 + Cevd.*(MUdemand - 1)./(pi .* MUdemand)) * 100; % [%]
        elseif strategy==2            
            % equivalent viscous damping on SDOF by Grant 2005
            %----------------------------------------------------------------------
            a=0.215; b=0.642; c=0.824; d=6.444; %Takeda
            % effective damping according to Grant 2005     
            EVD = (0.05*MUdemand.^(-0.378) + a.*(1-1./(MUdemand.^b)).*(1+1./((Teff(2:end)+c).^d))) * 100; % [%]           
        elseif strategy==3       
            EVD=EVDsdofnew;            
        else
            warning('input is wrong!Check strategy and associate varargin')
        end
        
        % reduction factor
        etha = (7 ./ (2+EVD) ).^alpha;
        
        % reduce spectrum
        INspectrum = [ADRSnew(2:end,1), etha .* ADRSnew(2:end,2:3)];
               
        % intersect demand and capacity
        [Tpp, PPx] = intersections(INspectrum(:,1), INspectrum(:,3), Teff(2:end), Cspectrum(2:end,1));
        
        [~, indexes, ~] = uniquetol(PPx(:));
        PPx=PPx(indexes); Tpp=Tpp(indexes);
        numberSolutions = length(PPx);
        
        PPy=interp1(Cspectrum(:,1), Cspectrum(:,2), PPx);
        PP=[PPx,PPy];
        PPgood= [PPx(1),PPy(1)];

    else
        % reduce spectrum
        INspectrum = ADRSnew;
        PPx=eldem;
        numberSolutions = length(PPx);
        PPy=interp1(Cspectrum(2:end-1,1), Cspectrum(2:end-1,2), PPx);
        
        PP=[PPx,PPy];
        PPgood= [PPx(1),PPy(1)];
    end

else
    
    PPgood=inf;
    PP=inf;
    numberSolutions = 0;
    
end

%% select the "good" solution
%if multiple solution, find the exact PP as the nearest to avgSd

if numberSolutions>1
    Tel=Teff(2); Trange=Tel:0.01:1.5*Tel;
    avgSd=geomean(interp1(ADRSspectrum(:,1), ADRSspectrum(:,3), Trange));
    indPP=knnsearch(PPx,avgSd);
    PPgood=[PPx(indPP), PPy(indPP)];
end


% figure
% hold on
% plot(Cspectrum(:,1),Cspectrum(:,2))
% plot(ADRSspectrum(:,1),ADRSspectrum(:,2))
% plot(ADRSnew(:,1),ADRSnew(:,2))
% plot(INspectrum(:,1), INspectrum(:,2))
% scatter(PPx,PPy)
% plot([0, Cspectrum(2,1)*1000], [0, Cspectrum(2,2)*1000])
% 
% figure
% hold on
% plot(Teff,Cspectrum(:,1))
% plot(ADRSspectrum(:,1),ADRSspectrum(:,3))
% plot(Teff,ADRSnew(:,3))
% plot(Teff(2:end), INspectrum(:,3))
% scatter(Tpp,PPx)


%% Plot

if strcmpi(plotter, 'plot')
    figure 
    hold on
    plot(ADRSspectrum(:,3), ADRSspectrum(:,2), 'k-.')
    plot(INspectrum(:,3), INspectrum(:,2), 'k')
    plot(Cspectrum(:,1), Cspectrum(:,2), 'k','linewidth',1.5)
    scatter(PPx,PPy,'filled','k')
    scatter(PPgood(1),PPgood(2),'filled','r')
    legend ('ADRSspectrum','INspectrum','Cspectrum','PP')
    xlim([0 PPx(end)*2])
end

end


