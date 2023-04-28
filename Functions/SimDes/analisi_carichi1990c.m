function [comb,soll]=analisi_carichi1990c(ldeck,l_marciapiede,l_sedestrad,h_girder,h_bent,...
    h_slab,Pp_bent,Pp_bentcap,tribl,tribl_sism,Ppdeck,Lx_pier,Ly_pier,categoria_ponte,Ss,sigma_r,d_bear,piertype,varargin)

%Varargin for seismic protection coefficient (seismic code 1986)
% Maximum number of optional inputs
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 1 optional inputs');
end

% set defaults for optional inputs
optargs = {1};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
 Iseismprot = optargs{:};

%% input elaboration

E_pila=(18000*(sigma_r*10)^0.5)*10^(-1); %kN/mq modulo di elasticità longitudinale del cls della pila

%bridge structural typology definition
if d_bear==0
    bridgetype='cont';
elseif d_bear~=0
    bridgetype='simplysupp';
else
    warning('check number of bearings')
end

if  strcmp(bridgetype, 'cont') 
    if tribl_sism==0
        beartype{1}='free';     
    else
        beartype{1}='fixed';
        tribl_sism=sum(tribl_sism);
    end
elseif strcmp(bridgetype, 'simplysupp')
    for i=1:2
        if tribl_sism(i)==0
            beartype{i}='free';
        else
             beartype{i}='fixed';
        end
    end
end

%% Permanent loads

R=(Ppdeck(1).*tribl); %kN, vector of reactions at the bearing/s
Npp=sum(R); %kN, sforzo normale dovuto ai carichi permanenti
Mpp_x=Npp*Ppdeck(2); %kN, momento flettente trasversale dovuto ai carichi permanenti
if strcmp(bridgetype, 'simplysupp')
    Mpp_y=R(1)*d_bear(1)-R(2)*d_bear(2); %kN, momento flettente trasversale dovuto ai carichi permanenti
elseif strcmp(bridgetype, 'cont')
    Mpp_y=0;
end

%% RUN for service, wind, traffic loads

%Call load analysis consistent to bridge load code 1980
[NsACC,MsxACC,MsyACC,Fx_traffic,My_traffic,Fy_wind,Mx_wind]=bridgeloadanalysis90(tribl,ldeck,l_marciapiede,l_sedestrad,h_girder,h_bent,h_slab,Lx_pier,categoria_ponte,beartype,d_bear);

%% Seismic loads
%Call seismic load analysis consistent with seismic code 1975
[Fx_seism,Fy_seism,Fv_seism, My_seism,Mx_seism]=seismloadanalysis75(Npp,h_bent,Pp_bent,Pp_bentcap,tribl_sism,Ppdeck,Lx_pier,Ly_pier,Ss,piertype,beartype,E_pila,Iseismprot);

%% COMBINAZIONI DI CARICO

%Solllecitazioni prodotte da tutte le forze in gioco sulla struttura nel
%seguente ordine che costistuisce il vettore: N,Mx,My,Vx,Vy.
%Si ricordi che tutte le forze e momenti sono state calcolate in valore
%assoluto, in modo da poterle considerare di verso concorde e poterle
%sommare per ottenere l'effetto piu sfavorevole in termini di sollecitazione

%A1 group (Permanent Load)
N_permanent=Pp_bent+Pp_bentcap+Npp;
P=[N_permanent Mpp_x sum(Mpp_y) 0  0];
   
W=[0 Mx_wind(end) 0 0 Fy_wind(end)];
    
A1=P+W;

%AII AND AIII group (Service loads)
emptycol=zeros(length(NsACC),1);
tempACC=[NsACC',MsxACC',MsyACC' emptycol emptycol];
tempWIND=[emptycol, [Mx_wind(1:end-1)';Mx_wind(1:end-1)'], emptycol, emptycol, [Fy_wind(1:end-1)';Fy_wind(1:end-1)']];
tempTRAFFIC=[emptycol emptycol My_traffic' Fx_traffic' emptycol];

A2=P+abs(tempACC)+0.6*tempWIND;
A3=P+abs(tempACC)+0.2*tempWIND+tempTRAFFIC;

%AIV is neglected because curve bridges are not considered

%AV 
tempP=abs(P.*ones(4,1));
seism=[Fv_seism(1) 0 My_seism Fx_seism 0;... %Fx_sisma
     Fv_seism(2) Mx_seism 0 0 Fy_seism;...%Fy_sisma
    -Fv_seism(1) 0 My_seism Fx_seism 0;...%Fx_sisma
    -Fv_seism(2) Mx_seism 0 0 Fy_seism];%Fy_sisma
A5=tempP+seism;

%OUTPUT allocation
comb.c1=A1;
comb.c2=A2;
comb.c3=A3;
comb.c5=A5;

soll=vertcat(A1,A2,A3,A5);

end
