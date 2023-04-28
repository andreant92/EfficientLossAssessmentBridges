function [Fx_seism,Fy_seism,Fv_seism, My_seism,Mx_seism]=seismloadanalysis37(Npp,ldeck,h_bent,Pp_bent,Pp_bentcap,tribl_sism,Pp,appoggio_testapila,cat_seism,varargin)

%% Optional input

% Maximum number of optional inputs
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 1 optional inputs');
end

% set defaults for optional inputs
optargs = {[0.1 0.05]};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
 Coeffs = optargs{:};

%% Azione sismica 
%I coefficienti devono essere inseriti in base all'anno es. RD 1935 C=0.1
%(1a) e C=0.05 (2a)

C=Coeffs(1);%rapporto tra le forze orizzontali ed i pesi corrispondenti su cui agiscono per le zone di 1 categ.
C1=Coeffs(2);%rapporto tra le forze orizzontali ed i pesi corrispondenti su cui agiscono per le zone di 2 categ.

%Si calcolano le azioni sismiche statiche
for i=1:length(appoggio_testapila)
    if strcmp(appoggio_testapila{i},'fixed')
        mass_seism_x(i)=(Pp(1)*tribl_sism(i)); %kN, massa sismica longitudinale afferente al testa pila
    elseif strcmp(appoggio_testapila{i},'free')
        mass_seism_x(i)=0; %kN ,massa sismica longitudinale afferente al testa pila
    end
end

mass_seismY=Npp+(1/3)*Pp_bent+Pp_bentcap; %kN, massa sismica trasversale afferente al testa pila

mass_seismX=(1/3)*Pp_bent+Pp_bentcap+sum(mass_seism_x);

if cat_seism==1
    Fx_seism=C*mass_seismX;%kN, Forza sismica longitudinale dell'impalcato applicata in testa alla pila
    Fy_seism=C*mass_seismY;%kN, Forza sismica trasversale dell'impalcato applicata in testa alla pila
    Fv_seism=0.14*mass_seismY; %kN,  Forza sismica verticale dell'impalcato applicata in testa alla pila
elseif cat_seism==2
    Fx_seism=C1*mass_seismX;%kN, Forza sismica longitudinale dell'impalcato applicata in testa alla pila
    Fy_seism=C1*+mass_seismY;%kN, Forza sismica trasversale dell'impalcato applicata in testa alla pila
    Fv_seism=0.14*mass_seismY; %kN,  Forza sismica verticale dell'impalcato applicata in testa alla pila
end
My_seism=(Fx_seism*h_bent);%kNm, momento flettente longitudinale generato dal sisma alla base della pila
Mx_seism=(Fy_seism*h_bent);%kNm, momento flettente trasversale generato dal sisma alla base della pila    

end