function [Fx_seism,Fy_seism,Fv_seism, My_seism,Mx_seism]=seismloadanalysis75(Npp,h_bent,Pp_bent,Pp_bentcap,tribl_sism,Pp,Lx_pier,Ly_pier,Ss,piertype,appoggio_testapila,E_pila,varargin)

%% Optional input

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
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
 Iseismprot = optargs{:};

%% Azioni sismiche 
%Si calcolano le azioni sismiche statiche
 for i=1:length(appoggio_testapila)
    if strcmp(appoggio_testapila{i},'fixed')
        mass_deck_x(i)=(Pp(1)*tribl_sism(i)); %kN, massa sismica longitudinale afferente al testa pila
    elseif strcmp(appoggio_testapila{i},'free')
        mass_deck_x(i)=0; %kN ,massa sismica longitudinale afferente al testa pila
    end
end
mass_seism=zeros(1,2);
mass_seism(2)=(Npp+(1/3)*Pp_bent+Pp_bentcap)/9.81; %t, massa sismica trasversale afferente al testa pila
mass_seism(1)=((1/3)*Pp_bent+Pp_bentcap+sum(mass_deck_x))/9.81; %t massa sismica in direzione longitudinale

C=max(0,(Ss-2)/100); %coefficiente d'intensità sismica

%Dynamic parameters calculation
if strcmp(piertype, 'RECT')
    J(1)=(Ly_pier*((Lx_pier)^3))/12; %m^4, momento d'inerzia
    J(2)=(Lx_pier*((Ly_pier)^3))/12; %m^4, momento d'inerzia
elseif strcmp(piertype, 'CIRC')
    J=ones(1,2).*pi*Ly_pier^4/64; %m^4, momento d'inerzia
end

K=3*E_pila*10^3.*J/((h_bent)^3);%rigidezza sezione di base pila
%Direzione X
T=2*pi*(mass_seism./K).^0.5; %s, periodo fondamentale del sistema ad 1 g.d.l.
for i=1:2
    if  T(i)>0.8 %s
        R(i)=0.862./(T(i).^0.66);% coefficiente di risposta della struttura
    elseif  T(i)<= 0.8
        R(i)=1;
    end
end

Fx_seism=C*Iseismprot*R(1)*mass_seism(1)*9.81; %kN, Forza sismica longitudinale dell'impalcato applicata in testa alla pila
Fy_seism=C*Iseismprot*R(2)*mass_seism(2)*9.81;  %kN, Forza sismica trasversale dell'impalcato applicata in testa alla pila
Fv_seism=2*Iseismprot*C*mass_seism*9.81; %kN,  Forza sismica verticale dell'impalcato applicata in testa alla pila
My_seism=(Fx_seism*h_bent);%kNm, momento flettente longitudinale generato dal sisma alla base della pila
Mx_seism=(Fy_seism*h_bent);%kNm, momento flettente trasversale generato dal sisma alla base della pila

end