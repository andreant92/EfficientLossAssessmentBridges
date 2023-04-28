function [longbars, trasvbars]=MTAsimdesignCIRC(year,D,Nsd,Msd_x,Msd_y,Vsd_x,Vsd_y,sigmac_r,sigmas_r, varargin)

%Progettazione a Pressoflessione e taglio della sezione rettangolare con MTA

%Convenzione assi: asse x è l'asse longitudinale dell'impalcato
%                  asse y è l'asse trasversale dell'impalcato

% Age is an index of year of construction, it determines the code to be
% used
%S= 1 pre-1957
%S= 2 || 3 || between 1957 and 1979
%S= 4 || 5 || 6 post-1980

%% Optional input

% Maximum number of optional inputs
numvarargs = length(varargin);
if numvarargs > 2
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

% set defaults for optional inputs
optargs = {0.07, 15};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[cover, n] = optargs{:};


%% Elaborazione degli input
%Calcolo sigma ammissibile del cls
% sigmac_r=sigmac_r*10;
if year <= 1972
        if sigmac_r*10>=250 %kg/cmq
            sigmac_adm= (75+(sigmac_r*10-225)/9); %kg/cmq tensione ammissibile del conglomerato
        elseif sigmac_r*10<250
            sigmac_adm=sigmac_r*10/3; %kg/cmq
        end
else
        sigmac_adm= 60+((sigmac_r*10-150)/4); %kg/cmq, tensione di compressione ammissibile del conglomerato
end

%calcolo sigmas_adm
if year >= 1939 && year<1972 %%RD 1939
        if sigmas_r*10>=2300 && sigmas_r*10<2700 %kg/cmq
            sigmas_adm= 1400; %kg/cmq tensione ammissibile del conglomerato
        elseif sigmas_r*10>=2700 && sigmas_r*10<3100
            sigmas_adm=1600; %kg/cmq
        elseif sigmas_r*10>=3100
            sigmas_adm=1800; %kg/cmq
        end
elseif year >= 1972 && year<1974 %%DM 1972
        if sigmas_r*10>=2200 && sigmas_r*10<3200 %kg/cmq
            sigmas_adm= 1200; %kg/cmq tensione ammissibile del conglomerato
        elseif sigmas_r*10>=3200 && sigmas_r*10<3800
            sigmas_adm=1600; %kg/cmq
        elseif sigmas_r*10>=3800 && sigmas_r*10<4100
            sigmas_adm=2200; %kg/cmq
        elseif sigmas_r*10>=4100 && sigmas_r*10<4400
            sigmas_adm=2400; %kg/cmq
        elseif sigmas_r*10>=4400
            sigmas_adm=2600; %kg/cmq
        end
elseif year >=1974 %%DM 1974
        if sigmas_r*10>=2200 && sigmas_r*10<3200 %kg/cmq
            sigmas_adm= 1200; %kg/cmq tensione ammissibile del conglomerato
        elseif sigmas_r*10>=3200 && sigmas_r*10<3800
            sigmas_adm=1600; %kg/cmq
        elseif sigmas_r*10>=3800 && sigmas_r*10<4400
            sigmas_adm=2200; %kg/cmq
        elseif sigmas_r*10>=4400
            sigmas_adm=2600; %kg/cmq
        end 
end

plotter='plot'; %write plot or char to plot results

% some conversions
sigmac_adm=sigmac_adm*10^2;%[kg/cmq to kN/mq]
sigmas_adm=sigmas_adm*10^2;%[kg/cmq to kN/mq]

%% Progetto dell'armatura longitudinale 

[longbars, ~]=LreinfdesignMTAcirc(year,D,Nsd,Msd_y,Msd_x,sigmac_adm,sigmas_adm,cover,n);
if any(longbars==0)
    trasvbars=0;
    warning('this geom/mech configuration is not adequate to this stress')
    return
end

%% Progetto a taglio, con forza di taglio
% limiti min e max da norma
if year<1972 %%RD 1939
        if sigmac_r*10<=120% valida  per  conglomerati di cemento idraulico normale(Portland), d'alto forno e pozzolanico.
            tau_adm=[4 14]; %kg/cmq, per conglomerai di cemento idraulico normale(Portland), d'alto forno e pozzolanico.
        elseif sigmac_r*10>=160 %valida per conglomerati di cemento ad alta resistenza ed alluminoso.
            tau_adm=[6 14]; %kg/cmq
        end
else
        tau_adm=[(4+((sigmac_r*10-150)/75)) (14+((sigmac_r*10-150)/35))]; %kg/cmq, per conglomerai di cemento idraulico normale(Portland), d'alto forno e pozzolanico.
end
tau_adm=tau_adm*10^2; %conversion in [kN/mq]

%fisso tre passi di staffatura usuali dell'anno di costruzione, questa scelta verifica la condizione che ci siano almeno 3 staffe al metro
passo_staffe=[.25 .20 .15];

%Definisco i passi di tentativo
index=passo_staffe(:)>10*longbars(1,1) | passo_staffe(:)>10*longbars(2,1) | passo_staffe(:)>0.5*min(Lx_pier,Ly_pier);
passo_staffe(index)=[];

%% Progetto a taglio,con forza di taglio lungo x

Ved(1)=max(Vsd_x);%kN
Ved(2)=max(Vsd_y);%kN

trasvbars=TreinfdesignMTA(year,Lx_pier,Ly_pier,Ved,tau_adm,sigmas_adm,passo_staffe,cover);

end
