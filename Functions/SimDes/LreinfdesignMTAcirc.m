function [longbars,RHOfinal]=LreinfdesignMTAcirc(year,D,Ns,Msy,Msx,sigmac_adm,sigmas_adm,varargin)

%Progetto a pressoflessione retta di una sezione rettangolare usando domini
%di resistenza semplificati e MTA
%% Optional input

% Maximum number of optional inputs
numvarargs = length(varargin);
if numvarargs > 3
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

% set defaults for optional inputs
optargs = {0, 0.08, 15};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[plotter, cover, n] = optargs{:};

%% Input definition

k=sigmas_adm/sigmac_adm;
lam=cover/(D/2);

Nad=Ns/((D/2)^2*pi*sigmac_adm);
Mady=Msy/((D/2)^3*pi*sigmac_adm);
Madx=Msx/((D/2)^3*pi*sigmac_adm);

[rho(1), resdom{1}, groupdom{1}] = MTAintdomADIMcirc(k, Nad, Mady, 'noplot', lam, n);
[rho(2), resdom{2}, groupdom{2}] = MTAintdomADIMcirc(k, Nad, Madx, 'noplot', lam, n);

rhoreq=max(rho); %rho required

clear Nad Madx Mady

if any(rho==0)
    RHOfinal=0;
    longbars=0;
    return
end
%% Check minimum

Ac_min=(max(Ns))/(0.7*sigmac_adm); %mq, area dicalcestruzzo strettamente necessaria per sezione
%soggetta a N centrato.

%%%%Limiti normativi
As_tot=rhoreq.*pi*(D/2)^2; %mq,armatura longitudinale totale della sezione

if year <1972
        % Calcolo moltiplicatore Ac_min attraverso l'interpolazione lineare
        m=interp1([1.999 2 8 8.0001],[1 0.008 0.005 0],Ac_min); %1.999 and 8.0001 in order to avoid unique point
        Asmin=m*Ac_min;
        Asmax=[];
elseif year >= 1972 && year<1974 %%DM 1972
        Asmin=max([0.006*Ac_min,0.003*pi*(D/2)^2]);
        Asmax=0.05*Ac_min;
elseif year >=1974 && year< 1980 %%DM 1974
        Asmin=max([0.006*Ac_min,0.003*pi*(D/2)^2]);
        Asmax=0.05*pi*(D/2)^2;
elseif year >= 1980 
        Asmin=max([0.008*Ac_min,0.003*pi*(D/2)^2]);
        Asmax=0.06*pi*(D/2)^2;
end

%Update Astot if it is not enough to respect minimum.
if As_tot<Asmin
    rhoreq=Asmin/(pi*(D/2)^2);
    As_tot=Asmin; %mq,armatura longitudinale totale della sezione
end

%Check on Astot if it is higher than the maximum
if As_tot>Asmax
    RHOfinal=0;
    longbars=0;
    return
end

%% Hypotesize real configration of bars

%trial configuration

fi=[16 18 20 22 24 26 28]*10^(-3);%[m] tipolgia armatura della sezione
int=[0.25 0.2 0.15 0.10];%[m] interferri usuali dell'anno di costruzione ossia 10,15,20,25 cm
% 0.075 e 0.05 corrispondono al doppio strato
A_fi=(pi.*(fi/2).^2); %[mq]area sezione ferri

nf=floor((2*pi*(D/2))./int);%[-] numero ferri dell'armatura longitudinale soggetta a trazione
%o compressione disposti lungo la dir index
RHOvect=[(A_fi'*nf)./(pi*(D/2)^2)]'; %[-]
%     [RHOsorted,~]=sort(RHOvect(:));
%Choose the likely configuration
diff=(RHOvect(:)-rhoreq);
good=find(diff>0); [~,id]=min(diff(good)); rightcomb=(good(id));

if isempty(good)
    warning('no realistic reinforcement configuration for this geom/mech sample')
    RHOfinal=0;
    longbars=0;
else
    [x,y]=find(RHOvect(rightcomb)==RHOvect);
    RHOfinal=RHOvect(rightcomb); %%
    longbars(:)=[fi(y), nf(x) ];
end


%% PLOT
if strcmp(plotter,'plot')
    Ms={Msy,Msx};
    tit={'dir long','dir trasv'};
    figure
    for j=1:2
        subplot(1,2,j)
        hold on
        for i=1:length(groupdom{j})
            plot( groupdom{j}{i}(2,:).*pi*(D/2)^3.*sigmac_adm, groupdom{j}{i}(1,:).*pi*(D/2)^2.*sigmac_adm,'color',[.5 .5 .5])
        end
        scatter(Ms{j}(1:end-4),Ns(1:end-4),15,'r','filled')
        scatter(Ms{j}(end-3:end),Ns(end-3:end),15,'b','filled')
        plot(resdom{j}(2,:).*pi*(D/2)^3.*sigmac_adm,resdom{j}(1,:).*pi*(D/2)^2*sigmac_adm,'k','Linewidth',3)
        title(tit{j}); xlabel('N[kN]'); ylabel('M[kNm]')
    end
end
end