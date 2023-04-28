function [DomainREF, DomainSIMPL] = MTAintdom(h,b, AsT, AsC, sigmac, sigmas, plotter, varargin)
% This function builds simplified interaction domain (M-N) for symmetrically 
% reinforced square section

%% example input
% n=15; %[coeff di omogeneizz] [-]
% c=0.04; [m]
% h=0.6; [m]
% b=0.4; [m]
% sigmac=8500; %[kN/mq]
% sigmas=215000;%[kN/mq]

%% Optional input

% Maximum number of optional inputs
numvarargs = length(varargin);
if numvarargs > 2
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

% set defaults for optional inputs
optargs = {0.05, 15};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[ n, c] = optargs{:};

%% Input elaboration

d=h-c; %clear height
Ac=b*h;%[mq]

%% Analytical domain costruction
% With reference to "Aurelio Ghersi, Dalle tensioni ammissibili
%agli stati limite: un approccio unitario"

% punto A
Na=-(AsT+AsC)*sigmas;
Ma=(AsT-AsC)*(h/2-c)*sigmas;

%punto B
Nb=-(AsT+c/(d)*AsC)*sigmas;
Mb=(AsT+c/(d)*AsC)*(h/2-c)*sigmas;

%punto C
Nc=-(-b*d/2*(sigmac^2/(sigmas/n+sigmac))+(AsT+c/d*AsC)*sigmas-AsC*(d-c)/d*n*sigmac);
Mc=b*d/2*sigmac^2/(sigmas/n+sigmac)*(h/2-sigmac/(sigmas/n+sigmac)*d/3)+...
    ((AsT-c/d*AsC)*sigmas+AsC*(d-c)/d*n*sigmac)*(h/2-c);

%punto D
Nd=-(-b*h/2*sigmac-(AsT*c/h+AsC*d/h)*n*sigmac);
Md=b*h^2/12*sigmac+(-AsT*c/h+AsC*d/h)*(h/2-c)*n*sigmac;

%punto E1
Ne1=-(-(b*h+n*AsT+n*AsC)*0.7*sigmac);
distOG=(n*AsT-n*AsC)/(b*h+n*AsT+n*AsC)*(h/2-c);
sigmacinf=-(0.2*h+distOG)/(0.5*h+distOG)*sigmac;
Me1=b*h^2/12*(sigmac+sigmacinf)+((-AsT*c/h+AsC*d/h)*n*sigmac+...
    (AsT*d/h-AsC*c/h)*n*sigmacinf)*(h/2-c);

%punto E
Ne=-(-(b*h+n*AsT+n*AsC)*0.7*sigmac);
Me=-n*(AsT-AsC)*(h/2-c)*0.7*sigmac;

DomainREF=[Ma, Mb, Mc, Md, Me1, Me; Na, Nb, Nc, Nd, Ne1, Ne];

%% Simplified domain

N=sigmac*(n*(AsC+AsT)+Ac); %[kN]
T=sigmas*(AsC+AsT);
M=AsT*sigmas*0.9*(h-c); %[kN]

DomainSIMPL=[0 M 0; -T 0 N];

%% Plot
if strcmpi(plotter, 'plot')
    figure
    hold on
    ref=plot(DomainREF(1,:),DomainREF(2,:), 'r-')
    simpl=plot(DomainSIMPL(1,:),DomainSIMPL(2,:), 'b-')
    scatter(DomainREF(1,:),DomainREF(2,:), 'r','filled')
    scatter(DomainSIMPL(1,:),DomainSIMPL(2,:), 'b','filled')
    legend ([ref,simpl],'Refined domain','Simplified domain')
end