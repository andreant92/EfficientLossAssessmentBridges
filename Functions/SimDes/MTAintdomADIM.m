function [rho, resdom, groupdom] = MTAintdomADIM(k, Nsad, Msad, plotter, varargin)
% This function builds simplified interaction domain (M-N) for symmetrically
% reinforced square section
addpath('C:\Users\netti\OneDrive - Politecnico di Bari\Documenti\MATLAB')
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
optargs = {0.03, 15};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[lam, n] = optargs{:};

%% Input elaboration
rangeRho=0.0002:0.0002:0.01;

%% Analytical domain costruction
% ADAPTED With reference to "Aurelio Ghersi, Dalle tensioni ammissibili
%agli stati limite: un approccio unitario"
% adimens def
Nrad=zeros(length(rangeRho),3); %preallocation
Mrad=zeros(length(rangeRho),3); %preallocation

for i=1:length(rangeRho)
    %A
    Nrad(i,3)=((1+lam)+2*rangeRho(i)*n); %[-]
    %%B
    Nrad(i,2)= n/(2*(k+n))-rangeRho(i)*(k+k*lam-n+n*lam); %[-]
    Mrad(i,2)= n/(2*(k+n))*((1+lam)/2-1/3*(n/(k+n)))+...
        rangeRho(i)/2*(1-lam)^2*(k+n);    %[-]
    %%C
    Nrad(i,1)=-2*rangeRho(i)*k;%[-]
end

%% RUN
iter=0;
conv=0;
distS=((Msad.^2)+(Nsad.^2)).^(1/2);% raggio vettore,distnza O-(Msd,Nsd)

while(conv==0 && iter<length(rangeRho))
    iter=iter+1;
    
    for i=1:length(Nsad)
        %Calcolo distanza
        sollx=[0 Msad(i)].*100; %retta passante per il punto i
        solly=[0 Nsad(i)].*100; %retta passante per il punto i
        [x5(i),y5(i)]=intersections(Mrad(iter,:),Nrad(iter,:),sollx,solly);% P=intersezione retta passante
        %per i punti (0,Nrd_max) e (Mrd_max,0) e retta passante per O(0,0) e (Msd,Nsd).
    end
    
    distR(iter,:)=((x5.^2)+(y5.^2)).^(1/2);%distanza punti O-P
    
    %Il ciclo terminerà quando conv=1, affinchè sia vero ènecessario verificare la seguente condizione:
    %tutti i punti  rappresentativi lo stato di sollecitazione della sezione devono trovarsi all'interno del
    %dominio d'interazione ottenuto dopo n cicli che li contiene .
    
    conv=all((distR(iter,:)-distS')>0);
end

%output
rho=rangeRho(iter);
resdom=[Nrad(iter,:);Mrad(iter,:)];
for i=1:length(Nrad)
    groupdom{i}=[Nrad(i,:);Mrad(i,:)];
end

%TRial plot
if strcmp(plotter, 'plot')
    figure
    hold on
    for i=1:length(Nrad)
        plot( Mrad(i,:), Nrad(i,:),'color',[.5 .5 .5])
    end
    scatter(Msad,Nsad,15,'r','filled')
    plot(resdom(2,:),resdom(1,:),'k','Linewidth',2)
    title(['k=',num2str(k,2),' n=',num2str(n), ' \lambda=', num2str(lam,2)])
end
% warning if no solution
if iter==length(rangeRho)
%     gg=msgbox('Warning:no solution for this geometry and stress combination','ERRORE', 'error');
    warning('no solution for this geometry and stress combination')
    rho=0;
    resdom=0;
end
end