function [rho, resdom, groupdom] = MTAintdomADIMcirc(k, Nsad, Msad, plotter, varargin)
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
optargs = {0.03, 15};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[lam, n] = optargs{:};

%% Input elaboration

anum=n/(k+n)*(2-lam); %adimensional neutral axis
FIe=1-(lam-0.01); FIe=round(FIe*10000)/10000; %sup limit of the circular crown of steel uniformly distributed
%Round to the 0.001 to avoid complex number in the solution of the follo
%wing integral
FIi=1-(lam+0.01); FIi=round(FIi*10000)/10000; %sup limit of the circular crown of steel uniformly distributed

rhofk=0.04*(1-lam); %dummy percentage of reinforcements

%Reinforcements
rangeRho=0.0002:0.0005:0.01; %Geometric percentage of reinforc to generate int domains & 0.0021%
rangeGamma=rangeRho/rhofk; %"Density" of geometric percentage of reinforc to generate int domains

%% Analytical domain costruction
% ADAPTED With reference to "Aurelio Ghersi, Dalle tensioni ammissibili
%agli stati limite: un approccio unitario"
% adimens def
Nrad=zeros(length(rangeRho),3); %preallocation
Mrad=zeros(length(rangeRho),3); %preallocation
    
syms c a
for i=1:length(rangeRho)

    %A
    Nrad(i,3)=1+rhofk*rangeGamma(i)*n; %[-]
    %%B    
    % Axial load
    temp=int((c+a-1)/a*(1^2-c^2)^0.5,c,1-a,1); % circular sector of concrete contribution to Nbalanc
    ccnum=2/pi*double(subs(temp,{a},{anum})); %[-]
    
    % circular crown of steel contribution to Nbalanc
    temp1=int((c+a-1)/a*((FIe)^2-c^2)^0.5,c,-FIe,FIe)-int((c+a-1)/a*((FIi)^2-c^2)^0.5,c,-FIi,FIi);
    snums=2/pi*rangeGamma(i)*n*double(subs(temp1,{a},{anum})); %[-]
    
    Nrad(i,2)=(ccnum)+(snums);  %[-]
    
    %Flexure
    % circular sector of concrete contribution to Mbalanc
    temp2=int((c+a-1)/a*c*(1^2-c^2)^0.5,c,1-a,1);
    mcnum=2/pi*double(subs(temp2,{a},{anum}));
    
    % circular crown of steel contribution to Mbalanc
    temp3=int((c+a-1)/a*c*(FIe^2-c^2)^0.5,c,-FIe,FIe)-int((c+a-1)/a*c*(FIi^2-c^2)^0.5,c,-FIi,FIi);
    msnum=2/pi*rangeGamma(i)*n*double(subs(temp3,{a},{anum}));
    
    Mrad(i,2)=(mcnum)+(msnum);  %[-]
    
    clear temp temp1 temp2 temp3
    
    %%C
    Nrad(i,1)=-rhofk*rangeGamma(i)*k;%[-]
end

% N=Nrad*pi*11000*1.3^2;
% M=Mrad*pi*11000*1.3^3;

%% RUN
iter=0;
conv=0;
distS=((Msad.^2)+(Nsad.^2)).^(1/2);% raggio vettore,distnza O-(Msd,Nsd)
x5=zeros(1,length(Nsad)); y5=zeros(1,length(Nsad)); 

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
groupdom=cell(1,length(Nrad));
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