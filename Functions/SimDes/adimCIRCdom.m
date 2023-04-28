function domain =adimCIRCdom(rho, knum, plotter, varargin)
% This function calculates interaction domain (adimensional) for circular
% cross section.

%% Optional input

% Maximum number of optional inputs
numvarargs = length(varargin);
if numvarargs > 2
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

% set defaults for optional inputs
optargs = {0.04, 15};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[lam, nmat] = optargs{:};

% %% Example 
% 
% % Data of a generic case
% nf=150; %number of long bars
% fi=0.01; %[m] diameter of long bars
% sigmas=255*10^3; sigmac=8.5*10^3; %[kN/mq] steel and concrete admissible stresses
% copr=0.08; %[m] Cover layer 
% nmat=15; 
% rnum=1.000; %[m] clear height
% dnum=2*rnum-copr; %[m] clear height
% As=nf*pi*(fi/2)^2; %[mq] area of steel reinforcements
% rho=As/(pi*rnum^2); %[-] geometric ratio of long reinf
% 
% % INPUT trasformation for this function
% knum=sigmas/sigmac; %[-] coeff accounting for mechanical parameters
% lam=copr/rnum; %[-] coeff accounting for cover layer
% anum=nmat/(knum+nmat)*(2-lam); %[-] coeff accounting for neutral axis depth (depending on mech parameters and cover layer)

%% Defaults
anum=nmat/(knum+nmat)*(2-lam); %[-] coeff accounting for neutral axis depth (depending on mech parameters and cover layer)
FIe=0.93; FIi=0.91;
rhog=((0.02)*2*0.92); %
gammanum=rho/rhog;,

%% Calculate interaction domain 

%%========== A) only tensile stress

nad(1)=-rhog*gammanum*knum;
mad(1)=0;

%%========== C) only compressive stress

nad(3)=1+rhog*gammanum*nmat;
mad(3)=0;

%%========== B) "rottura" bilanciata

syms c a 
% Axial load
Ncls=int((c+a-1)/a*(1^2-c^2)^0.5,c,1-a,1);
Nclsnum=2/pi*double(subs(Ncls,{a},{anum}));

Nsteel=int((c+a-1)/a*((FIe)^2-c^2)^0.5,c,-FIe,FIe)-int((c+a-1)/a*((FIi)^2-c^2)^0.5,c,-FIi,FIi);
Nsteelnum=2/pi*gammanum*nmat*double(subs(Nsteel,{a},{anum}));

nad(2)=Nclsnum+Nsteelnum

%Flexure

Mcls=int((c+a-1)/a*c*(1^2-c^2)^0.5,c,1-a,1);
Mclsnum=2/pi*double(subs(Mcls,{a},{anum}));

Msteel=int((c+a-1)/a*c*(FIe^2-c^2)^0.5,c,-FIe,FIe)-int((c+a-1)/a*c*(FIi^2-c^2)^0.5,c,-FIi,FIi);
% SnumC=2*sigmac*double(subs(temp1,{x R},{xnum, rnum}))
Msteelnum=2/pi*gammanum*nmat*double(subs(Msteel,{a},{anum}));

mad(2)=Mclsnum+Msteelnum
%OUTPUT
domain=[nad;mad];

%% Plot
if strcmpi(plotter, 'plot')
    figure
    hold on
    simpl=plot(domain(2,:),domain(1,:), 'r-')
    scatter(domain(2,:),domain(1,:), 'r','filled')
    legend ([simpl],'Simplified domain')
end

% %% TEST
% 
% DOM(1,:)=domain(1,:)*pi*sigmac*rnum^3;
% DOM(2,:)=domain(2,:)*pi*sigmac*rnum^2;

