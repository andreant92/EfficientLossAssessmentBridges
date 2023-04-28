function [deltaSDOF, VbaseSDOF, csiSDOF, ethaSDOF, massSDOF, TeffSDOF, EDPds, CrMemb, displBRIDGEsub, shearBRIDGEsub ] = DBAssTiterSpeed(FDpiers, FDbears, Mdeck, ...
    Mpiers, DSthresh, varargin)

%% Optional input

% Maximum number of optional inputs
numvarargs = length(varargin);
if numvarargs > 2
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

% set defaults for optional inputs
optargs = {'noplot', 1*ones(1,length(FDpiers)),'modal'};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
[ plotter, deltamax ] = optargs{:};

%% INPUT ARRANGEMENT

for pr=1:length(FDpiers)
    %Avoid perfect elasto-plastic situation
    FDbears{pr}(end,end)=FDbears{pr}(end,end)*1.001;
    FDpiers{pr}(end,end)=FDpiers{pr}(end,end)*1.001;
end
%Pre-allocation
deltastep=0.0001;
nmaxiter=100;
bilinearisation=0;

%% RUN
for pr=1:length(FDpiers)
    
    deltacontrol=0:deltastep:deltamax(pr); %Incremental displ load step
    coeffshape=[Mdeck(pr)/(Mdeck(pr)+Mpiers(pr)) Mpiers(pr)/(Mdeck(pr)+Mpiers(pr))]; %initial shape of the deformed profile
    masses=[Mdeck(pr) Mpiers(pr)]/9.81; %Masses vector
    FDsubass={FDbears{pr}, FDpiers{pr}}; %FD cell array
    deltaSUB=zeros(2,length(deltacontrol));
    shearSUB=zeros(2,length(deltacontrol));
    
%     deform1=coeffshape/max(coeffshape)'*deltacontrol(2);
%     deform1=deform1';
%     for d=2:length(deltacontrol)
%         
%         deltac=deltacontrol(d);
%         %Tentative profile
%         deform1=normalize(deform1)*deltac;
%         
%         conv=0; toll=0.001; iter=0; %Initialize Control parameters for the while-loop
%         while conv==0 && iter<nmaxiter
%             
%             iter=iter+1;
%             deform=deform1;
%             displMemb(1)=deform(1)-deform(2); %displ of the bear
%             displMemb(2)=deform(2); %displ of the pier
%             
%             %Find effective stiffnesses
%             for i=[1 2]
%                 Vi(i)=interp1(FDsubass{i}(:,1), FDsubass{i}(:,2), displMemb(i),'linear','extrap');
%                 effStiff(i)=Vi(i)/displMemb(i);
%             end
%             
%             %------------%%%%%%%%%%%   MODAL ANALYSIS
%             %Mass Matrix
%             M=[masses(1) 0; 0 masses(2)];
%             %Stiffness MAtrix
%             K=[effStiff(1) -effStiff(1); -effStiff(1) sum(effStiff)];
%             % Perform Modal analysis
%             [per, modalshapes, masspart, ~] = MODALanalysisGeneral(K, M);
%             %------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%             deform1=modalshapes(:,end).*deltac; %define new deformed shape
%             %             Vmemb(1)=interp1(FDsubass{1}(:,1), FDsubass{1}(:,2),deform1(1)-deform1(2),'linear','extrap');
%             %             Vmemb(2)=interp1(FDsubass{2}(:,1), FDsubass{2}(:,2), deform1(2),'linear','extrap');
%             
%             conv=all(abs(deform1-deform)<toll); %Calculate convergence
%         end
%         
%         Vmemb(1)=interp1(FDsubass{1}(:,1), FDsubass{1}(:,2),deform1(1)-deform1(2),'linear','extrap');
%         Vmemb(2)=interp1(FDsubass{2}(:,1), FDsubass{2}(:,2), deform1(2),'linear','extrap');
%         
%         deltaSUB(:,d)=deform1;
%         shearSUB(:,d)=[Vmemb(2), Vmemb(1)];
%     end

    
%%%%%%%%%%     Check strategy
%%%%%%         comparison with sdof model

    curvePier=FDpiers{pr};curveBear=FDbears{pr};    
    %Curve elab (avoiding softening)
    if curvePier(end,2)<=curvePier(2,2)
        curvePier(:,2)=[0 curvePier(2,2) curvePier(2,2)*1.01];
    end
    
    forces=0:min([curvePier(end,2);curveBear(end,2)])/10000:min([curvePier(end,2);curveBear(end,2)]);
    dpier=interp1(curvePier(:,2),curvePier(:,1),forces);    
    dbear=interp1(FDbears{pr}(:,2),FDbears{pr}(:,1),forces);
    
    %delta top
    dtop=dbear+dpier;
    
%     dbear=interp1(FDbears{pr}(:,2),FDbears{pr}(:,1),800);
%     dpier=interp1(FDpiers{pr}(:,2),FDpiers{pr}(:,1),800,'linear','extrap');    
    deltaSUBtemp=[dtop;dpier];
    shearSUBtemp=[forces;forces];

    %Re-sampling
    deltaSUB(1,:)=deltacontrol;
    deltaSUB(2,:)=interp1( deltaSUBtemp(1,:),deltaSUBtemp(2,:), deltacontrol,'linear','extrap');
    shearSUB(1,:)=interp1( deltaSUBtemp(1,:),shearSUBtemp(1,:), deltacontrol,'linear','extrap');
    shearSUB(2,:)=interp1( deltaSUBtemp(1,:),shearSUBtemp(2,:), deltacontrol,'linear','extrap');

%     figure
%     hold on
%     plot(FDpiers{pr}(:,1),FDpiers{pr}(:,2),'displayname','pier')
%     plot(FDbears{pr}(:,1),FDbears{pr}(:,2),'displayname','bear')
%     %plot(dtop,forces,'displayname','strategy Cardone')
%     plot(deltaSUB(1,:),shearSUB(1,:),'displayname','strategy ANdrea')
%     %plot(deltaSUB2(1,:),shearSUB2(1,:),'displayname','strategy CArdone')
%     legend
%     ylim([0 max(shearSUB(1,:))*1.5])
%         
    %% DS definition
    
    deltaMEMB=[deltaSUB(1,:)-deltaSUB(2,:); deltaSUB(2,:)]; %displ of each member %[bear; bent];
    
    %Identify DS on the Individual Pier Capacity Curve
    for ds=1:size(DSthresh{1},2)
          %Bears
        DStopSub(1,ds)=interp1(deltaMEMB(1,:),deltaSUB(1,:),DSthresh{1}(pr,ds));%,'linear','extrap');
        %Piers
        DStopSub(2,ds)=interp1(deltaMEMB(2,:),deltaSUB(1,:),DSthresh{2}(pr,ds));%,'linear','extrap');
        %Critical memb
        [DStopDispl{pr}(ds),indmemb(ds)]=min(abs(DStopSub(:,ds)));
        CrMemb{pr,ds}=strrep(strrep(num2str(indmemb(ds)),'1','BEAR'),'2', 'PIER');

%         %Bears
%         ind(1,ds)=find(deltaMEMB(1,:)<=DSthresh{1}(pr,ds),1,'last');
%         %Piers
%         ind(2,ds)=find(deltaMEMB(2,:)<=DSthresh{2}(pr,ds),1,'last');
%         %Critical memb
%         [dsstep{pr}(ds),indmemb(ds)]=min(ind(:,ds));
%         CrMemb{pr,ds}=strrep(strrep(num2str(indmemb(ds)),'1','BEAR'),'2', 'PIER');
      end
    
    %% eqSDOF characterisation
    
%     deltaTemp=[deltaSUB(1,:); deltaMEMB(1,:)]; 
    deltaSDOF{pr} = sum(masses'.*(deltaSUB.^2))./sum(masses'.*deltaSUB);
    VbaseSDOF{pr} =  shearSUB(1,:);
    
%     figure
%     hold on 
%     plot(deltaSDOF{pr},VbaseSDOF{pr})
%     plot(deltaSUB(1,:), shearSUB(1,:))
    
    %EDP calculation
    EDPds(pr,:)=interp1(deltaSUB(1,:), deltaSDOF{pr}, DStopDispl{pr},'linear','extrap');
    
    %bilinearisation
    if bilinearisation==1
        temp=[deltaSDOF{pr};VbaseSDOF{pr}];
        bilin=bilinearization(temp', EDPds(pr,1),  EDPds(pr,4), 'plot');
        deltaSDOF{pr} = bilin(:,1);
        VbaseSDOF{pr} = bilin(:,2);
    end
    
    % pier ductility demand
    mu{pr} = max(ones(2,1), deltaMEMB ./ [FDsubass{1}(2,1); FDsubass{2}(2,1)]); % [adim]  Domanda di duttilità sui subassemblies
    
    % pier and abutments effective damping
    csiPier{pr} = 0.05 + (0.444.*((mu{pr}(2,:)-1)./(mu{pr}(2,:).*pi))); %EVD piers (Takeda)
    csiBear{pr} = 0.05 + (0.670.*((mu{pr}(1,:)-1)./(mu{pr}(1,:).*pi))); %EVD bearings (EPP)
    
    %Subsystem effective damping
    csiSDOF{pr}=(csiBear{pr}.*shearSUB(2,:).*deltaMEMB(1,:)+csiPier{pr}.*shearSUB(1,:).*deltaMEMB(2,:))./...
        (shearSUB(2,:).*deltaMEMB(1,:)+shearSUB(1,:).*deltaMEMB(2,:));
    ethaSDOF{pr}=(0.07./(0.02+csiSDOF{pr})).^0.5;
    
    %Effective mass
    massSDOF{pr}=sum(masses'.*deltaSUB)./ deltaSDOF{pr};
    KeffSDOF{pr}= VbaseSDOF{pr}./deltaSDOF{pr};
    TeffSDOF{pr}= 2*pi*(massSDOF{pr}./KeffSDOF{pr}).^0.5;
    
    displBRIDGEsub{pr}=deltaSUB;
    shearBRIDGEsub{pr}=shearSUB;
  
end

%% plot
if strcmp(plotter,'plot')
    figure
    hold 'on'
    for pr=1:length(FDpiers)
        plot(deltaSDOF{pr},VbaseSDOF{pr})
        EDPdsSH(pr,:)=interp1(deltaSDOF{pr}(2:end),VbaseSDOF{pr}(2:end),EDPds(pr,:));
        scatter(EDPds(pr,:), EDPdsSH(pr,:) )
    end
    xlim([0 max(EDPds(2:3,end))*1.5])
        ylim([0 max(EDPdsSH(2:3,end))*1.5])  
end
end

function  bilinear = bilinSUBASS( curva, DISPmechanism, PLOT )
% bilinEURO Calcola un'approssimazione bilineare secondo EuroCodice 8.

%% Function
ind=find(curva(1,:)<=DISPmechanism);
%area sottesa dalla curva di partenza (fino a DISPmechanism)
E = trapz(curva(1,ind),curva(2,ind));

% snervamento della bilineare
forceY = max(curva(2,:));
deltaY = 2 * (DISPmechanism - E / forceY);

% curva bi-lineare
bilinear = [   0               0;                deltaY  forceY; curva(end,1)  	forceY];

%% plot di controllo

if strcmp(PLOT,'plot')
    figure
    plot(curva(1,:),curva(2,:))
    hold on
    plot(bilinear(:,1),bilinear(:,2))
end

end