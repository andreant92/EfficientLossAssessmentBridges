function [deltaSDOF, VbaseSDOF, csiSDOF, ethaSDOF, massSDOF, EDPds, CrMemb, displBRIDGEsub, shearBRIDGEsub] = DBAssSeries(FDpiers, FDbears, Mdeck, Mpiers, DSthresh, varargin)

%% Optional input
numvarargs = length(varargin);
if numvarargs > 2
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 2 optional inputs');
end
% set defaults for optional inputs
optargs = {[0 0 0 0], 'noplot'};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[ backfInt, plotter ] = optargs{:};

%% INPUT ARRANGEMENT

for pr=1:length(FDpiers)
    %Avoid perfect elasto-plastic situation
    FDbears{pr}(end,end)=FDbears{pr}(end,end)*1.002;
    FDpiers{pr}(end,end)=FDpiers{pr}(end,end)*1.002;
end

%% SDOF model calculation
% RUN

for pr=1:length(FDpiers)
    
    masses=[Mdeck(pr) Mpiers(pr)]/9.81; %Masses vector
    
    curvePier=FDpiers{pr};curveBear=FDbears{pr};
    %Curve elab (avoiding softening)
    if curvePier(end,2)<=curvePier(2,2)
        curvePier(:,2)=[0 curvePier(2,2) curvePier(2,2)*1.01];
    end
    
    forcestep=0:min([curvePier(end,2);curveBear(end,2)])/1000:min([curvePier(end,2);curveBear(end,2)]);
    dpier=interp1(curvePier(:,2),curvePier(:,1),forcestep);
    dbear=interp1(FDbears{pr}(:,2),FDbears{pr}(:,1),forcestep);
    
    step=0.001; delta=0:step:max(dpier+dbear)*3; %%%% ho aggiunto *3
    deltaSUB(1,:)=delta;
    shearSUB=interp1(dpier+dbear,forcestep,deltaSUB(1,:),'linear','extrap');
    deltaSUB(2,:)=interp1(dpier+dbear,dpier,deltaSUB(1,:),'linear','extrap');
    
    clear delta forcestep
    %% DS definition
    
    deltaMEMB=[deltaSUB(1,:)-deltaSUB(2,:); deltaSUB(2,:)]; %displ of each member
    
    %Identify DS on the Individual Pier Capacity Curve
    for ds=1:size(DSthresh{1},2)
          %Bears
        DStopSub(1,ds)=interp1(deltaMEMB(1,:),deltaSUB(1,:),DSthresh{1}(pr,ds),'linear','extrap');
        %Piers
        DStopSub(2,ds)=interp1(deltaMEMB(2,:),deltaSUB(1,:),DSthresh{2}(pr,ds),'linear','extrap');
        %Critical memb
        [DStopDispl{pr}(ds),indmemb(ds)]=min(abs(DStopSub(:,ds)));
        CrMemb{pr,ds}=strrep(strrep(num2str(indmemb(ds)),'1','BEAR'),'2', 'PIER');
    end
    
    %% eqSDOF characterisation
    
    deltaSDOF{pr} = sum(masses'.*(deltaSUB.^2))./sum(masses'.*deltaSUB);
    VbaseSDOF{pr} =  shearSUB;
    
    %EDP calculation
    EDPds(pr,:)=interp1(deltaSUB(1,:), deltaSDOF{pr}, DStopDispl{pr},'linear','extrap');
    
    % pier ductility demand
    mu{pr} = max(ones(2,1), deltaMEMB ./ [FDbears{pr}(2,1); FDpiers{pr}(2,1)]); % [adim]  Domanda di duttilità sui subassemblies
    
    % pier and abutments effective damping
    csiPier{pr} = 0.05 + (0.444.*((mu{pr}(2,:)-1)./(mu{pr}(2,:).*pi))); %EVD piers (Takeda)
    
    if backfInt(pr)==0
        csiBear{pr} = 0.05 + (0.670.*((mu{pr}(1,:)-1)./(mu{pr}(1,:).*pi))); %EVD bearings (EPP)
    elseif backfInt(pr)==1
        csiBear{pr} = 0.05*ones(1,length(mu{pr})); %EVD bearings (EPP)
    end
    
    %Subsystem effective damping
    csiSDOF{pr}=(csiBear{pr}.*shearSUB.*deltaMEMB(1,:)+csiPier{pr}.*shearSUB.*deltaMEMB(2,:))./...
        (shearSUB.*deltaMEMB(1,:)+shearSUB.*deltaMEMB(2,:));
    ethaSDOF{pr}=(0.07./(0.02+csiSDOF{pr})).^0.5;
    
    %Effective mass
    massSDOF{pr}=sum(masses'.*deltaSUB)./ deltaSDOF{pr};
        
    displBRIDGEsub{pr}=deltaSUB;
    shearBRIDGEsub{pr}=shearSUB;
    
    clear deltaSUB shearSUB
end

%% plot
if strcmp(plotter,'plot')
    %plot1
    mkrDS={'^','d','o'};
%     style={'--','-','-','-','--','--'}
    col=linspace(0.2,0.8, length(FDpiers));
    figure
    hold 'on'
    for pr=1:length(FDpiers)
        plot(deltaSDOF{pr},VbaseSDOF{pr},'color',col(pr)*ones(1,3),'linewidth',1.5)
        for ds=1:length(EDPds(pr,:))
            scatter(EDPds(pr,ds), interp1(deltaSDOF{pr}(2:end),VbaseSDOF{pr}(2:end),EDPds(pr,ds),'linear','extrap'),30,col(pr)*ones(1,3),mkrDS{ds},'filled')
        end
        ylim([0 6000])
    end
    DS1=scatter(NaN,NaN,30,col(1)*ones(1,3),mkrDS{1},'filled','displayname','DS1')
    DS2=scatter(NaN,NaN,30,col(1)*ones(1,3),mkrDS{2},'filled','displayname','DS2')
    DS3=scatter(NaN,NaN,30,col(1)*ones(1,3),mkrDS{3},'filled','displayname','DS3')

 PIER=plot(NaN,NaN,'color',col(2)*ones(1,3),'linewidth',1.5,'Displayname','Pier Subass')
 ABUT=plot(NaN,NaN,'color',col(2)*ones(1,3),'linewidth',1.5,'Displayname','Abut Subass')

legend([DS1,DS2,DS3, PIER,ABUT])
    
    %plot2

    for pr=1:length(FDpiers)
        figure        
        hold on
        plot(FDpiers{pr}(:,1),FDpiers{pr}(:,2),'k--','displayname','Pier')
        plot(FDbears{pr}(:,1),FDbears{pr}(:,2),'k-.','displayname','Bearing system')
        %         plot(displVECT{pr}(1,:),VbaseSDOF{pr},'displayname','strategy ANdrea')
        plot(deltaSDOF{pr},VbaseSDOF{pr},'k-','linewidth',1.5,'displayname','2DoF Model')
        legend
        ylim([0 max(VbaseSDOF{pr}(1,:))*1.5])
    end
end
end

