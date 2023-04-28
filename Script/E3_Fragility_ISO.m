%% 3) Fragility Analysis for simply supported bridges

% V1.0 Updated 25-04-2023
% Author: dr. Andrea Nettis. PhD
% Polytechnic University of Bari. Bari, Italy

% This script is aimed to perform the fragility analysis and loss assessment for
% simply-supported girder bridges (ISO)

%% 
clc
clear
close all 

%Settings
path=cd;
addpath(genpath('./Functions'))

bridgename='B2';

%Insert folder for input and output
folderpath=['./',bridgename]; 
fileinputDBA=[folderpath, '/SimplPerfAss'];
fileinputData=[folderpath, '/SubassembliesInput'];
fileinputHazard=[folderpath,'/HazCurves'];

fileoutput=[folderpath,'/FragilityOutput.mat'];

load(fileinputDBA, 'PerfDisplBearCSMl', 'PerfDisplBearCSMt', 'PerfDisplDeckCSMl', 'PerfDisplDeckCSMt', 'PerfDisplPierCSMl', 'PerfDisplPierCSMt',...
    'DSthreshT','DSthreshL','DBAres')
load(fileinputData)

%Other input
DLRmatrix=[0 .03 .08 .25 1];
% UnitaryCost=930;
% 
% tempData=load(fullfile(folderpath,'vars','BridgeLHSData.mat'),'Data');
% TotalArea(br)=tempData.Data{9,3}*tempData.Data{6,3};
% 
% TotreplCost=TotalArea*UnitaryCost; %[euro]


%% IM definition

%trial Cosenza
load('Spectra.mat')

%%%%% MODAL ANALYSIS to detect 1st period
for cs=1:length(DBAres{1}.teffSDOF)
    for pr=1:length(DBAres{1}.teffSDOF{cs})
        perMatrixT(cs,pr)=DBAres{1}.teffSDOF{cs}{pr}(2);
    end
end

if perMatrixT(1,1)<.15 || perMatrixT(1,end)<.05
    perMatrixT=perMatrixT(:,2:end-1);
end
perMinT=round(100*min(min(perMatrixT))/5)*5/100;
perMinTUp=round(100*min(min(1.5*perMatrixT))/5)*5/100;
perMaxT=round(100*max(max(perMatrixT))/5)*5/100;
rangPerT=[perMinT,perMinTUp]; 

for cs=1:length(DBAres{2}.teffSDOF)
    perMatrixL(cs)=DBAres{2}.teffSDOF{cs}(2);
end
perMinL=round(100*min(perMatrixL)/5)*5/100; 
perMinLUp=round(100*min(1.5*perMatrixL)/5)*5/100; 

rangPerL=[perMinL, perMinLUp]; 

for gm=1:length(ADRSspectra)    
    avgSa(1,gm)=geomean (interp1(ADRSspectra{gm}(:,1), ADRSspectra{gm}(:,2), rangPerT(1):0.05:rangPerT(end)));
    avgSa(2,gm)=geomean (interp1(ADRSspectra{gm}(:,1), ADRSspectra{gm}(:,2), rangPerL(1):0.05:rangPerL(end)));
end
IMcell=avgSa';

clear perMatrixL perMatrixT

%% Arrange DS thresholds and DCR evaluation 

 for sample=1:size(DSthreshL,1)
    %% Find capacity (Matrix with damage states thresholds)
    for ds=1:size(DSthreshL{sample,1},2)
        MatrixDS{ds}=zeros(2,size(DSthreshL{sample,1},1),2);
        for pr=1:size(DSthreshL{sample,1},1)
            MatrixDS{ds}(1,pr,1)=DSthreshT{sample,1}(pr,ds);        %bear transv
            MatrixDS{ds}(2,pr,1)=DSthreshT{sample,2}(pr,ds);        %pier transv
            MatrixDS{ds}(1,pr,2)=DSthreshL{sample,1}(pr,ds);        %bear long
            MatrixDS{ds}(2,pr,2)=DSthreshL{sample,2}(pr,ds);        %pier long
        end
    end
    
    for gm=1:length(ADRSspectra)
        %%
        % 3d matrix having [dim1 (row1): bear deform, row2: pier deform
        % dim2: number of subassemblies, dim3: direction
        
        MatrixEDP=zeros(2,size(DSthreshT{sample,1},1),2);
        % Find peformance edp (displacement) in trasnverse direction
        MatrixEDP(1,:,1)=(PerfDisplBearCSMt{sample}(gm,:));
        MatrixEDP(2,:,1)=PerfDisplPierCSMt{sample}(gm,:);
        
        % Find performance edp (displacement) in longitudinal direction
        MatrixEDP(1,:,2)=(PerfDisplBearCSMl{sample}(gm,:));
        MatrixEDP(2,:,2)=PerfDisplPierCSMl{sample}(gm,:);
        
        %Find DCR (onedirectional)
        for ds=1:size(DSthreshL{sample,1},2)
            DCRmatrix{ds}=MatrixEDP./MatrixDS{ds};
            
            %Multidir CDR 1- Transverse
            temp=DCRmatrix{ds}(:,:,1);
            multidirDCRcsm{1,sample}(ds,gm)=max(max(temp));
            [CrMembTcsm{sample}{gm}(1,ds), CrMembTcsm{sample}{gm}(2,ds)]=find(temp==multidirDCRcsm{1,sample}(ds,gm),1);
            %2- Longitudinal
            temp=DCRmatrix{ds}(:,:,2);
            multidirDCRcsm{2,sample}(ds,gm)=max(max(temp));
            [CrMembLcsm{sample}{gm}(1,ds), CrMembLcsm{sample}{gm}(2,ds)]=find(temp==multidirDCRcsm{2,sample}(ds,gm),1);
            
        end
    end
 end

%% Fragility Analysis
x = 0:0.01:5; %IM interval

for cs=1:size(DSthreshL,1)
    for dr=1:2
        IM=IMcell(:,dr);
        for ds=1:4
            EDPcsm=multidirDCRcsm{dr,cs}(ds,:);
            [~, alphaCSMdcr{cs,dr}(ds),  betaCSMdcr{cs,dr}(ds), regrparCSMdcr{dr,ds}]=cloudFRAGls(IM,EDPcsm',1,'noplot');
        end
    end
end

%% Calculate fragility fractiles
x = 0:0.01:5; %IM interval
nds=4;

for dr=1:2
    temp=reshape([alphaCSMdcr{:,dr}],nds,size(DSthreshL,1));
    alphaMatrixCSM(:,:,dr)=temp'; clear temp
    temp=reshape([betaCSMdcr{:,dr}],nds,size(DSthreshL,1));
    betaMatrixCSM(:,:,dr)=temp'; clear temp      
    
    %%%% Fractiles on probabilities
    %CSM
    for ds=1:size(alphaMatrixCSM,2)
        for cs=1:size(alphaMatrixCSM,1)
            probabilitiesCSM(cs,:)=normcdf(log(x/alphaMatrixCSM(cs,ds,dr))/betaMatrixCSM(cs,ds,dr));
        end
        FragFractilesCSM{dr,ds}=prctile(probabilitiesCSM,[16 50 84]);
        FragilityBridge{dr,ds}=probabilitiesCSM;
    end   
end
clear probabilitiesCSM probabilitiesTH

% for dr=1:2   
%     %CSM
%     for ds=1:size(alphaMatrixCSM,2)
%         for cs=1:size(alphaMatrixCSM,1)
%             FragilityBridge{dr,ds}(cs,:)=normcdf(log(x/alphaMatrixCSM(cs,ds,dr))/betaMatrixCSM(cs,ds,dr));
%         end
%     end    
% end
%% Plot - Fragility curves

dsCol=[.5 .5 .5; .2 .2 .8; .8 .2 .2; .2 .8 .2];
correctColors=[1, .5];
lnStyle={':','-',':'};
lnWidth=[1 1.2 1.2 1.2];
FS=7;
xlabelStr='AvgSa [g]';
x = 0:0.01:5; %IM interval
xLim=1;
%dr=1;

figure1=figure('PaperPosition',[.25 .25 17 7]);
[a,posAxis]=tight_subplot(1,2,[0.05 0.075], [.1 .18],[.075 .075]);


drText={'Transverse dir.','Longitudinal dir.'};
imText={'$AvgSa(T_{T,min}-1.5T_{T,min})$','$AvgSa(T_{L}-1.5T_{L})$'};

titleDir={'Transv. dir.','Long. dir.'};

for dr=1:size(FragFractilesCSM,1)
    for ds=1:size(FragFractilesCSM,2)
        
        hold(a(dr),'on')
        
        % Boundaries CSM
        Xdispl = [x, fliplr(x)];
        inBetween = [FragFractilesCSM{dr,ds}(1,:),...
            fliplr(FragFractilesCSM{dr,ds}(3,:))];
        area=fill(Xdispl, inBetween,dsCol(ds,:)*1,'FaceAlpha',.2,'Parent',a(dr));
        set(area,'Edgecolor','w','Edgealpha',0)        
      
        % Median CSM
        plot(x,  FragFractilesCSM{dr,ds}(2,:),...
            lnStyle{2},'color',dsCol(ds,:),'Linewidth',lnWidth(1),'Parent',a(dr))
      
        xlim(a(dr),[0 xLim])
        ylim(a(dr),[0 1])
        set(a(dr),'FontSize',FS,...
            'Xtick',0:0.5:xLim,'Xticklabelmode','auto', 'YGrid','on','Ytick',[0:0.25:1],...
            'Yticklabelmode','auto','XGrid','off')
    end
    title(a(dr),titleDir{dr})

xlabel(a(dr),imText{dr},'Fontsize',8,'Interpreter','Latex');
ylabel(a(dr),'P [ DCR_{DS}>1 | IM ]','Fontsize',8)

end

dr=1;
DS1=plot(NaN,NaN,lnStyle{2},'color',dsCol(1,:),'Linewidth',lnWidth(1),'Displayname','DS1','Parent',a(dr));
DS2=plot(NaN,NaN,lnStyle{2},'color',dsCol(2,:),'Linewidth',lnWidth(1),'Displayname','DS2','Parent',a(dr));
DS3=plot(NaN,NaN,lnStyle{2},'color',dsCol(3,:),'Linewidth',lnWidth(1),'Displayname','DS3','Parent',a(dr));
DS4=plot(NaN,NaN,lnStyle{2},'color',dsCol(4,:),'Linewidth',lnWidth(1),'Displayname','DS4','Parent',a(dr));

l1=legend(a(dr),[DS1,DS2,DS3,DS4],'Location','southeast');
set(l1,'Position',[0.4 0.90 0.7 0.04],'Orientation','horizontal','NumColumns', 4,'Fontsize',8)
box(l1,'off')

print(figure1,fullfile(folderpath,'fragility_curves'),'-dpng','-r600');

%% Loss calculation
load(fileinputHazard)
for dr=1:2
    HazCurve=HazCurves{dr};

    for br=1:size(FragilityBridge{dr,1},1)
        for ds=1:4
            FragilityVect(ds,:)=FragilityBridge{dr,ds}(br,:);
            FragilityCurve{ds}=[x', FragilityBridge{dr,ds}(br,:)'];
        end
        LossCurveTemp=VulnCurveCalculator(DLRmatrix,  FragilityVect, x);
        LossRatioCurve{dr}(br,:)=LossCurveTemp(:,2);

        RiskCSM{dr}(br,:)=RiskCalculatorNum(HazCurve,FragilityCurve,'nolin','noplot');
        [EAL{dr}(br), lossExceedanceCurve{cs,dr}{br}]=EALcalculator(HazCurve,[LossCurveTemp(:,1),LossCurveTemp(:,2)],'noplot');
    end

    LossRatioCurvesFractiles{dr}=prctile(LossRatioCurve{dr},[16 50 84]);
    frctEAL(:,dr)=prctile(EAL{dr},[16 50 84]);
end

%% Plot - Fragility curves

vCol=[.2 .2 .2];
lnStyle={':','-',':'};
lnWidth=[1 1.2 1.2 1.2];
FS=7;
xlabelStr='AvgSa [g]';
x = 0:0.01:5; %IM interval
xLim=1.5;

figure1=figure('PaperPosition',[.25 .25 17 7]);
[a,posAxis]=tight_subplot(1,2,[0.1 0.1], [.15 .1],[.15 .1]);

drText={'Transverse dir.','Longitudinal dir.'};
imText={'$AvgSa(T_{T,min}-1.5T_{T,min})$','$AvgSa(T_{L}-1.5T_{L})$'};

titleDir={'Transv. dir.','Long. dir.'};

for dr=1:size(FragFractilesCSM,1)
    for ds=1:size(FragFractilesCSM,2)
        
        hold(a(dr),'on')
        
        % Boundaries 
        Xdispl = [x, fliplr(x)];
        inBetween = [LossRatioCurvesFractiles{dr}(1,:),...
            fliplr(LossRatioCurvesFractiles{dr}(3,:))];
        area=fill(Xdispl, inBetween,[.8 .8 .8],'FaceAlpha',.2,'Parent',a(dr));
        set(area,'Edgecolor','w','Edgealpha',0)    

             % Median CSM
        plot(x, LossRatioCurvesFractiles{dr}(2,:),...
            lnStyle{2},'color',vCol,'Linewidth',lnWidth(1),'Parent',a(dr))
        % Median CSM
        plot(x, LossRatioCurvesFractiles{dr}(1,:),...
            lnStyle{1},'color',vCol,'Linewidth',lnWidth(1),'Parent',a(dr))
        % Median CSM
        plot(x, LossRatioCurvesFractiles{dr}(3,:),...
            lnStyle{1},'color',vCol,'Linewidth',lnWidth(1),'Parent',a(dr))
        
         
        xlim(a(dr),[0 xLim])
        ylim(a(dr),[0 1])
        set(a(dr),'FontSize',FS,...
            'Xtick',0:0.5:xLim,'Xticklabelmode','auto', 'YGrid','on','Ytick',[0:0.25:1],...
            'Yticklabelmode','auto','XGrid','off')
    end
    title(a(dr),titleDir{dr})

xlabel(a(dr),imText{dr},'Fontsize',8,'Interpreter','Latex');
ylabel(a(dr),'Loss Ratio [ % C_{repl} ]','Fontsize',8)

end

VulnCurve=plot(NaN,NaN,lnStyle{2},'color',vCol,'Linewidth',lnWidth(1),'Displayname','Loss Curve (50th)','Parent',a(2));
VulnCurvePrc=plot(NaN,NaN,lnStyle{1},'color',vCol,'Linewidth',lnWidth(1),'Displayname','Loss Curve (16-84th)','Parent',a(2));

l=legend(a(2),[VulnCurve,VulnCurvePrc],'Location','southeast');
set(l,'Position',[-0.03 0.945 0.7 0.04],'Orientation','horizontal','NumColumns', 4,'Fontsize',8)
box(l,'off')

print(figure1,fullfile(folderpath,'loss_curves'),'-dpng','-r600');

%% save
save(fileoutput,'FragFractilesCSM','LossRatioCurvesFractiles',...
    'frctEAL','EAL','LossRatioCurve')

