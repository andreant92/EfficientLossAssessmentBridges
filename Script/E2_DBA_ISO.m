%% 2) DBA Assessment for simply supported bridges

% V1.0 Updated 25-04-2023
% Author: dr. Andrea Nettis. PhD
% Polytechnic University of Bari. Bari, Italy

% This script is aimed to perform the simplified performance assessment for
% simply-supported girder bridges (ISO)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%Settings
path=cd;
addpath(genpath('./Functions'))

bridgename='B2';

%Insert folder for input and output
folderpath=['./',bridgename]; 
fileoutput=[folderpath, '/SimplPerfAss'];
fileinput=[folderpath,'/SubassembliesInput.mat'];

fileSpectra='./Spectra.mat';

load(fileinput);
load(fileSpectra);

%% Analysis in the TRANSVERSE DIRECTION

dr=1; %Index for result allocation

%%%%%% Modification of consitutive laws
FDabutTemp=[0 0; 0.1 10^9; 0.2 2*10^9]; %Dummy fd law (rigid behavior)
PPt=cell(1,length(FDpier{dr}));
DSthreshT=cell(length(FDpier{dr}),2);

for cs=length(FDpier{dr}):-1:1 %1:length(FDpier{dr}) %length(FDpier{dr}):-1:1
    %%
    
    BearMatrix=[TableRED.abTypeT(cs,1), TableRED.bearTypeT(cs,:), TableRED.abTypeT(cs,2)];

    %%%%%%%% Initialize - Calculate F-Dlaws of the subassemblies
    FDbearEq=FDbear{dr}(cs,:); 
    FDpierEq=cell(1,length(FDbearEq));
    FDpierEq{1}=FDabutTemp; FDpierEq{end}=FDabutTemp;
    for pr=1:length(FDpier{dr}(cs,:))
        FDpierEq{pr+1}=FDpier{dr}{cs,pr};
    end

    tic
    disp(['Processing Realisation #',num2str(cs)])
    %%%%%%%%%%% 1 - find FDlaw for each pier+bearing subassembly
    for pr=1:length(FDbearEq)
        DSthreshT{cs,2}(pr,:)=[FDpierEq{pr}(2,1),FDpierEq{pr}(2,1)+0.5*(FDpierEq{pr}(3,1)-FDpierEq{pr}(2,1)), ...
            FDpierEq{pr}(2,1)+0.66*(FDpierEq{pr}(3,1)-FDpierEq{pr}(2,1)) FDpierEq{pr}(end,1)];
        %[FDpier{dr}{cs,pr}(2,1),FDpier{dr}{cs,pr}(end,1)];
        
        if BearMatrix(pr)==1 %General type bearing
            DSthreshT{cs,1}(pr,:)=[FDbearEq{pr}(2,1)*10, 10*(FDbearEq{pr}(2,1)+0.5*(FDbearEq{pr}(3,1)- FDbearEq{pr}(2,1))),...
                10*(FDbearEq{pr}(2,1)+0.66*(FDbearEq{pr}(3,1)- FDbearEq{pr}(2,1))), 10*FDbearEq{pr}(end,1)];
        elseif BearMatrix(pr)==2 %Neoprene
            DSthreshT{cs,1}(pr,:)=[FDbearEq{pr}(2,1), FDbearEq{pr}(2,1)+0.33*(.25- FDbearEq{pr}(2,1)),...
                .25, .40];
        end
    end

    %%%%%%%%%% 2 - Pushover of the series system pier+bear
    deckmassT=sum(TableRED.DeckMassesT{cs});
    piermass=[10 0.3*TableRED.pierMasses(cs,:)+TableRED.capbeamMasses(cs,:) 10]; %Dummy mass on abutments
    [deltaSDOF{cs}, shearSDOF{cs}, csiSDOF{cs}, ethaSDOF{cs}, massSDOF{cs},  teffSDOF{cs}, EDPds{cs}, crmemb{cs}, displbridge ]=DBAssTiterSpeed...
        (FDpierEq, FDbearEq,deckmassT,piermass,DSthreshT(cs,:),'noplot');

    displBridgeT{cs}=displbridge;

    %Re-arrangement
    displ=reshape([deltaSDOF{cs}{:}],[],length(deltaSDOF{cs}));
    shear=reshape([shearSDOF{cs}{:}],[],length(shearSDOF{cs}));
    mass=reshape([massSDOF{cs}{:}],[],length(massSDOF{cs}));
    csi=reshape([csiSDOF{cs}{:}],[],length(csiSDOF{cs}));
    
    %Capcurve calculation (IPM)
    for pr=1:length(deltaSDOF{cs})
        capcurve{pr}=[displ(:,pr) shear(:,pr)./mass(:,pr)];
        teffT(cs, pr)= 2*pi * (  capcurve{pr}(20,1)./( capcurve{pr}(20,2)) ).^0.5;
        deltaY(pr)=EDPds{cs}(pr,1);
    end     
    
    %%%%%%%%%% 3 - DCR calculation performing modCSM for each subassembly 
     [DCR{dr,cs}, PPt{cs}, crSUBt{cs}] =bridgeSSdcr(ADRSspectra,capcurve,deltaY,csi,EDPds{cs});
    toc
    
    %Extract results (transverse)
    for gm=1:length(PPt{cs})
        for pr=1:size(PPt{cs},1)
            PPsdof=PPt{cs}{pr,gm}(1);
            PerfDisplDeckCSMt{cs}(gm,pr)=interp1(deltaSDOF{cs}{pr}(2:end),displBridgeT{cs}{pr}(1,(2:end)),PPsdof);
            PerfDisplDeckCSMt{cs}(isnan(PerfDisplDeckCSMt{cs}))=0.0001; %Avoiding NaN for further elaborations
            PerfDisplPierCSMt{cs}(gm,pr)=interp1(deltaSDOF{cs}{pr}(2:end),displBridgeT{cs}{pr}(2,(2:end)),PPsdof);
            PerfDisplPierCSMt{cs}(isnan(PerfDisplPierCSMt{cs}))=0.0001; %Avoiding NaN for further elaborations
            PerfDisplBearCSMt{cs}(gm,pr)=PerfDisplDeckCSMt{cs}(gm,pr)-PerfDisplPierCSMt{cs}(gm,pr);
            PerfDisplBearCSMt{cs}(isnan(PerfDisplBearCSMt{cs}))=0.0001; %Avoiding NaN for further elaborations

        end
    end
        
end

clear capcurve deltaY

%% Analysis in the LONGITUDINAL DIRECTION

dr=2; 

PPl=cell(1,length(FDpier{dr}));

for cs=length(FDpier{dr}):-1:1
    %%
    tic
    disp(['Processing case ', num2str(cs)])
    %%%%%%%%%% 1 - Input elaboration
    BearMatrix=[TableRED.abTypeL(cs,1), TableRED.bearTypeL(cs,:), TableRED.abTypeL(cs,2)];

    FDabutTemp=[0 0; 0.1 10^9; 0.2 2*10^9];
    %Initialize and arrange force displ relations
    FDbearEq=FDbear{dr}(cs,:);
    FDpierEq=cell(1,length(FDbearEq));
    FDpierEq{1}=FDabutTemp; FDpierEq{end}=FDabutTemp;
    for pr=1:length(FDpier{dr}(cs,:))
        FDpierEq{pr+1}=FDpier{dr}{cs,pr};
    end
    
    %Define DS thresholds %%NO fix point (Modify it if necessary)
    for pr=1:length(FDbearEq)
        DSthreshL{cs,2}(pr,:)=[FDpierEq{pr}(2,1),FDpierEq{pr}(2,1)+0.5*(FDpierEq{pr}(3,1)-FDpierEq{pr}(2,1)),...
            FDpierEq{pr}(2,1)+0.66*(FDpierEq{pr}(3,1)-FDpierEq{pr}(2,1)), FDpierEq{pr}(end,1)];
        %[FDpier{dr}{cs,pr}(2,1),FDpier{dr}{cs,pr}(end,1)];
        %[FDpier{dr}{cs,pr}(2,1), FDpier{dr}{cs,pr}(2,1)+0.5*(FDpier{dr}{cs,pr}(3,1)-FDpier{dr}{cs,pr}(2,1)), FDpier{dr}{cs,pr}(2,1)+0.66*(FDpier{dr}{cs,pr}(end,1)- FDpier{dr}{cs,pr}(2,1)), FDpier{dr}{cs,pr}(end,1)];
        
        if BearMatrix(pr)==1
            DSthreshL{cs,1}(pr,:)=[10*FDbearEq{pr}(2,1), 10*(FDbearEq{pr}(2,1)+0.5*(FDbearEq{pr}(3,1)- FDbearEq{pr}(2,1))),...
                10*(FDbearEq{pr}(2,1)+0.66*(FDbearEq{pr}(3,1)- FDbearEq{pr}(2,1))), 10*FDbearEq{pr}(end,1)];
        elseif BearMatrix(pr)==2
            DSthreshL{cs,1}(pr,:)=[FDbearEq{pr}(2,1), FDbearEq{pr}(2,1)+0.33*(.25- FDbearEq{pr}(2,1)),...
                .25, .40];
            %[FDbear{dr}{cs,pr}(2,1),FDbear{dr}{cs,pr}(2,1)+0.5*(FDbear{dr}{cs,pr}(3,1)- FDbear{dr}{cs,pr}(2,1)),FDbear{dr}{cs,pr}(end,1), Inf];
        end
    end
    
    %Define fixed points (fixed or quasi-fixed bearings)
    fixpoint=(max(TableRED.fixL{cs}));

    %Define masses  %%NO fix point (Modify it if necessary)
%     temp=sum(TableRED.DeckMassesL{cs});
%     deckmassL=temp(fixpoint);
    deckmassL=sum(TableRED.DeckMassesT{cs});
    piermass=[10 0.3*TableRED.pierMasses(cs,:)+TableRED.capbeamMasses(cs,:) 10]; %Dummy mass on abutments
    
    AbutBackfillInter=1; 
    AbutBkfVect=zeros(1,length(FDbearEq));
    if AbutBackfillInter==1
        %Modify the fd law of the bearing of the first/last abutment to account for abut-backfill
        %interaction for push actions
        if BearMatrix(end)==2 %for neoprene bearing (put the abtbkf FDlaw in the 'end' pier
            
            FDabut{cs}(end,2)=FDabut{cs}(end,2)*1.01;
            
            AbutBkfVect(end)=1;
            displ=unique(sort([FDbearEq{end}(:,1)',FDabut{cs}(:,1)']));
            tempBear=interp1(FDbearEq{end}(:,1),FDbearEq{end}(:,2),displ,'linear','extrap');
            tempAbBkf=interp1(FDabut{cs}(:,1),FDabut{cs}(:,2),displ,'linear','extrap');
            FDbearEq{end}=[displ', tempBear'+tempAbBkf'];
            
            %plot(displ', tempBear'+tempAbBkf')
%             [temp,~,ic]=unique(FDbearEq{end}(:,2))
%             FDbearEq{end}(ic,2)
            
            DSthreshL{cs,1}(end,1)=min( DSthreshL{cs,1}(end,1),FDabut{cs}(2,1));
            DSthreshL{cs,1}(end,2)=min( DSthreshL{cs,1}(end,2),FDabut{cs}(3,1));
            DSthreshL{cs,1}(end,3)=min( DSthreshL{cs,1}(end,3),FDabut{cs}(4,1));
            DSthreshL{cs,1}(end,4)=5; %min( DSthreshL{cs,1}(end,3),FDabut{cs}(4,1));
                
%             figure
%             hold on
%             plot(displ,tempBear)
%             plot(displ,tempAbBkf)
%             
%             if FDbearEq{end}(2,1)<FDabut{cs}(2,1)
%                 FDbearEq{end}= [0 0; FDbearEq{end}(2,1) FDbearEq{end}(2,2), ...
%                     FDabut{cs}(2,1) max(interp1(FDbearEq{end}(:,1),FDbearEq{end}(:,2),FDabut{cs}(2,1)),0.01);... %to avoid unique points in interp1
%                     FDabut{cs}(3,1) interp1(FDbearEq{end}(:,1),FDbearEq{end}(:,2),FDabut{cs}(3,1))+ FDabut{cs}(3,2);...
%                     %FDabut{cs}(3,1) FDabut{cs}(3,2);...
%                     FDabut{cs}(4,1) interp1(FDbearEq{end}(:,1),FDbearEq{end}(:,2),FDabut{cs}(4,1))+ FDabut{cs}(4,2)];
%                 %Damage state modification
%                 DSthreshL{cs,1}(end,2)=min( DSthreshL{cs,1}(end,2),FDabut{cs}(3,1));
%                 DSthreshL{cs,1}(end,3)=min( DSthreshL{cs,1}(end,3),FDabut{cs}(4,1));
%                 
%             elseif FDbearEq{end}(2,1)>=FDabut{cs}(2,1)
%                 FDbearEq{end}= [0 0;  ...
%                     FDabut{cs}(2,1) max(interp1(FDbearEq{end}(:,1),FDbearEq{end}(:,2),FDabut{cs}(2,1)),0.01);...%to avoid unique points in interp1
%                     %FDabut{cs}(3,1) interp1(FDbearEq{end}(:,1),FDbearEq{end}(:,2),FDabut{cs}(3,1))+ FDabut{cs}(3,2)];
%                     FDabut{cs}(3,1) interp1(FDbearEq{end}(:,1),FDbearEq{end}(:,2),FDabut{cs}(3,1))+ FDabut{cs}(3,2);...
%                     FDabut{cs}(4,1) interp1(FDbearEq{end}(:,1),FDbearEq{end}(:,2),FDabut{cs}(4,1))+ FDabut{cs}(4,2)];
%                 %Damage state modification
%                 DSthreshL{cs,1}(end,1)=min( DSthreshL{cs,1}(end,2),FDabut{cs}(2,1));
%                 DSthreshL{cs,1}(end,2)=min( DSthreshL{cs,1}(end,2),FDabut{cs}(3,1));
%                 DSthreshL{cs,1}(end,3)=min( DSthreshL{cs,1}(end,2),FDabut{cs}(4,1));
%             end
            
        elseif BearMatrix(end)==1 && BearMatrix(1)==1
            %Ab-back interact correction is provided for the movable
            %abutment, dummyfixpoint is made up of ones on the bent and 1
            %or 0 on the abutments to find the movable abutment obtaining
            %the right index to correct FDbear and DSthresh
            dummyfixpoint=ones(1,length(fixpoint));dummyfixpoint(1)=fixpoint(1);dummyfixpoint(end)=fixpoint(end);
            ind=find(dummyfixpoint==0);
            AbutBkfVect(ind)=1;
            FDbearEq{ind}= [0 0;  ...
                FDabut{cs}(2,1) max(interp1(FDbearEq{ind}(:,1),FDbearEq{ind}(:,2),FDabut{cs}(2,1)),0.01);...%to avoid unique points in interp1
                %FDabut{cs}(3,1) interp1(FDbearEq{ind}(:,1),FDbearEq{ind}(:,2),FDabut{cs}(3,1))+ FDabut{cs}(3,2)];
                FDabut{cs}(3,1) FDabut{cs}(3,2);...
                FDabut{cs}(4,1) FDabut{cs}(4,2)];
                
            %Damage state modification
            DSthreshL{cs,1}(ind,1)=FDabut{cs}(2,1);
            DSthreshL{cs,1}(ind,2)=FDabut{cs}(3,1);
            DSthreshL{cs,1}(ind,3)=FDabut{cs}(4,1);
        end
    end

    %%%%% Check plot
    CheckPlotter=0;
    if CheckPlotter==1
     figure
     subplot(1,2,1)
     hold on
     for i=1:6
         plot(FDbearEq{i}(:,1),FDbearEq{i}(:,2),'displayname',['sub',num2str(i)])
     end
     legend
     subplot(1,2,2)
     hold on
     for i=1:6
         plot(FDpierEq{i}(:,1),FDpierEq{i}(:,2),'displayname',['sub',num2str(i)])
     end
     legend
     axis([0 .2 0 20000])
    end
    
    %Exclude fixed points and correct it to evaluate loss of support
    goodPoint=1:length(FDbearEq); suppLossPoint=zeros(1,length(FDbearEq));
    if BearMatrix(end)==1 && fixpoint(end)~=0
        suppLossPoint(end)=1;
%         FDbearEq{end}=[0 0; 0.2 0.005; 0.8 0.02]; %Free abutment (no friction)
        DSthreshL{cs,1}(end,:)=[Inf Inf Inf 0.4]; %LOss of support a 0.6m
        goodPoint(end)=[];
    elseif BearMatrix(1)==1 && fixpoint(1)~=0
        suppLossPoint(1)=1;
        goodPoint(1)=[];
%         FDbearEq{1}=[0 0; 0.2 0.005; 0.8 0.02]; %Free abutment (no friction)
        DSthreshL{cs,1}(1,:)=[Inf Inf Inf 0.4]; %LOss of support a 0.6m
    end

    %Modify DS matrix to exclude fixed abutments
    DSthreshRED{1}=DSthreshL{cs,1}(goodPoint,:);
    DSthreshRED{2}=DSthreshL{cs,2}(goodPoint,:);
    
    %%%%%%%%%% 2 - Capacity curve attainment
    [deltaSDOFbr{cs}, VbaseSDOFbr{cs}, massSDOFbr{cs}, kSDOFbr{cs}, tSDOFbr{cs}, ethaSDOFbr{cs}, csiSDOFbr{cs}, edpDSbr{cs}, CrMembBr{cs}, displBridgeL{cs}] = ...
        DBAssL(FDpierEq(goodPoint), FDbearEq(goodPoint), deckmassL(goodPoint), piermass(goodPoint), DSthreshRED, 'noplot',AbutBkfVect(goodPoint) );

    %Capcurve calculation (IPM)
    capcurve=[deltaSDOFbr{cs} VbaseSDOFbr{cs}./massSDOFbr{cs}];
    teffL(cs, 1)= 2*pi * (  capcurve(5,1)./( capcurve(5,2)) ).^0.5;
    teffL(cs, 2)= 2*pi * (  capcurve(end-3,1)./( capcurve(end-3,2)) ).^0.5;
    
    deltaY=edpDSbr{cs}(1);
    %plot(capcurve(:,1),capcurve(:,2))


%%%%%%%%%% 3 - DCR calculation performing modCSM for the SDOF model of the bridge
    [DCR{dr,cs},PPl{cs},crSUBl{cs}]=bridgeSSdcr(ADRSspectra,{capcurve},deltaY,csiSDOFbr{cs}, edpDSbr{cs});
    toc
    
    %Accounting for the excluded points
    indexes=find(suppLossPoint);
    temp=[displBridgeL{cs}{1}(1,:); zeros(1,length(displBridgeL{cs}{1}))];
    for i=1:length(indexes)
        temp2=[displBridgeL{cs}(1:indexes(i)-1), temp, displBridgeL{cs}(indexes(i):end)];
        displBridgeL{cs}=temp2;
    end
    
    %Extract results (longitudinal)    
     for gm=1:length(PPt{cs})
        PPsdof=PPl{cs}{gm}(1);
        for pr=1:size(displBridgeL{cs},2)
            PerfDisplDeckCSMl{cs}(gm,pr)=interp1(deltaSDOFbr{cs}(2:end),displBridgeL{cs}{pr}(1,(2:end)),PPsdof,'linear','extrap');
            PerfDisplPierCSMl{cs}(gm,pr)=interp1(deltaSDOFbr{cs}(2:end),displBridgeL{cs}{pr}(2,(2:end)),PPsdof,'linear','extrap');
            PerfDisplBearCSMl{cs}(gm,pr)=PerfDisplDeckCSMl{cs}(gm,pr)-PerfDisplPierCSMl{cs}(gm,pr);
        end
    end
    
end

clear capcurve deltaY AbutBackfillInter CheckPlotter deckmassL piermass FDabutTemp FDbearEq FDpierEq cs dr pr

% %% Global DCR evaluation: CSM
% 
% for cs=size(TableRED,1):-1:1   
%     %% Find capacity (Matrix with damage states thresholds)
%     for ds=1:size(DSthreshL{cs,1},2)
%         MatrixDS{ds}=zeros(2,size(DSthreshL{cs,1},1),2);
%         for pr=1:size(DSthreshL{cs,2},1)
%             MatrixDS{ds}(1,pr,1)=DSthreshT{cs,1}(pr,ds);        %bear transv
%             MatrixDS{ds}(2,pr,1)=DSthreshT{cs,2}(pr,ds);        %pier transv
%             MatrixDS{ds}(1,pr,2)=DSthreshL{cs,1}(pr,ds);        %bear long
%             MatrixDS{ds}(2,pr,2)=DSthreshL{cs,2}(pr,ds);        %pier long
%         end
%     end
%     %% Find demand and DCR
%     
%     for gm=1:length(ADRSspectra)
%         %%
%         % 3d matrix having [dim1 (row1): bear deform, row2: pier deform
%         % dim2: number of subassemblies, dim3: direction
%         
%         MatrixEDP=zeros(2,size(DSthreshT{cs,1},1),2);
%         % Find peformance edp (displacement) in trasnverse direction
%         MatrixEDP(1,:,1)=(PerfDisplBearCSMt{cs}(gm,:));
%         MatrixEDP(2,:,1)=PerfDisplPierCSMt{cs}(gm,:);
%                 
%         % Find performance edp (displacement) in longitudinal direction        
%         MatrixEDP(1,:,2)=(PerfDisplBearCSMl{cs}(gm,:));
%         MatrixEDP(2,:,2)=PerfDisplPierCSMl{cs}(gm,:);
%       
%         %Find DCR (multidirectional)
%         for ds=1:size(DSthreshL{cs,1},2)
%             DCRmatrix{ds}=MatrixEDP./MatrixDS{ds};
%             
%             %Multidir CDR 1- Transverse 
%             temp=DCRmatrix{ds}(:,:,1);
%             multidirDCRcsm{1,cs}(ds,gm)=max(max(temp));
%             [CrMembTcsm{cs}{gm}(1,ds), CrMembTcsm{cs}{gm}(2,ds)]=find(temp==multidirDCRcsm{1,cs}(ds,gm),1);
%             %2- Longitudinal
%             temp=DCRmatrix{ds}(:,:,2);
%             multidirDCRcsm{2,cs}(ds,gm)=max(max(temp));
%             [CrMembLcsm{cs}{gm}(1,ds), CrMembLcsm{cs}{gm}(2,ds)]=find(temp==multidirDCRcsm{2,cs}(ds,gm),1);
% 
%         end
%         
%         %Find DCR (onedirectional)
%         for ds=1:size(DSthreshL{cs},2)
%             DCRmatrix{ds}=MatrixEDP./MatrixDS{ds};
%             %bearings, square shaped - strength interaction domain
%             onedirDCRmatrix{ds}(1,:)= max(DCRmatrix{ds}(1,:,:),[],3);
%             %piers, elliiptic shaped - strength interaction domain
%             onedirDCRmatrix{ds}(2,:)= (sum(DCRmatrix{ds}(2,:,:).^2,3)).^0.5;
%             [~,critdir{cs}(gm,ds)]=max(max(max(DCRmatrix{ds},[],2),[],1));
%             onedirDCR{cs}(gm,ds)=max(onedirDCRmatrix{ds}(:));
%         end
%     end
% end

%% Plot
%%%% control plot, all the capacity curves are shown for transverse and
%%%% longitudinal directions

%transverse direction
figure
hold on
col=linspace(0.2, 0.8,length(deltaSDOF{1}))'*ones(1,3); 
for cs=1:length(deltaSDOF)
    for pr=1:length(deltaSDOF{cs})
%         taglicheck{pr}(cs,:)=shearSDOF{cs}{pr};
%         [Maxsh(pr),csmax(pr)]=max(taglicheck{pr}(:,end))
        plot(deltaSDOF{cs}{pr}, shearSDOF{cs}{pr},'color', col(pr,:))
        axis([0 0.6 0 8000])
    end
end
xlabel('Displacement [m]'); ylabel('Shear [kN]')
% P1=plot(NaN,NaN,'color', col(2,:),'displayname','Sub. 1'); P2=plot(NaN,NaN,'color', col(3,:),'displayname','Sub. 2');
% legend([P1,P2])
title('Transv direction: Variation in FD law')
print(gcf,fullfile(folderpath,'capcur_transvdir'),'-dpng','-r600');

%longitudinal direction
figure
hold on
col=[.5 0 0];
for cs=1:length(deltaSDOF)
    plot(deltaSDOFbr{cs}, VbaseSDOFbr{cs},'color', col)
    axis([0 0.6 0 15000])
end
xlabel('Displacement [m]'); ylabel('Shear [kN]')
title('Longitudinal direction: Variation in FD laws')
print(gcf,fullfile(folderpath,'capcur_longdir'),'-dpng','-r600');


%% salva
DBAres{1}=table(deltaSDOF', shearSDOF', csiSDOF', ethaSDOF', massSDOF', teffSDOF', EDPds', crmemb',displBridgeT','VariableNames',...
    {'deltaSDOF', 'shearSDOF', 'csiSDOF', 'ethaSDOF', 'massSDOF', 'teffSDOF', 'EDPds', 'crMemb', 'displBridge'});
DBAres{2}=table(deltaSDOFbr', VbaseSDOFbr', csiSDOFbr', ethaSDOFbr', massSDOFbr',tSDOFbr', edpDSbr',displBridgeL','VariableNames',...
    {'deltaSDOF', 'shearSDOF', 'csiSDOF', 'ethaSDOF', 'massSDOF', 'teffSDOF', 'EDPds', 'displBridge'});

% CSMres=table(multidirDCRcsm(1,:)',multidirDCRcsm(2,:)',onedirDCR',CrMembTcsm',CrMembLcsm','VariableNames',...
%     {'DCRt','DCRl','onedirDCR','CrMembT','CrMembL'});
% 
% CSMres=table(multidirDCRcsm(1,:)',multidirDCRcsm(2,:)',onedirDCR',CrMembTcsm',CrMembLcsm','VariableNames',...
%     {'DCRt','DCRl','onedirDCR','CrMembT','CrMembL'});

salva=1;

if salva == 1
    %save(fileoutput)
    %save(fileoutput,'CSMres','DBAres')
    %save(fileoutput,'DBAres')
    save(fileoutput,'PerfDisplBearCSMl', 'PerfDisplBearCSMt', 'PerfDisplDeckCSMl', 'PerfDisplDeckCSMt', 'PerfDisplPierCSMl', 'PerfDisplPierCSMt',...
    'DSthreshT','DSthreshL','DBAres')
else
    disp ('no save active!!')
end
