%% 1) Input data and LHS

% V1.0 Updated 25-04-2023
% Author: dr. Andrea Nettis. PhD
% Polytechnic University of Bari. Bari, Italy

% This script is aimed to import the data form and to perform the LHS to
% generate bridge model realisations

%% %%%%%%%%%%%%%%%%%%%%% INPUT 
clc
clear
close all

%Settings
path=cd;
addpath(genpath('./Functions'))

bridgename='B2'; %Insert an id code for the bridge
filenameXLSX=['./xlsx/',bridgename,'.xlsx']; %Insert the path of the input xlsx
folderpathoutput=['./',bridgename]; %Insert output folder
mkdir(folderpathoutput)

%Additional input
n=10; %Number of samples from LHS
corrMatClass=1; %Correlate material properties with mat classes?

%% ========================== RUN 
%% Import the deterministic data
opts = spreadsheetImportOptions("NumVariables", 1);

%================= 1st IMPORT
% Specify sheet and range
opts.Sheet = "MAIN";
opts.DataRange = "C19:C19";
% Specify column names and types
opts.VariableNames = "VarName3";
opts.VariableTypes = "string";
% Specify variable properties
opts = setvaropts(opts, "VarName3", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "VarName3", "EmptyFieldRule", "auto");

% Import the data
SS = readtable(filenameXLSX, opts, "UseExcel", false);
% Convert to output type
SS = table2cell(SS);

% Clear temporary variables
clear opts

%================= 2nd IMPORT
opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "ISO";
opts.DataRange = "C5:G60";

% Specify column names and types
opts.VariableNames = ["VarName", "VarType", "VarValue", "detProb", "probNumb"];
opts.VariableTypes = ["char", "char", "char", "categorical", "double"];

% Specify variable properties
opts = setvaropts(opts, ["VarName", "VarType", "VarValue"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName", "VarType", "VarValue", "detProb"], "EmptyFieldRule", "auto");

% Import the data
DataTable = readtable(filenameXLSX, opts, "UseExcel", false); emptyrow=[];
for row=1:length(DataTable.VarName)
    if isempty(DataTable.VarName{row})
        emptyrow=[emptyrow, row];
    end
end
DataTable(emptyrow,:)=[];   
% Clear temporary variables
clear opts numIdx

% ======================== Deterministic Data arrangement
DataDET = DataTable(DataTable.detProb == 'DET',1:3);
% Convert to output type
Data = table2cell(DataDET);
numIdx = cellfun(@(x) ~isnan(str2double(x)), Data);
Data(numIdx) = cellfun(@(x) {str2double(x)}, Data(numIdx));

%Create variables in the Workspace 
for i=1:length(Data)
    if DataDET{i,2}=="string"
        eval(strcat(Data{i,1},"='",Data{i,3},"'"));
%         eval(strcat(Data{i,1},"tab='",Data{i,3},"'"))
    elseif Data{i,2}=="double"
        t=convertStringsToChars(Data{i,1}); t1=num2str(Data{i,3});
        t1 = strrep(t1, ',', '.');
        eval([t,'=[',t1,']']);
    end
end

clear i t1 t DataDET connectiontab connection

%% Import probabilistic data and perform LHS

% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 8);

% Specify sheet and range
opts.Sheet = "PD";
opts.DataRange = "A3:H15";

% Specify column names and types
opts.VariableNames = ["ProbIndex","Distribution", "MeanOrLower", "StdOrUpper", "Min", "Max", "Discr", "Trasp"];
opts.VariableTypes = ["double","char", "double", "double", "double", "double", "char", "double"];

% Specify variable properties
opts = setvaropts(opts, ["Distribution", "Discr"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Distribution", "Discr"], "EmptyFieldRule", "auto");

% Import the data
PDtable = readtable(filenameXLSX, opts, "UseExcel", false);

% Clear temporary variables
clear opts numIdx raw

%Delete empty lines
PDtable(isnan(PDtable.ProbIndex),:)=[];

%=========== Definition of INPUT probability distribution

minabs=[]; maxabs=[];

pdMAT=cell(1,size(PDtable,1));
for i=1:size(PDtable,1)
    
    %steel tensile strenght [MPa] %FeB32k by Verderame (2001)
    %concrete strenght [Mpa] from Zelaschi (2016)
    %DECK concrete strenght [Mpa] from Zelaschi (2016)
    
    pdMAT{i} = makedist(PDtable.Distribution{i},PDtable.MeanOrLower(i),PDtable.StdOrUpper(i));
    minabs=[minabs PDtable.Min(i)];
    maxabs=[maxabs PDtable.Max(i)];
    
%     % UNLOCK for a check-plot
%     figure
%     hold on
%     subplot(1,2,1)
%     x=1:10;
%     plot(x,pdf(pdMAT{i},x))
%     subplot(1,2,2)
%     plot(x,cdf(pdMAT{i},x))
end

for i=1:length(pdMAT)
    minNORM(i)=cdf(pdMAT{i},minabs(i));
    maxNORM(i)=cdf(pdMAT{i},maxabs(i));
end

%% Arrange INPUT data for bridge samples generation

% Sampling LH
%Generation of latin hypercube between min e max
x=lhsdesign_modified(n, minNORM, maxNORM); % generate latin hypercube samples (MATLAB Function, see MATLAB documentation for more information)

samplesLHS=zeros(n,length(pdMAT));                                             % preallocation for the matrix
for i=1:length(pdMAT)
    %%LH sampling
    prob=x(:,i);
    probtype=pdMAT{i}.DistributionName;
    samplesLHS(:,i) = icdf(pdMAT{i},prob);                   % map latin hypercube samples to values using inverse cumulative distribution functions
    
    %Convert Discrete random variables
    if length(str2num(PDtable.Discr{i}))>1
        %Discretization
        values=str2num(PDtable.Discr{i});
        samplesLHS(:,i)=round(samplesLHS(:,i));
        samplesLHS(:,i)=interp1(unique(samplesLHS(:,i)),values,samplesLHS(:,i));
    end
%     figure
%     hist(samplesLHS(:,i),50)
end

%Probabilistic Matrix Data
DataProb=DataTable(DataTable.detProb=='PD',:);

% Import Spreadsheet with probabilistic variable transposition
opts = spreadsheetImportOptions("NumVariables", 7);

% Specify sheet and range
opts.Sheet = "TRASP";
opts.DataRange = "B3:H20";
% Specify column names and types
opts.VariableNames = ["PDindex", "Name", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8"];
opts.VariableTypes = ["double", "char", "char", "char", "char", "char", "char"];

% Specify variable properties
opts = setvaropts(opts, ["Name", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Name", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8"], "EmptyFieldRule", "auto");

% Import the data
TRtable = readtable(filenameXLSX, opts, "UseExcel", false);
TRtable(isnan(TRtable.PDindex),:)=[];

%Trasposition if necessary (e.g. transverse stiffness configuration
for i=1:size(PDtable,1)
   if PDtable.Trasp(i)==0
       index=DataProb.probNumb==PDtable.ProbIndex(i);
       eval([DataProb.VarName{index},'=samplesLHS(:,i)']);
   elseif PDtable.Trasp(i)==1
        indprob=PDtable.ProbIndex(i);    
        rows=find(TRtable.PDindex==indprob);
        maxind=max(str2num(PDtable.Discr{i}));
        CellTemp = table2cell(TRtable(rows(:),3:maxind+2));
        TRdata=cell(length(samplesLHS),size(CellTemp,1));
        for cs=1:length(samplesLHS)
            TRdata(cs,:)=CellTemp(:,samplesLHS(cs,indprob));
        end

        for r=1:length(rows)
            temp=TRdata(:,r);
            %Create the correspinding var
            k=strcmp(TRtable.Name{rows(r)},DataProb.VarName);
            if strcmp(DataProb.VarType(k),'double')
                for l=1:length(temp)
                    tempNum(l)=str2num(cell2mat(temp(l)));
                end
                eval([TRtable.Name{rows(r)},'=tempNum']);
                clear temp
            else
                eval([TRtable.Name{rows(r)},'=temp']);
                clear temp
            end              
        end
        clear rows TRdata %TRtable tempnum
   end    
end
    
clear rows
close all

%% Correlation between design class and mechanical properties
%Activate only if both material design class and mechanical properties are
%simulated 

if corrMatClass==1
    %Correlate the value obtained in LHS Sampling with the mean
    %corresponding to the design class
    
    %Find vector of fc fy sigmac and sigmas
    colfc=DataProb.probNumb(strcmp(DataProb.VarName,'fc'));
    colsigmac=DataProb.probNumb(strcmp(DataProb.VarName,'sigmac'));
    colfy=DataProb.probNumb(strcmp(DataProb.VarName,'fy'));
    colsigmas=DataProb.probNumb(strcmp(DataProb.VarName,'sigmas'));
    
    %find the avg values assuming that sigmac and sigmas are the
    %characteristic values (5%)
    sigmaCmean= samplesLHS(:,colsigmac)./(1/PDtable.StdOrUpper(colfc)-1.64).*(1/PDtable.StdOrUpper(colfc));
    sigmaSmean= samplesLHS(:,colsigmas)./(1/PDtable.StdOrUpper(colfy)-1.64).*(1/PDtable.StdOrUpper(colfy));
    
    %Scale the sampled value to the average
    fc=samplesLHS(:,colfc).*sigmaCmean;
    fy=samplesLHS(:,colfy).*sigmaSmean;
end

clear colfc colfy colsigmac colsigmas sigmaCmean sigmaSmean corrMatClass

%% Elaborations on superstructure-substructure connections

%%%%%%%%% Fixity definition on abutments
%Transposition in n size
[abx_fixity,aby_fixity,abx_type,aby_type]=checkOrdVars(n,abx_fixity,aby_fixity,abx_type,aby_type);

%Xdir (Longitudinal)
for a=1:length(abx_fixity)
    abxfix(a,:)=str2num(strrep(strrep(strrep(abx_fixity{a}, 'NE', '1'),'FX', '1'),'FR', '0'));
    abxtype(a,:)=str2num(strrep(strrep(strrep(abx_type{a}, 'N', '2'),'X', '1'),'S', '0'));
end

%Ydir (Transverse)
for a=1:length(aby_fixity)
    abyfix(a,:)=str2num(strrep(strrep(aby_fixity{a}, 'FX', '1'),'FR', '0')); %abut trasv
    abytype(a,:)=str2num(strrep(strrep(strrep(aby_type{a}, 'N', '2'),'X', '1'),'S', '0'));
end

%%%%%%%%% Fixity definition on piers
[bearx_fixity, beary_fixity,bearx_type,beary_type ]=checkOrdVars(n,bearx_fixity,beary_fixity,bearx_type,beary_type);

%Xdir (Longitudinal)
for  a=1:length(bearx_fixity)
    if strcmp(SS,'HYP')
        bearFixLong(:,a)=str2num(strrep(strrep(bearx_fixity{a}, 'FX', '1'),'FR', '0')); %bear long
    elseif strcmp(SS,'ISO')
        temp=str2num(strrep(strrep(strrep(strrep(bearx_fixity{a}, 'XR', '1'),... %bear long
            'RX', '2'),...
            'XX', '3'),...
            'RR', '0'));
        t1=[0 1 2 3]; t2=[0 0; 1 0; 0 1; 1 1]';
        for i=1:length(H_piers)
            ind=t1==temp(i);
            bearFixLong(:,i,a)=t2(:,ind);
        end
    end
    bearTypeLong(:,a)=str2num(strrep(strrep(strrep(bearx_type{a}, 'N', '2'),'X', '1'),'S', '0')); %bear long
end

%Y dir (Transverse)
for  a=1:length(beary_fixity)
    if strcmp(SS,'HYP')
        bearFixTransv(:,a)=str2num(strrep(strrep(beary_fixity{a}, 'FX', '1'),'FR', '0')); %bear long
    elseif strcmp(SS,'ISO')
        temp=str2num(strrep(strrep(strrep(strrep(bearx_fixity{a}, 'XR', '1'),... %bear long
            'RX', '2'),...
            'XX', '3'),...
            'RR', '0'));
        t1=[0 1 2 3]; t2=[0 0; 1 0; 0 1; 1 1]';
        for i=1:length(H_piers)
            ind=t1==temp(i);
            bearFixTransv(:,i,a)=t2(:,ind);
        end
    end
        bearTypeTransv(:,a)=str2num(strrep(strrep(strrep(beary_type{a}, 'N', '2'),'X', '1'),'S', '0')); %bear long
end

%% Vars elaboration for table generation
  
% Seismic DATA
for  a=1:length(year)
    if year(a)<1975
        seismDATA(a)=cat_seism(a);
    elseif year(a)>=1975
        seismDATA(a)=coeff_seism(a);
    end
end
clear cat_seism coeff_seism

%Capbeam
for a=1:size(Dim_capbeam,1)
    if length(Dim_capbeam(a,:))==3
        Hcapbeam(a)=Dim_capbeam(a,3);
    elseif length(Dim_capbeam(a,:))==5
        Hcapbeam(a)=Dim_capbeam(a,3)+Dim_capbeam(a,4);
    end
end

%Period of construction
for i=1:length(year)
    %Year elaboration
    if year(i) >= 1957 && year(i)<1961
        buildPer(i)=1;
    elseif year(i) >= 1962 && year(i)<1974
        buildPer(i)=2;
    elseif year(i) >= 1975 && year(i)<1979
        buildPer(i)=3;
    elseif year(i) >= 1980 && year(i)<1985
        buildPer(i)=4;
    elseif year(i) >= 1985 && year(i)<1989
        buildPer(i)=5;
    elseif year(i)>=1990 && year(i)<2005
        buildPer(i)=6;
    end
end

%Period of construction

if length(buildPer)==1
    buildPer=buildPer*ones(1,n);
end

roadDim=(horzcat(W_deck, W_walk, W_road)); %Dimensions of the deck
Hdeck=(horzcat(H_girder, H_slab)); %Height of the slab and girder
DesMat=(horzcat(sigmac, sigmas)); %Design material mechanical proprerties
G1vect=(areagird.*numgird+H_slab.*W_deck)*25; %[kN/m]%G1 deck
%Pier dimensions
if D_pier>0
    pierDim=D_pier;
else
    pierDim(:,:,1)=Lx_pier;
    pierDim(:,:,2)=Ly_pier;
end

% StructType=checkOrdVars(n,SS);
[year,buildPer,seismData, roadDim,Hdeck,catTraff, DesMat, lengthSpans, G2vect, G1vect,pierDim, Hpier, pierType,distBear,dCapbeam,Hbear,StructType,NumGird] = checkOrdVars(n,...
    year,buildPer,seismDATA,roadDim,Hdeck,cat_traffic, DesMat,L_spans, G2_deck,G1vect,pierDim,H_piers,pier_type,d_bear,Dim_capbeam,H_bear,SS,numgird);
[bearHeight] = checkOrdVars(n,bear_height);
[bearStiffFixed,bearStiffNeo, gapAbut, heightBkw, widthBkw, kBkw] = checkOrdVars(n,bearstiff_FIX,bearstiff_NEO, gap_ab, height_bkw, width_bkw, k_bkw);

clear W_deck W_walk W_road DataProb DataTable PDTable CellTemp D_Pier areagird cat_traffic S seismDATA cat_traffic L_spans G2_deck Dim_capbeam ind index indprob k L_bridge maxind...
 max_abs Hcapbeam opts pdMAT pier_typeTab t1 t2 TRtable values x probtype gap_ab height_bkw width_bkw bearstiff_NEO bearstiff_FIX

%% Tributary length definition

% tributary length for vertical load
tribl=cell(n,1);
triblTOT=cell(n,1);
for a=1:length(samplesLHS)
    triblTOT{a}=horzcat([0; lengthSpans(a,1)/2], zeros(2,length(H_piers)),[lengthSpans(a,end)/2; 0]);
    for i=1:length(H_piers)
        triblTOT{a}(:,i+1)=[lengthSpans(a,i)/2;lengthSpans(a,i+1)/2]; %1st row tributary length on the left of the ab/pier
        %2nd row tributary length on the right of the ab/pier
    end
    tribl{a}=triblTOT{a}(:,2:end-1);
end

% tributary length for horizontal load
tribhor=cell(n,1);
tribhorTOT=cell(n,1);
for a=1:length(samplesLHS)
    if strcmp(SS,'HYP')
        %Need review (include free abutments in transv direction)
        fixity=[abfix(1,1),ones(1,length(bearFixLong)) , abfix(1,2);...
            abfix(2,1), bearFixLong, abfix(2,2)];
        ind=find(fixity(2,:)==1);
        tribhorTOT=zeros(1,length(fixity));
        for i=1:length(fixity)
            if fixity(2,i)==1
                tribhorTOT(i)= L_bridge/length(ind); %HP: infinite axial stiffness of the deck
            else
                tribhorTOT(i)=0;
            end
        end
        tribhor=tribhorTOT(2:end-1);
        warning('The routine for risk assessment for continuous-deck (HYP) structures is under development!')
        
    elseif strcmp(SS,'ISO')
        
        fixitylong=horzcat([NaN; abxfix(a,1)], bearFixLong(:,:,a), [abxfix(a,2); NaN]);
        tribhor{a}=zeros(2,length(H_piers));
        %associo per ogni pila la tributary length alle azioni orizzontali
        for i=2:length(fixitylong)-1
            if fixitylong(1,i)==1 && fixitylong(2,i-1)==1
                tribhor{a}(1,i-1)=lengthSpans(a,i-1)/2;
            elseif fixitylong(1,i)==1 && fixitylong(2,i-1)==0
                tribhor{a}(1,i-1)=lengthSpans(a,i-1);
            end
            if fixitylong(2,i)==1 && fixitylong(1,i+1)==1
                tribhor{a}(2,i-1)=lengthSpans(a,i)/2;
            elseif fixitylong(2,i)==1 && fixitylong(1,i+1)==0
                tribhor{a}(2,i-1)=lengthSpans(a,i);
            end
        end
        %tribhorTOT contiene anche tribhor delle spalle
        tribhorTOT{a}=zeros(2,length(fixitylong));
        tribhorTOT{a}(:,2:end-1)= tribhor{a};
        %spalla sx
        if fixitylong(2,1)==1 && fixitylong(1,2)==1
            tribhorTOT{a}(2,1)=lengthSpans(a,1)/2;
        elseif fixitylong(2,1)==1 && fixitylong(1,2)==0
            tribhorTOT{a}(2,1)=lengthSpans(a,1);
        end
        %spalla dx
        if fixitylong(1,end)==1 && fixitylong(2,end-1)==1
            tribhorTOT{a}(1,end)=lengthSpans(a,end)/2;
        elseif fixitylong(1,end)==1 && fixitylong(2,end-1)==0
            tribhorTOT{a}(1,end)=lengthSpans(a,end);
        end
        
    end
end

%% Pier - specific data allocation

PierDATA=cell(n,1);

for a=1:n
    PierDATA{a}=struct([]);
    for i=1:length(Hpier(a,:))
        PierDATA{a}(i).Hbent=Hpier(a,i);
        if strcmp(pierType{a}, 'RECT')
            Apier=pierDim(a,i)*pierDim(a,i); %[mq]
            PierDATA{a}(i).sect(1)=pierDim(a,i,1);
            PierDATA{a}(i).sect(2)=pierDim(a,i,2);
        elseif strcmp(pierType{a}, 'CIRC')
            Apier=pi*(pierDim(a,i)/2)^2; %[mq]
            PierDATA{a}(i).sect=pierDim(a,i);
        end
        
        PierDATA{a}(i).PPbent=(Hpier(a,i)-(dCapbeam(a,3)+dCapbeam(a,4)))*Apier*25; %[kN]
        PierDATA{a}(i).tribl=tribl{a}(:,i); %[m]
        PierDATA{a}(i).tribhor=tribhor{a}(:,i); %[m]
        
        if strcmp(SS,'HYP')
            PierDATA{a}(i).d_bear=0; %[m]
        elseif strcmp(SS,'ISO')
            PierDATA{a}(i).d_bear=distBear(a,:);
        end
        
        PierDATA{a}(i).PPcapbeam=((dCapbeam(a,1)*dCapbeam(a,3))+...
            ((dCapbeam(a,1)+dCapbeam(a,2))*dCapbeam(a,4)/2))*dCapbeam(a,end)*25;%[kN]
        
        PierDATA{a}(i).PPdeck=[G1vect(a)+G2vect(a), 0]; %[kN/m] [weigth, trasv dir eccentricity]
    end
end

%% RUN load analysis

for a=1:n
    tic
    for i=1:length(PierDATA{a})
        [PierDATA{a}(i).comb,PierDATA{a}(i).soll]=loadanalysisbridge(buildPer(a),roadDim(a,1),roadDim(a,2),roadDim(a,3),Hdeck(a,1),PierDATA{a}(i).Hbent,Hdeck(a,2),...
            PierDATA{a}(i).PPbent,PierDATA{a}(i).PPcapbeam,PierDATA{a}(i).tribl,PierDATA{a}(i).tribhor,PierDATA{a}(i).PPdeck,...
            catTraff(a),seismData(a),PierDATA{a}(i).d_bear,DesMat(a,1),pierType{a},PierDATA{a}(i).sect);
    end
    toc
end

%% RUN Simulated design
% if SD==1
    for a=1:n
        tic
        for i=1:length(PierDATA{a})
            [PierDATA{a}(i).longbars, PierDATA{a}(i).trasvbars]=MTAsimdesign(year(a),PierDATA{a}(i).sect,...
                pierType{a},PierDATA{a}(i).soll(:,1),PierDATA{a}(i).soll(:,2),PierDATA{a}(i).soll(:,3),...
                PierDATA{a}(i).soll(:,4),PierDATA{a}(i).soll(:,5),DesMat(a,1),DesMat(a,2));
        end
        toc
    end
% else
%     for a=1:n
%         tic
%         mat=samplesLHS(a,[3 4])
%         for i=1:length(PierDATA{a})
%             temp=InputReinf(:,[1:2]);
%             ind=find(all(mat'==temp'));
%             PierDATA{a}(i).longbars=InputReinf(ind,[3:4]);
%             PierDATA{a}(i).trasvbars=[0.01, 0.2]; %dummy
%         end
%         toc
%     end
% end

% substitute transv reinforcement by LHS sampling
try
    for a=1:n
        for i=1:length(PierDATA{a})
            PierDATA{a}(i).trasvbars(2)=step_trsv_bars(a);
            PierDATA{a}(i).trasvbars(1)=diam_trsv_bars(a);
        end
    end
end

salva=1;
if salva==1
    save([folderpathoutput,'/LHSData'])
end

%% Generate output and export

%Masses and effective heights
LongBars=cell(n,1);TransvBars=cell(n,1);

for a=1:n
    %Temp variables to help the calculations
    capbeamMasses=[PierDATA{a}.PPcapbeam];
    pierMasses=[PierDATA{a}.PPbent];
    deckMassesL=((G2vect(a)+G1vect(a)).*tribhorTOT{a});%for long dir
    deckMassesT=((G2vect(a)+G1vect(a)).*triblTOT{a});%for long dir
    heightCapbeam=sum(dCapbeam(:,3:4),2);
    
    %Calculate/Allocate masses
    massesL{a}=deckMassesL;
    massesT{a}=deckMassesT;
    piermasses(a,:)=pierMasses;
    capbeammasses(a,:)=capbeamMasses;
    
    %LongBars allocation and modification
    if strcmp(pierType{a}, 'CIRC')
        try            
            temp=reshape([PierDATA{a}.longbars]',[2,length(pierMasses)]);
            [Astot,crpier]=max(temp(2,:).*pi.*(temp(1,:)/2).^2); %Fond the maximum As and associate it to all the piers
            LongBars{a}=ones(size(temp)).*temp(:,crpier);
            TransvBars{a}=reshape([PierDATA{a}.trasvbars]',[2,length(pierMasses)]);
            gd(a)=1; %good case?=1
        catch
            gd(a)=0; %bar case?=0 , no reinforcements are present
        end
    else
        gd(a)=1; %This section for RECTsection is not implemented
        tempR(:,1,1)=PierDATA{a}.longbars(1,:)'; tempR(:,1,2)=PierDATA{a}.longbars(2,:)';
        LongBars{a}=tempR;
        TransvBars{a}=PierDATA{a}.trasvbars';
    end
    %Fixity of devices
    fixitylongCell{a}=horzcat([NaN; abxfix(a,1)], bearFixLong(:,:,a), [abxfix(a,2); NaN]);
    fixitytransvCell{a}=horzcat([NaN; abyfix(a,1)], ones(2,length(deckMassesT)-2), [abyfix(a,2); NaN]);
end

%Generate Table
Table=table(       StructType, year,  roadDim,  Hdeck,  lengthSpans, [G1vect, G2vect],pierDim,  pierType,   massesL',     massesT',      piermasses, capbeammasses,   Hpier,  heightCapbeam,[fc, fy],LongBars,   TransvBars,...
    'VariableNames',{'Struct','year','roadDim','Hdeck','lengthSpans','Loads',        'PierDim','PierType','DeckMassesL','DeckMassesT','pierMasses','capbeamMasses', 'hPier', 'hCapbeam', 'Mat','LongBars', 'TransvBars',});

Table.NumGird=NumGird;
Table.fixL=fixitylongCell'; 
Table.fixT=fixitytransvCell';
Table.bearTypeL=bearTypeLong';
Table.bearTypeT=bearTypeTransv';
Table.abTypeL=abxtype;
Table.abTypeT=abytype;
Table.bearHeight=bearHeight;
Table.backwall=[heightBkw widthBkw, kBkw];
Table.gapAbut=gapAbut;
Table.hBear=Hbear;
if bearStiffFixed>0
    Table.bearStiffFixed=bearStiffFixed;
end
if bearStiffFixed>0
    Table.bearStiffNeo=bearStiffNeo;
end

%Delete bad cases
TableRED=Table;
TableRED(gd==0,:)=[];

salvaTable=1;
if salvaTable==1
    save([folderpathoutput,'/Tables'],'Table','TableRED','PierDATA')
end

% clear fixitylongCell fixity fixitytransvCell bear_height bearDimAb bearDimPier bearFixLong bearFixTransv bearHeight bearStiffFixed bearStiffNeo...
%     bearTypeLong bearTypeTransv buildPer capbeammasses catTraff dCapbeam DesMat distBear fc fy G1vect G2vect...
%     gapAbut gd heightBkw widthBkw heightCapbeam lengthSpans LongBars massesL massesT NumGird pierDim piermasses...
%     pierType seismData sigmac sigmas StructType TransvBars tribhor tribhorTOT tribl triblTOT roadDim

%% Generation of simplified beam-SDOF models: FD laws attainment

cover=0.04; %[m]
for a=1:size(TableRED,1)
    %Calcuate effective heights
    Mcapbeam=TableRED.capbeamMasses(a,:);
    Mpier=TableRED.pierMasses(a,:);
    MdeckT=sum(TableRED.DeckMassesT{a}(:,2:end-1));
    MdeckL=sum(TableRED.DeckMassesL{a}(:,2:end-1));
    Hbear=TableRED.hBear(a,:);
    Hpier=TableRED.hPier(a,:);
    Hcapbeam=TableRED.hCapbeam(a);
    Hdeck=TableRED.Hdeck(a,:);
    
    if strcmp(StructType{a},'ISO') %if ISO the eff length is 0-Hcapbeam (later you can add the lengthof the bearings)
        EffHeightT(a,:)=((Mpier*0.3+Mcapbeam).*(Hpier-Hcapbeam/2)+...
            MdeckT.*(Hpier+Hbear+0.7*sum(Hdeck)))./(Mpier*0.3+Mcapbeam+MdeckT);
        EffHeightL(a,:)= (Hpier-Hcapbeam/2);
    elseif  strcmp(StructType{a},'HYP') %if CONT the eff length is calculated as Pinto (2009)
        EffHeightT(a,:)=((Mpier*0.3+Mcapbeam).*(Hpier-Hcapbeam/2)+...
            MdeckT.*(Hpier+Hbear+0.7*sum(Hdeck)))./(Mpier*0.3+Mcapbeam+MdeckT);
        EffHeightL(a,:)=((Mpier*0.3+Mcapbeam).*(Hpier-Hcapbeam/2)+...
            MdeckT.*(Hpier+Hbear))./(Mpier*0.3+Mcapbeam+MdeckT);
        warning('The routine for risk assessment for continuous-deck (HYP) structures is under development!')
    end
    EffHeight{a}=[EffHeightT(a,:);EffHeightL(a,:)];
end


for dr=1:2
    bilinfd=cell(size(TableRED,1),1);
    fd=cell(size(TableRED,1),1);
    buckl=cell(size(TableRED,1),1);
    mc=cell(size(TableRED,1),1);
    bilinmc=cell(size(TableRED,1),1);
    
    for a=1:size(TableRED,1)
        disp(['Processing bridge realisation #', num2str(a),' dir=',num2str(dr)])
        tic
        if strcmp(TableRED.PierType(a),'CIRC')
            parfor j=1:length(TableRED.PierDim(a,:))
                geometry=[TableRED.PierDim(a,j) cover EffHeight{a}(dr,j)]*1000; %[mm]
                longreinf=[TableRED.LongBars{a}(2,j) TableRED.LongBars{a}(1,j)*1000]; %[- mm]
                transvreinf=[TableRED.TransvBars{a}(1,j)*1000 TableRED.TransvBars{a}(2,j)*1000]; %[mm mm]
                material=[TableRED.Mat(a,1) TableRED.Mat(a,2) TableRED.Mat(a,2) 200000 0 ];
                P=PierDATA{a}(j).soll(1,1)
                [mc{a,j},bilinmc{a,j},bilinfd{a,j},fd{a,j},buckl{a,j},~,~]=CUMBIAcircfun(geometry,longreinf,...
                    transvreinf,material,P,'plot',0); %% check buckling effects! and update the bilinear until buckl is reached                
            end
        elseif strcmp(TableRED.PierType(a),'RECT')
            parfor j=1:size(TableRED.PierDim(a,:,:),2)
                matr=[2,1]; %for index inversion

                geometry=[TableRED.PierDim(a,j,dr) TableRED.PierDim(a,j,matr(dr)) cover EffHeight{a}(dr,j)]*1000; %[mm]
                ReinfLayer=distrReinf(TableRED.PierDim(a,j,dr),TableRED.LongBars{a}(:,j,dr), TableRED.LongBars{a}(:,j,matr(dr)), cover);
                longreinf=ReinfLayer %[TableRED.LongBars{a}(2,j,dr) TableRED.LongBars{a}(1,j,dr)*1000]; %[- mm]
                transvreinf=[TableRED.TransvBars{a}(1,j)*1000 TableRED.TransvBars{a}(2,j)*1000, TableRED.TransvBars{a}(3,j) TableRED.TransvBars{a}(4,j)]; %[mm mm]
                material=[TableRED.Mat(a,1) TableRED.Mat(a,2) TableRED.Mat(a,2) 200000 0 ];
                P=PierDATA{a}(j).soll(1,1)
                [mc{a,j},bilinmc{a,j},bilinfd{a,j},fd{a,j},buckl{a,j},~,~]=CUMBIArectfun(geometry,longreinf,...
                    transvreinf,material,P,1500,'noplot',0,'single', 'biaxial','n'); %% check buckling effects! and update the bilinear until buckl is reached
            end
        end
        toc
    end
    TableFDcumb{dr} = table(mc,bilinmc,bilinfd,fd,buckl,...
        'VariableNames',{'mc','bilinmc','bilinfd','fd','buckl'});
end


%% Bearings and abutments Definition

%Assumptions:
%same pier type for each subassembly, only the type of the fixed one is specified,...
% if there are free displ bearings, it is accounted in fixity condition
% (1 or 0)
%same pier dimensions are assumed for each subassemly (to change)

%Longitudinal direction
for a=1:size(TableRED,1)
    NumbBear=TableRED.NumGird(a);
    fixitylongCellRED{a}=TableRED.fixL{a};
    deckMasses=(TableRED.DeckMassesT{a}); %deckmassesT is for simulate the vertical load
    %Temp vector with index of bearing type
    BearMatrix=[TableRED.abTypeL(a,1), TableRED.bearTypeL(a,:), TableRED.abTypeL(a,2)];
    
    for pr=1:length(BearMatrix)
        if BearMatrix(pr)==2
%             bearDim=TableRED.bearDimPier(a,:);
            k=(TableRED.bearStiffNeo(a)/TableRED.bearHeight(a))*TableRED.NumGird(a); %Stiffness of a single neoprene bearing (Cardone) [kN/m]
            Fr(:,pr)=min(0.45*deckMasses(:,pr),1.5*0.4*k); %Strength of a single line of neoprene bearing (Cardone), 0.7 is the friction coefficient by Cardone 2011, 0.4 for Cardone 2013 [kN]
            DeltaY(:,pr)= Fr(:,pr)/k; DeltaUl(:,pr)=0.4*ones(2,1);
            Fu(:,pr)=Fr(:,pr); %Strength at collapse
        elseif BearMatrix(pr)==1
            %
            k=min([10^12,0.1*10^9*TableRED.NumGird(a)])*ones(2,1); %Excluded as a mechanism of collapse
            Fr(:,pr)=min([2000000,500000*TableRED.NumGird(a)])*ones(2,1); %Strength and stiffness of a single line of fixed bearings (Cardone), if you want inelastic response of fixed bearings
            DeltaY(:,pr)= Fr(:,pr)./k; DeltaUl(:,pr)=0.4*ones(2,1); %DeltaUt=0.15*ones(2,length(DeltaY));
            Fu(:,pr)=Fr(:,pr)./DeltaY(:,pr).*DeltaUl(:,pr); %Strength at collapse
        end
    end
    
    %output
    DeltaUl(1,1)=NaN; DeltaUl(2,end)=NaN;
    DeltaY(1,1)=NaN; DeltaY(2,end)=NaN; 
    FbearYL{a}=Fr.*fixitylongCellRED{a}; 
    DeltaYbearL{a}=DeltaY; 
    DeltaUbearL{a}=DeltaUl;
    FbearUL{a}=Fu.*fixitylongCellRED{a};
end
clear Fr Fy DeltaY DeltaUl

%Transverse direction
for a=1:size(TableRED,1)
    NumbBear=TableRED.NumGird(a);
    fixitytransvCellRED{a}=TableRED.fixT{a};
    deckMasses=(TableRED.DeckMassesT{a}); %deckmassesT is for simulate the vertical load
    %Temp vector with index of bearing type
    BearMatrix=[TableRED.abTypeT(a,1), TableRED.bearTypeT(a,:), TableRED.abTypeT(a,2)];
    
    for pr=1:length(BearMatrix)
        if BearMatrix(pr)==2
            k=(TableRED.bearStiffNeo(a)/TableRED.bearHeight(a))*TableRED.NumGird(a); %Stiffness of a single neoprene bearing (Cardone) [kN/m]
            Fr(:,pr)=min(0.45*deckMasses(:,pr),1.5*0.4*k); %Strength of a single line of neoprene bearing (Cardone), 0.7 is the friction coefficient by Cardone 2011, 0.4 for Cardone 2013 [kN]
            DeltaY(:,pr)= Fr(:,pr)/k; DeltaUt(:,pr)=0.4*ones(2,1);
            Fu(:,pr)=Fr(:,pr); %Strength at collapse
        elseif BearMatrix(pr)==1
            k=min([10^12,0.125*10^9*TableRED.NumGird(a)])*ones(2,1);
            Fr(:,pr)=min([2000000,500000*TableRED.NumGird(a)])*ones(2,1); %Strength and stiffness of a single line of fixed bearings (Cardone)
            DeltaY(:,pr)= Fr(:,pr)./k; DeltaUt(:,pr)=0.4*ones(2,1); %DeltaUt=0.15*ones(2,length(DeltaY));
            Fu(:,pr)=Fr(:,pr)./DeltaY(:,pr).*DeltaUt(:,pr); %Strength at collapse
        end
    end
    
    %output
    DeltaUt(1,1)=NaN; DeltaUt(2,end)=NaN; %Abutments
    DeltaY(1,1)=NaN; DeltaY(2,end)=NaN; 
    FbearYT{a}=Fr.*fixitytransvCellRED{a}; DeltaYbearT{a}=DeltaY; DeltaUbearT{a}=DeltaUt;
    FbearUT{a}=Fu.*fixitytransvCellRED{a};
end
clear Fr Fy DeltaY DeltaUt

%Table output
 TableBEAR{2} = table(fixitylongCellRED',FbearYL',FbearUL',DeltaYbearL',DeltaUbearL',...
        'VariableNames',{'fix','FbearY','FbearU','DeltaY','DeltaU'});  
 TableBEAR{1} = table( fixitytransvCellRED',FbearYT',FbearUT',DeltaYbearT',DeltaUbearT',...
        'VariableNames',{'fix','FbearY','FbearU','DeltaY','DeltaU'}); 
 
clear FbearYL FbearUL FbearYT FbearUT DeltaYbearL DeltaUbearL DeltaYbearT DeltaUbearT

%% Modify force-displacement relation for lap-splice, buckling, pdeltaeffects

%%%%%%%%%% 1) Find F-D relation for each pier

%Include lap-splice effect?
lapsplice=0;
% Modified capacity curves with lap splice
if lapsplice==1    
    for cs=1:length(TableRED.PierDim)
        if strcmp(TableRED.PierType{cs}, 'CIRC')
            for pr=1:size(TableRED.PierDim,2)
                D=TableRED.PierDim(cs,pr)*10^3;
                step=TableRED.TransvBars{cs}(2,pr)*10^3; %step of transv hoops
                dbh=TableRED.TransvBars{cs}(1,pr)*10^3; %step of transv hoops
                fy=TableRED.Mat(cs,2); %Maximum bar stress to be transferred
                fc=TableRED.Mat(cs,1);
                fsh=0.0015*210000;
                nbars=TableRED.LongBars{cs}(2,pr);dbl=TableRED.LongBars{cs}(1,pr)*10^3;
                p=min(pi*(TableRED.PierDim(cs,pr)*10^3-cover)/(2*nbars)+2*(dbl+cover),2*2^0.5*(cover+dbl)); %perimeter of each bar %Priestley 1996 pag 398
                ls=0.48*dbl*fy/fc^0.5; %Optimum splice length
                rhoh(cs,pr)=1.4*pi*(dbl/2)^2/(p*ls)*fy/fsh; %optimum amount of transv bars
                ds   = D - 2*cover + dbh;
                Ash=0.25*pi*(dbh^2)
                rhos(cs,pr)=4*Ash/(ds*step);
                
                %M0 moment of unconfined section... to be continued
                
            end
        end
    end
end
            
%%%% Buckling - Modify force-displacement relationships
buckling=1;
if buckling==1
    for dr=1:2
        FD=cell(length(TableFDcumb{dr}.fd),size(TableFDcumb{dr}.fd,2));
        for cs=1:length(TableFDcumb{dr}.fd)
            for pr=1:size(TableFDcumb{dr}.fd,2)
                
                bucklDispl=TableFDcumb{dr}.buckl{cs,pr}(2,1); %Only Berry and Eberhard 2005
                FDtemp=TableFDcumb{dr}.bilinfd{cs,pr}; %Allocate previous-calculated FD
                
                if isnan(bucklDispl)==0
                    bucklForce=interp1(FDtemp(:,1),FDtemp(:,2),bucklDispl,'linear','extrap');
                    if bucklDispl<TableFDcumb{dr}.bilinfd{cs,pr}(2,1)
                        FD{cs,pr}=[0 0 ; bucklDispl bucklForce];
                    elseif bucklDispl>TableFDcumb{dr}.bilinfd{cs,pr}(2,1) && bucklDispl<TableFDcumb{dr}.bilinfd{cs,pr}(3,1)
                        FD{cs,pr}=[FDtemp(1:2,:); bucklDispl bucklForce];
                    end
                else                   
                    FD{cs,pr}=FDtemp;
                end
                clear FDtemp               
            end
        end
        FDpier{dr}=FD;
        clear FD
    end
end

%%%% P-Delta effect - Modify force-displacement relationships
pdelta=0;
if pdelta==1
    for dr=1:2
        for cs=1:length(TableFDcumb{dr}.fd)
            deckMasses=sum(TableRED.DeckMassesT{a});
            for pr=1:size(TableFDcumb{dr}.fd,2)
                %calculate coefficient of pdelta
                theta=FDpier{dr}{cs,pr}(3,1)*deckMasses(pr)./(FDpier{dr}{cs,pr}(3,2)*EffHeight{cs}(dr,pr)); %theta at collapse
                if theta>=0.1
                    Vreduct=FDpier{dr}{cs,pr}(:,1)*deckMasses(pr)./EffHeight{cs}(dr,pr); %Allocate previous-calculated FD
                    FDpier{dr}{cs,pr}(:,2)=[FDpier{dr}{cs,pr}(:,2)-Vreduct];
                end               
            end
        end
    end
end
    
%Find FDlaws for each bearing
for dr=1:2
    for cs=1:length(FDpier{dr})
        for pr=1:length(TableRED.lengthSpans(cs,:))+1
            NotNaNind=find(~isnan(TableBEAR{dr}.FbearY{cs}(:,pr)),pr);
            FDbear{dr}{cs,pr}=[0 0;...
                min(TableBEAR{dr}.DeltaY{cs}(:,pr)) sum(TableBEAR{dr}.FbearY{cs}(NotNaNind,pr));...
                min(TableBEAR{dr}.DeltaU{cs}(:,pr)) sum(TableBEAR{dr}.FbearU{cs}(NotNaNind,pr))];
        end
    end
end

%Find FDlaws for abutment backfill
for cs=1:size(TableRED,1)
    if size(TableRED.backwall(cs,:),2)==2 | TableRED.backwall(cs,3)==[]
        kAbBkw=11.5*10^3*TableRED.backwall(cs,2)*TableRED.backwall(cs,1)/1.7;
    else
        %Nielson 2005
        %the ultimate def (delta/Hwall) is proportional to the stiffness k
        deltaUlt=(0.0334+0.000231.*TableRED.backwall(cs,3)).*TableRED.backwall(cs,1);
        kAbBkw=TableRED.backwall(cs,3)*TableRED.backwall(cs,1)/1.7; %[kN/cm per meter of wall]
        %Assuming that also the strength is proportional to the stiffness k       
        %coeff=(TableRED.backwall(cs,3)-115)/(288-115); %[kPa]
        Press=0.277*TableRED.backwall(cs,3)+183.10;
        FmaxAbBkw=Press*TableRED.backwall(cs,1)^2/1.7; %[kN per meter of wall]
        
    end
    %FabBkw=TableRED.backwall(cs,1)^2*TableRED.backwall(cs,2)/1.7*239;
    %FabBkw=TableRED.backwall(cs,1)*TableRED.backwall(cs,2)*239; %Nielson
    FDabut{cs}=[0 0; TableRED.gapAbut(cs) 0.001 ; TableRED.gapAbut(cs)+FmaxAbBkw/(kAbBkw*100) FmaxAbBkw*TableRED.backwall(cs,2); ...
        TableRED.gapAbut(cs)+deltaUlt FmaxAbBkw*TableRED.backwall(cs,2)];
    

end

%% saving subassembly data
fileout=[folderpathoutput,'/SubassembliesInput'];
salvaSUB=1;

% for  cs=1:size(TableRED,1)
% coeffHeff{cs}=EffHeight{cs}./H_piers;
% end

if salvaSUB == 1
    save(fileout,'TableRED', 'TableFDcumb','TableBEAR','FDpier','FDbear','FDabut','bilinmc') %,'coeffHeff')
else
    disp ('no save active!!')
end
