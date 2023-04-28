function [DCR,PP,crSub]=bridgeSSdcr (ADRSspectra,capcurve,deltaY,csiSDOF,EDPds)

%% CDR for simply supported bridges
%This function performes CSM assessment (modified version) for each
%subassembly cap curve and ADRSspectra to calculate CDR for different limit
%state 
for pr=1:length(capcurve)
    ind=capcurve{pr}(:,2)>=0;
    capcurve{pr}=capcurve{pr}(ind,:);
    csisdof{pr}=csiSDOF(ind,pr);
end


for pr=1:length(capcurve)  
    parfor gm=1:length(ADRSspectra) %parfor working also
        PP{pr,gm} = modCSM( ADRSspectra{gm}, capcurve{pr}, deltaY(pr), 'noplot' , 0.5, 0.0005, 3, 1, csisdof{pr}*100 );
        EDPgm(pr,gm)=PP{pr,gm}(1);
    end
end

%% Define DCR as in Borzi,2015

%Find DCR for each gm as the maximum among all the elements (pier-bearings)
for gm=1:length(ADRSspectra)
    if length(EDPgm(:,gm))==1
        DCR(:,gm)=(EDPgm(:,gm)./EDPds);
        crSub(:,gm)=1;
    else
        [DCR(:,gm),crSub(:,gm)]=max(EDPgm(:,gm)./EDPds);
    end
end
