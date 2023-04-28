function [comb,soll]=loadanalysisbridge(S,W_deck,W_walk,W_road,H_girder,Hbent,H_slab,...
    PPbent,PPcapbeam,tribl,tribhor,PPdeck,cat_traffic,seismDATA,d_bear,sigmac,piertype,sect)

%% This function perform load analysis of bridges according to the Italian codes from 1939 to 1990

%Input arrangement

%Piertype
if strcmp(piertype,'RECT')
    Lx=sect(1);
    Ly=sect(2);
elseif strcmp(piertype,'CIRC')
    Lx=sect(1);
    Ly=sect(1);
end

%Seismic Data
if S==1 || S==2 || S==3
    cat_seism=seismDATA;
elseif S==4 || S==5 || S==6
    coeff_seism=seismDATA;
end
    
%% RUN
switch S
    case 1 % S=1 Carichi Ponti R.D. 1962 (Da cambiare) - Sismica R.D. 1927 (e aggiornamenti)
        %Load analysis
        [comb,soll]=analisi_carichi1962c(W_deck,W_walk,W_road,H_girder,Hbent,H_slab,...
            PPbent,PPcapbeam,tribl,tribhor,PPdeck,Lx,cat_traffic,cat_seism,d_bear);
        
    case 2  % S=2 Carichi Ponti DM1962 - Sismica R.D. 1927 (e aggiornamenti)
        %Load analysis
        [comb,soll]=analisi_carichi1962c(W_deck,W_walk,W_road,H_girder,Hbent,H_slab,...
            PPbent,PPcapbeam,tribl,tribhor,PPdeck,Lx,cat_traffic,cat_seism,d_bear);
        
    case 3 % S=3 Carichi Ponti DM1962 - Sismica DM1975
        %Load analysis
        [comb,soll]=analisi_carichi1975c(W_deck,W_walk,W_road,H_girder,Hbent,H_slab,...
            PPbent,PPcapbeam,tribl,tribhor,PPdeck,Lx,Ly,cat_traffic,cat_seism,sigmac,d_bear,piertype);
        
    case 4 %S=4 Carichi Ponti DM1980 - Sismica DM1975
        %Load analysis
        [comb,soll]=analisi_carichi1980c(W_deck,W_walk,W_road,H_girder,Hbent,H_slab,...
            PPbent,PPcapbeam,tribl,tribhor,PPdeck,Lx,Ly,cat_traffic,coeff_seism,sigmac,d_bear,piertype,1,'noplot');
        
    case 5 %S=5 Carichi Ponti DM1980 - Sismica DM1986
        %Load analysis
        [comb,soll]=analisi_carichi1980c(W_deck,W_walk,W_road,H_girder,Hbent,H_slab,...
            PPbent,PPcapbeam,tribl,tribhor,PPdeck,Lx,Ly,cat_traffic,coeff_seism,sigmac,d_bear,piertype,1.4,'noplot');
        
    case 6  %S=6 Carichi Ponti DM1990 - Sismica DM1986
        [comb,soll]=analisi_carichi1990c(W_deck,W_walk,W_road,H_girder,Hbent,H_slab,...
            PPbent,PPcapbeam,tribl,tribhor,PPdeck,Lx,Ly,cat_traffic,coeff_seism,sigmac,d_bear,piertype,1.4);
end