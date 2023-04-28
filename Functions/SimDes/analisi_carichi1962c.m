function  [comb,soll]=analisi_carichi1962c(ldeck,l_marciapiede,...
    l_sedestrad,h_girder,h_bent,h_slab,Pp_bent,Pp_bentcap,tribl,tribl_sism,Ppdeck,Lx_pier,categoria_ponte,cat_sism,d_bear,varargin)

%% INPUT elaboration

%bridge structural typology definition
if d_bear==0
    bridgetype='cont';
elseif d_bear~=0
    bridgetype='simplysupp';
else
    warning('check number of bearings')
end

if  strcmp(bridgetype, 'cont') 
    if tribl_sism==0
        beartype{1}='free';     
    else
        beartype{1}='fixed';
        tribl_sism=sum(tribl_sism);
    end
elseif strcmp(bridgetype, 'simplysupp')
    for i=1:2
        if tribl_sism(i)==0
            beartype{i}='free';
        else
             beartype{i}='fixed';
        end
    end
end

%% Permanent loads

R=(Ppdeck(1).*tribl); %kN, vector of reactions at the bearing/s
Npp=sum(R); %kN, sforzo normale dovuto ai carichi permanenti
Mpp_x=Npp*Ppdeck(2); %kN, momento flettente trasversale dovuto ai carichi permanenti
if strcmp(bridgetype, 'simplysupp')
    Mpp_y=R(1)*d_bear(1)-R(2)*d_bear(2); %kN, momento flettente trasversale dovuto ai carichi permanenti
elseif strcmp(bridgetype, 'cont')
    Mpp_y=0;
end

%% RUN
%Call load analysis consistent to bridge load code 1962
[NsACC,MsxACC,MsyACC,Fx_traffic,My_traffic,Fy_wind,Mx_wind]=bridgeloadanalysis62(tribl,ldeck,l_marciapiede,l_sedestrad,h_girder,h_bent,h_slab,Lx_pier,categoria_ponte,beartype,d_bear,Npp);

%% Seismic loads
%Call seismic load analysis consistent with seismic code 1975
[Fx_seism,Fy_seism,Fv_seism, My_seism,Mx_seism]=seismloadanalysis37(Npp,ldeck,h_bent,Pp_bent,Pp_bentcap,tribl_sism,Ppdeck,beartype,cat_sism);

%% COMBINAZIONI DI CARICO

%Solllecitazioni prodotte da tutte le forze in gioco sulla struttura nel
%seguente ordine che costistuisce il vettore: N,Mx,My,Vx,Vy.
%Si ricordi che tutte le forze e momenti sono state calcolate in valore
%assoluto, in modo da poterle considerare di verso concorde e poterle
%sommare per ottenere l'effetto piu sfavorevole in termini di sollecitazione

%G1 group (Permanent Load with wind)
N_permanent=Pp_bent+Pp_bentcap+Npp;
P=[N_permanent Mpp_x sum(Mpp_y) 0  0];  
W=[0 Mx_wind(end) 0 0 Fy_wind(end)];

G1=P+W;

%G2 AND G3 group (Service loads)
emptycol=zeros(length(NsACC),1);
tempACC=[NsACC',MsxACC',MsyACC' emptycol emptycol];
tempWIND=[emptycol, [Mx_wind(1:end-1)';Mx_wind(1:end-1)'], emptycol, emptycol, [Fy_wind(1:end-1)';Fy_wind(1:end-1)']];
tempTRAFFIC=[emptycol emptycol My_traffic.*ones(length(NsACC),1) Fx_traffic.*ones(length(NsACC),1) emptycol];

G2=P+abs(tempACC)+tempWIND;
G3=P+abs(tempACC)+tempWIND+tempTRAFFIC;

%G6 (Seismic loads) 
tempP=abs(P.*ones(4,1));
seism=[Fv_seism 0 My_seism Fx_seism 0;... %Fx_sisma
     Fv_seism Mx_seism 0 0 Fy_seism;...%Fy_sisma
    -Fv_seism 0 My_seism Fx_seism 0;...%Fx_sisma
    -Fv_seism Mx_seism 0 0 Fy_seism];%Fy_sisma
G6=tempP+seism;

%OUTPUT allocation
comb.c1=G1;
comb.c2=G2;
comb.c3=G3;
comb.c4=G6;

soll=vertcat(G1,G2,G3,G6)

end


