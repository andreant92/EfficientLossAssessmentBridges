function [NsACC,MsxACC,MsyACC,Fx_traffic,My_traffic,Fy_wind,Mx_wind]=bridgeloadanalysis90(tributL,ldeck,l_marciapiede,l_sedestrad,h_girder,h_bent,h_slab,Lx_pier,categoria_ponte,appoggio_testapila,dist_appoggi)

%% Analisi dei carichi, anno 1990

%convenzione assi: asse x è l'asse longitudinale dell'impalcato
%                  asse y è l'asse trasversale dell'impalcato
% La larghezza di ingombro convenzionale è stabilita per ciascun colonna in
% 3,50 m.
% Devono essere presi in considerazione i seguenti carichi mobili:

% In direzione trasversale il segno dell'eccentricità dev'essere definito in segno in modo coerente
% alla distribuzione trasversale dei carichi accidentali più gravosa.
% es.
% se concentro i carichi mobili a sx, e l'eccentricità dei carichi
% permanenti è discorde (a dx), quest'ultima dovrà essere negativa.

%% Load definition
Q1a=600;%kN % q1a) 	mezzo convenzionale da 60 t a tre assi.

q1b=30; %kN/m % q1b) 	carico ripartito pari a 3 t/m disposto, ai fini del calcolo delle
%strutture principali, lungo l’asse di una corsia d’ingombro;

q1e=4*l_marciapiede;%kN/m % q1e) 	carico della folla uniformemente ripartito in superficie pari a 0,4 t/m.

% Per i ponti di 1° categoria si devono considerare:
% una colonna di carico costituita da un solo mezzo q1a
% e, al di fuori dell’ingombro di questo, da uno o più tratti di carico q1b, disposti, ai fini
% del calcolo delle strutture principali, lungo l’asse della corsia nel modo più sfavorevole;
% una seconda colonna di carico analoga alla precedente, ma con carichi pari rispettivamente al
% 50% di q1a e al 50% di q1b; altre colonne di carico analoghe alle precedenti, ma con carichi
% pari rispettivamente al 35% di q1a ed al 35% di q1b; il carico q1e sui marciapiedi.

% Per i ponti di 2° categoria si devono considerare:
% una colonna di carico analoga a quella prevista per i ponti di 1° categoria,
% ma con carichi pari rispettivamente al 75% di q1a ed al 50% di q1b;
% una seconda colonna di carico analoga alla precedente, ma con carichi pari rispettivamente
% al 50% di q1a ed al 50% di q1b; altre colonne di carico analoghe alle precedenti, ma con carichi
% pari rispettivamente al 35% di q1a ed al 35% di q1b e il carico q1e sui marciapiedi.

% Qualora la larghezza della carreggiata suddetta contenga più di 4 colonne di 3,50 m devono
% prevedersi, in luogo di una sola colonna formata da q1A+q1B, due colonne cosi formate tra loro
% non contigue.

% Definizione delle corsie di carico
if l_sedestrad>=7
    n_corsie= 1+floor((l_sedestrad-3.5)/3.5);
elseif l_sedestrad>5.5 || l_sedestrad<7
    n_corsie = 2;
else
    n_corsie=1; 
end

%Incremento dinamico di ciascun carico mobile (media tra le 2 campate)

fi_d=(min(1.4,max(1,1.4-(tributL.*2-10)/150)));
% 
% Q1a_d=(fi-1).*Q1a+Q1a; %kN
% q1b_d=(fi-1).*q1b+q1b;%kN/m
% q1e_d=(fi-1).*q1e+q1e; %kN/m
%bridge structural typology definition
if length(dist_appoggi)==1
    bridgetype='cont';
elseif length(dist_appoggi)==2
    bridgetype='simplysupp';
else
    warning('check number of bearings')
end

%% Service loads

%initialization and transverse load distribution definition (maximum load and
%maximum eccentricity)

QdistMC=zeros(2,(n_corsie+2));
QdistMC(:,1)=q1e.*tributL.*fi_d;
QdistMC(:,end)=q1e.*tributL.*fi_d;
QtandMC=zeros(2,(n_corsie+2));

if categoria_ponte==1  
    % vettore di carichi agenti sull'impalcato  %kN/m 
    %Qvectdist per carichi distribuiti
    %Qvectconc per carichi accidentali
    %MC per Massimo Carico - ME per Massima Eccentricità
    
    QdistMC(:,2)=q1b.*(tributL-7.5).*fi_d;
    QtandMC(:,2)=Q1a.*fi_d;
    try
        QdistMC(:,3)= QdistMC(:,2)*0.5;
        QtandMC(:,3)= QtandMC(:,2)*0.5;
        QdistMC(:,4:end-1)=QdistMC(:,2)*0.35;
        QtandMC(:,4:end-1)= QtandMC(:,2)*0.35;
    end
    
elseif categoria_ponte==2
    
    QdistMC(:,2)=0.5*q1b.*(tributL-7.5).*fi_d;
    QtandMC(:,2)=0.75*Q1a.*fi_d;
    try
        QdistMC(:,3)= QdistMC(:,2)*0.5;
        QtandMC(:,3)= QtandMC(:,2)*0.5;
        QdistMC(:,4:end-1)=QdistMC(:,2)*0.35;
        QtandMC(:,4:end-1)= QtandMC(:,2)*0.35;
    end  
end

QdistME=zeros(2,(n_corsie+2));
QdistME(:,1:2)=QdistMC(:,1:2);
QtandME=zeros(2,(n_corsie+2));
QtandME(:,2)=QtandMC(:,2);

%Definition of longitudinal load distribution 
if strcmp(bridgetype,'simplysupp')
    
    %Distributed loads
    MATRIXloadD={QdistMC,QdistMC,... %carico entrambe le campate MC
        [QdistMC(1,:);zeros(1,length(QdistMC))],...%carico solo la prima campata MC
        [zeros(1,length(QdistMC));QdistMC(1,:)],...
        QdistME,QdistME,... %carico entrambe le campate ME
        [QdistME(1,:);zeros(1,length(QdistME))],...%carico solo la prima campata ME
        [zeros(1,length(QdistME));QdistME(1,:)]}; %carico solo la seconda campata ME;
    
    %Tandem loads    
%     [~,ind]=max(abs(dist_appoggi)); %ind è la campata dove mettere carico tandem (quando le campate sono entrambe caricate)
%     tempMC=zeros(2,length(QtandMC));
%     tempMC(ind,:)=QtandMC(ind,:);
%     tempME=zeros(2,length(QtandME));
%     tempME(ind,:)=QtandME(ind,:);
    
    MATRIXloadC={ [QtandMC(1,:);zeros(1,length(QtandMC))],...
        [zeros(1,length(QtandMC));QtandMC(2,:)],...
        [QtandMC(1,:);zeros(1,length(QtandMC))],...
        [zeros(1,length(QtandMC));QtandMC(2,:)],...
        [QtandME(1,:);zeros(1,length(QtandME))],...
        [zeros(1,length(QtandME));QtandME(2,:)],...
        [QtandME(1,:);zeros(1,length(QtandME))],...
        [zeros(1,length(QtandME));QtandME(2,:)]};

elseif strcmp(bridgetype,'cont')
    
    %Distributed loads
    MATRIXloadD={sum(QdistMC),... %carico entrambe le campate MC
        sum(QdistME)}; %carico entrambe le campate ME};
    
    %Tandem loads
    MATRIXloadC={max(QtandMC),... %massimo tra QtandMC per corsia 1 (fi(1)) e per corsia 2 (fi(2))
        max(QtandME)};  %massimo tra QtandME per corsia 1 (fi(1)) e per corsia 2 (fi(2))
end

%Definisco vettore eccentricità (dist risultante carichi da asse pila)
for m=1:n_corsie %vettore eccentricità del carico q1B
     tempECC(m)=ldeck/2-l_marciapiede-3.5/2-3.5*(m-1);
end

l_marc2=(ldeck-l_marciapiede-l_sedestrad); %lunghezza marciapiede2
eccmarc=[ldeck/2-l_marciapiede/2,-ldeck/2+l_marc2/2];

ecc=[eccmarc(1), tempECC, eccmarc(2)];
clear eccmarc ind m 

%Calcolo sollecitazioni

for i=1:length(MATRIXloadD)
    %Normal stress
    NsACC(i)=sum([MATRIXloadD{i}(:);MATRIXloadC{i}(:)]);
    %transverse moment
    momtempx=[MATRIXloadD{i}.*ecc;MATRIXloadC{i}.*ecc];
    MsxACC(i)=sum(momtempx(:));

    %longit moment    
    if strcmp(bridgetype,'simplysupp')
    momtempy=[sum(MATRIXloadD{i}(1,:))*dist_appoggi(1);...
        sum(MATRIXloadD{i}(2,:))*-dist_appoggi(2);...
        sum(MATRIXloadC{i}(1,:))*dist_appoggi(1);...
        sum(MATRIXloadC{i}(2,:))*-dist_appoggi(2)];
    else
            momtempy=0;
    end
    
    MsyACC(i)=abs(sum(momtempy));
end

clear momtempx momtempy

% %% AZIONE FRENANTE
% %Si calcola, la forza di frenatura e il momento longitudinale prodotto
% %dalla stessa;si assume agente nella direzione dell’asse della strada ed al livello della
% %superficie stradale, con intensità pari ad 1/10 della singola colonna di carico più pesante per
% %ciascuna carreggiata.
% 
% for i=1:length(dist_appoggi)
%     if strcmp(appoggio_testapila{i},'fixed')
%         if categoria_ponte==1
%             Fx_trafT(i)=(1/10*(Q1a+q1b*(tributL(i)*2-15))*fi_d(i)); %kN, forza di frenatura %dovrei mettere lunghezza di competenza orizzontale
%         elseif categoria_ponte==2
%             Fx_trafT(i)=(1/10*(0.75*Q1a+0.5*q1b*(tributL(i)*2-15))*fi_d(i)); %kg, forza di frenatura
%         end
%     elseif strcmp(appoggio_testapila{i},'free')
%         Fx_trafT(i)=0;
%         My_traffic(i)=0;
%     else
%         warning('error in appoggio_testapila')
%     end
% end
% Fx_traffic=sum(Fx_trafT);
% My_traffic=(Fx_traffic*h_bent); %kNm,momento flettente long. alla base della pila

%% AZIONE FRENANTE
%Si calcola, la forza di frenatura e il momento longitudinale prodotto
%dalla stessa;si assume agente nella direzione dell’asse della strada ed al livello della
%superficie stradale, con intensità pari ad 1/10 della singola colonna di carico più pesante per
%ciascuna carreggiata.
for k=1:size(MATRIXloadC,2)
    for i=1:length(dist_appoggi)
        if strcmp(appoggio_testapila{i},'fixed')
            if categoria_ponte==1
                %Caso 1: Forza frenante legata alla colonna di carico più pesante
                Fx_trafTd(i,k)=1/10*(MATRIXloadC{k}(i,2)/fi_d(i)+MATRIXloadD{k}(i,2)/fi_d(i)); %kN, forza di frenatura %dovrei mettere lunghezza di competenza orizzontale per ponti continui
                %Caso 2: Forza frenante legata a tutti i carichi tandem
                %contemporaneamente presenti
                Fx_trafTc(i,k)=0.2*sum(MATRIXloadC{k}(i,:)/fi_d(i)); %kN, forza di frenatura
            elseif categoria_ponte==2
                Fx_trafTd(i,k)=1/10*(MATRIXloadC{k}(i,2)/fi_d(i)+MATRIXloadD{k}(i,2)/fi_d(i)); %kN, forza di frenatura %dovrei mettere lunghezza di competenza orizzontale per ponti continui
                Fx_trafTc(i,k)=0.15*sum(MATRIXloadC{k}(i,:)/fi_d(i)); %kN, forza di frenatura
            end
        elseif strcmp(appoggio_testapila{i},'free')
            Fx_trafTd(i,k)=0;
            Fx_trafTc(i,k)=0;
        else
            warning('error in appoggio_testapila')
        end
    end
end
%Prendo il massimo tra caso1 e caso2
if strcmp(bridgetype,'simplysupp')
    FxT_D=sum(Fx_trafTd);FxT_C=sum(Fx_trafTc); %somma i contributi delle campate dx e sx se è semplic appogg
else
    FxT_D=Fx_trafTd;FxT_C=Fx_trafTc;
end
for k=1:length(FxT_D)
    Fx_traffic(k)=max([FxT_D(k),FxT_C(k)]);
end
My_traffic=(Fx_traffic*h_bent); %kNm,momento flettente long. alla base della pila

%% AZIONE VENTO
%Si calcola la forza del vento in  direzione trasversale in
%condizioni di ponte carico,considerando un carico di superficie pari a
%250kg/mq agente su una superficie di altezza pari a h_travi+h_pavimentaz.+
%3m (generica altezza di un mezzo).
    
Fy_wind_deck=2.5*tributL*(h_slab+h_girder);%kN,
Fy_wind_acc=2.5*tributL*(3);%kN,
%forza vento relativa all'impalcato applicata a metà altezza della superficie su cui agisce.
Fy_wind_bent= 2.5*Lx_pier*h_bent; %kN, forza vento che agisce sulla pila e sul pulvino applicata a metà
%altezza di quella totale

Fy_windT(1)=sum(Fy_wind_deck)+Fy_wind_bent+sum(Fy_wind_acc); %kN, se carico accidentale è presente su entrambe
%le campate
if strcmp(bridgetype,'simplysupp')
    Fy_windT(2)=sum(Fy_wind_deck)+Fy_wind_bent+sum(Fy_wind_acc); %kN, carico accidentale su entrambe le campata
    Fy_windT(3)=sum(Fy_wind_deck)+Fy_wind_bent+Fy_wind_acc(1); %kN, carico accidentale su prima campata
    Fy_windT(4)=sum(Fy_wind_deck)+Fy_wind_bent+Fy_wind_acc(2); %kN  carico accidentale su seconda campata
end
Fy_wind=[Fy_windT,...
    sum(Fy_wind_deck)+Fy_wind_bent]; %kN  carico accidentale non presente

%Le forze del vento che agiscono sull'impalcato producono momenti flettenti
%trasversali alla base della pila.
Mx_windT(1)=sum(Fy_wind_deck)*(h_bent+(h_slab+h_girder)/2)+sum(Fy_wind_acc)*(h_bent+h_slab+h_girder+1.5)+Fy_wind_bent*h_bent/2;%kNm
if strcmp(bridgetype,'simplysupp')
    Mx_windT(2)=sum(Fy_wind_deck)*(h_bent+(h_slab+h_girder)/2)+sum(Fy_wind_acc)*(h_bent+h_slab+h_girder+1.5)+Fy_wind_bent*h_bent/2;%kNm
    Mx_windT(3)=sum(Fy_wind_deck)*(h_bent+(h_slab+h_girder)/2)+(Fy_wind_acc(1))*(h_bent+h_slab+h_girder+1.5)+Fy_wind_bent*h_bent/2;%kNm
    Mx_windT(4)=sum(Fy_wind_deck)*(h_bent+(h_slab+h_girder)/2)+(Fy_wind_acc(2))*(h_bent+h_slab+h_girder+1.5)+Fy_wind_bent*h_bent/2;%kNm
end
Mx_wind=[Mx_windT,...
    (sum(Fy_wind_deck)*(h_bent+(h_slab+h_girder)/2)+Fy_wind_bent*h_bent/2)];%kN  carico accidentale non presente

clear Fy_windT My_windT

end