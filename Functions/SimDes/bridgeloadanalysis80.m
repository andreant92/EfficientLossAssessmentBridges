function [NsACC,MsxACC,MsyACC,Fx_traffic,My_traffic,Fy_wind,Mx_wind]=bridgeloadanalysis80(tributL,ldeck,l_marciapiede,l_sedestrad,h_girder,h_bent,h_slab,Lx_pier,categoria_ponte,appoggio_testapila,dist_appoggi,Npp)

%% Analisi dei carichi, anno 1980
%convenzione assi: asse x è l'asse longitudinale dell'impalcato
%  asse y è l'asse trasversale dell'impalcato

%La larghezza di ingombro convenzionale è stabilita per ciascuna colonna di carico in m 3,50.
%Definizione categorie ponte.
%Per i ponti di I categoria si devono considerare:
%una colonna di carichi q1A
%una colonna di carichi q1B
%altre ulteriori colonne di carichi q1B compatibili con la larghezza della carreggiata di intensità ridotta del 30%;
% carico q1F sui marciapiedi.

%Per i ponti di II categoria si devono considerare:
%due colonne di carichi q1B
%altre ulteriori colonne di carichi q1B compatibili con la larghezza della carreggiata di intensità ridotta del 30%;
% carico q1F sui marciapiedi.


% Quando diano luogo a condizioni di carico più gravoso si devono sostituire:
% — alla colonna di carico q1A il carico q1C per i ponti di I categoria;
% — ad una colonna di carico q1B il carico q1D per i ponti di Il categoria.

%% Load definition

%Calcolo dei carichi q1A e q1B rappresentate da un carico ripartito disposto,
%lungo l’asse delle corsie di ingombro.
if mean(tributL*2)<=40%m
    q1A=(2.89+ 52/(sum(tributL)))*10; %kN/m
elseif mean(tributL*2)>=40 && mean(tributL*2)<=400 %m
    q1A=(4.35-(sum(tributL))/250)*10; %kN/m
elseif mean(tributL*2)>400 %m
    q1A=2.75*10; %kN/m
end

if mean(tributL*2)<=15%m
    q1B=(0.40+27/(sum(tributL)))*10;  %kN/m
elseif mean(tributL*2)>=15 && mean(tributL*2)<=400 %m
    q1B=(2.23-(sum(tributL))/500)*10; %kN/m
elseif mean(tributL*2)>400 %m
    q1B=1.43*10; %kN/m
end

%Per quanto riguarda la disposizione trasversale dei carichi,
%le condizioni peggiori si hanno nelle  2 casi seguenti: max carico(per
%ottenere il max sforzo normale e la condizione di max eccentricità

% Definizione delle corsie di carico
if l_sedestrad>=7
    n_corsie= 1+floor((l_sedestrad-3.5)/3.5);
elseif l_sedestrad>5 || l_sedestrad<7
    n_corsie = 2;
else
    n_corsie=1; 
end

Pp=Npp/(sum(tributL));
%Incremento dinamico di ciascun carico mobile (media tra le 2 campate)
if categoria_ponte==1
    fi_d=1.4-((Pp/(q1B*n_corsie+0.4*l_marciapiede)+1)*0.002*sum(tributL));%incremento dinamico del carico mobile
elseif categoria_ponte==2
    fi_d=1.4-((Pp/(q1A+q1B*(n_corsie-1)+0.4*l_marciapiede)+1)*0.002*sum(tributL));%incremento dinamico del carico mobile
end

%bridge structural typology definition
if dist_appoggi==0
    bridgetype='cont';
elseif dist_appoggi~=0
    bridgetype='simplysupp';
else
    warning('check number of bearings')
end

%Other loads
q1F=4*l_marciapiede;% kN/m %carico folla
Q1C=550;% kN %carico sostitutivo di q1a
Q1D=310;% kN %carico sostitutivo di q1b

%% Service loads


%initialization and transverse load distribution definition (maximum load and
%maximum eccentricity)

%Definisco vettore di carichi distribuiti  Qdist(una componente su ogni corsia)
%es. se ci sono 3 corsie avrò [qfolla q1corsia q2corsia q3corsia qfolla] 
%MC= max carico, ME=max ecc
QdistMC=zeros(2,(n_corsie+2));
QdistMC(:,1)=q1F.*tributL.*fi_d; %[kN]
QdistMC(:,end)=q1F.*tributL.*fi_d;  %[kN]

if categoria_ponte==1
    % vettore di carichi agenti sull'impalcato  %kN/m 
    %Qvectdist per carichi distribuiti
    %Qvectconc per carichi accidentali
    %MC per Massimo Carico - ME per Massima Eccentricità
    
    QdistMC(:,2)=max(q1A.*tributL.*fi_d,Q1C.*fi_d); %è errato perchè prende Q1c su entrambe le corsie
    %aggiungi un if per togliere l'eventuale Q1c su entrambe le campate
    try
        QdistMC(:,3)= q1B.*tributL.*fi_d;
        QdistMC(:,4:end-1)=0.7*q1B.*tributL.*fi_d;
    end
    
    %Quando diano luogo a condizioni di carico più gravoso si devono sostituire:
    %— alla colonna di carico q1A il carico q1C per i ponti di I categoria;
    %— alla colonna di carico q1B il carico q1D per i ponti di II categoria;
        
elseif categoria_ponte==2
        
    QdistMC(:,2)=max(q1B.*tributL.*fi_d,Q1D.*fi_d);
    try
        QdistMC(:,3)= q1B.*tributL.*fi_d;
        QdistMC(:,4:end-1)=0.7*q1B.*tributL.*fi_d;
    end
    
%     QdistMCbis=QdistMC;
%     QdistMCbis(:,2)=Q1D.*fi_d;
end

QdistME=zeros(2,(n_corsie+2));
QdistME(:,1:2)=QdistMC(:,1:2);

%Definition of longitudinal load distribution 
if strcmp(bridgetype,'simplysupp')
    
    %Distributed loads
    MATRIXloadD={QdistMC,QdistMC,... %carico entrambe le campate MC
        [QdistMC(1,:);zeros(1,length(QdistMC))],...%carico solo la prima campata MC
        [zeros(1,length(QdistMC));QdistMC(1,:)],...
        QdistME,QdistME,... %carico entrambe le campate ME
        [QdistME(1,:);zeros(1,length(QdistME))],...%carico solo la prima campata ME
        [zeros(1,length(QdistME));QdistME(1,:)]}; %carico solo la seconda campata ME;

elseif strcmp(bridgetype,'cont')
    
    %Distributed loads
    MATRIXloadD={sum(QdistMC),... %carico entrambe le campate MC
        sum(QdistME)}; %carico entrambe le campate ME};
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
    NsACC(i)=sum(MATRIXloadD{i}(:));
    %transverse moment
    momtempx=MATRIXloadD{i}.*ecc;
    MsxACC(i)=sum(momtempx(:));
    %longit moment
    
    if strcmp(bridgetype,'simplysupp')
        momtempy=[sum(MATRIXloadD{i}(1,:))*dist_appoggi(1);...
            sum(MATRIXloadD{i}(2,:))*(-dist_appoggi(2))];
    else
        momtempy=0;
    end
    
    MsyACC(i)=abs(sum(momtempy));
end

clear momtempx momtempy

%% AZIONE FRENANTE
%Si calcola, la forza di frenatura e il momento longitudinale prodotto
%dalla stessa;si assume agente nella direzione dell’asse della strada ed al livello della
%superficie stradale, con intensità pari ad 1/10 della singola colonna di carico più pesante per
%ciascuna carreggiata. Tale carico viene paragonato col 25% del carico
%concentrato Q1D potenzialmente presente.

LoadC=zeros(2,length(MATRIXloadD)); %Vettore di carichi concentrati da paragonare con MatrixLoadD per il carico frenante
if strcmp(appoggio_testapila{1},'fixed')
    LoadC(1,:)=[0.25*Q1D 0 0.25*Q1D 0 0.25*Q1D 0 0.25*Q1D 0];
end
if strcmp(appoggio_testapila{2},'fixed')
    LoadC(2,:)=[ 0 0.25*Q1D 0 0.25*Q1D 0 0.25*Q1D 0 0.25*Q1D];
end
%LoadC is like MATRIXloadC of bridgeloadanalysis90 in traffic
%load calulation
%LoadC is  #supports row and #loadcases column

for k=1:length(MATRIXloadD)
    for i=1:length(dist_appoggi)
        if strcmp(appoggio_testapila{i},'fixed')
            %Caso 1: Forza frenante legata alla colonna di carico più pesante
            Fx_trafTd(i,k)=1/10*(MATRIXloadD{k}(i,2)/fi_d); %kN, forza di frenatura %dovrei mettere lunghezza di competenza orizzontale per ponti continui
            %Caso 2: Forza frenante legata a tutti i carichi tandem
            %contemporaneamente presenti
%             Fx_trafTc(i,k)=0.25*Q1D; %kN, forza di frenatura
        elseif strcmp(appoggio_testapila{i},'free')
            Fx_trafTd(i,k)=0;
%             Fx_trafTc(i,k)=0;
        else
            warning('error in appoggio_testapila')
        end
    end
end
%Prendo il massimo tra caso1 e caso2
if strcmp(bridgetype,'simplysupp')
    FxT_D=sum(Fx_trafTd);FxT_C=sum(LoadC); %somma i contributi delle campate dx e sx se è semplic appogg
else
    FxT_D=Fx_trafTd;FxT_C=LoadC;
end
for k=1:length(FxT_D)
    Fx_traffic(k)=max([FxT_D(k),FxT_C(k)]);
end
My_traffic=(Fx_traffic*h_bent); %kNm,momento flettente long. alla base della pila


% %% AZIONE FRENANTE
% %Si calcola, la forza di frenatura e il momento longitudinale prodotto dalla stessa;
% % si assume agente nella direzione dell’asse della strada ed al livello della
% % superficie stradale, d'intensità pari ad  1/10 del sovraccarico costituito da q1B.
% if strcmp(bridgetype, 'simplysupp')   
%     for i=1:length(dist_appoggi)
%         if strcmp(appoggio_testapila{i},'fixed')
%             if categoria_ponte==1
%                 Fx_trafT(i)=(1/10*(q1B*(tributL(i)*2)).*fi_d); %kN, forza di frenatura
%             elseif categoria_ponte==2
%                 Fx_trafT(i)=(1/10*(q1B*(tributL(i)*2)).*fi_d); %kN, forza di frenatura
%             end
%         elseif strcmp(appoggio_testapila{i},'free')
%             Fx_trafT(i)=0;
%             My_traffic(i)=0;
%         else
%             warning('error in appoggio_testapila')
%         end
%     end
% else
%     Fx_trafT=(1/10*(q1B*sum(tributL)).*fi_d); %kN, forza di frenatura
% end
% Fx_traffic=sum(Fx_trafT);
% My_traffic=(Fx_traffic*h_bent); %kNm,momento flettente long. alla base della pila

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