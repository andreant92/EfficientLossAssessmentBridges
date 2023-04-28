function [NsACC,MsxACC,MsyACC,Fx_traffic,My_traffic,Fy_wind,Mx_wind]=bridgeloadanalysis62(tributL,ldeck,l_marciapiede,l_sedestrad,h_girder,h_bent,h_slab,Lx_pier,categoria_ponte,appoggio_testapila,dist_appoggi,Npp)

%% Analisi dei carichi, anno 1962
% convenzione assi: asse x è l'asse longitudinale dell'impalcato
% asse y è l'asse trasversale dell'impalcato

%Inserisco la matrice A, di cui:

%==== 1st col: si riferisce alla lunghezza campate(m); 
%le restanti colonne fanno riferimento a carichi ripartiti longitudinalmente. 
%==== 2nd col: 1o sdc 'colonna indefinita di autcararri da 12 t' [t/m], %(corsia 3.11 m)
%==== 3rd col: 2o sdc 'rullo compressore isolato da 18 t'(t/m) %(corsia 3.11 m)
%==== 4st col: 4o sdc 'treno indefinito di carichi miilitari da 61,5 t'in %(t/ml) %(corsia 3.5 m)
%==== 5st col: 5o sdc 'treno indefinito di carichi militari da 32 t' in (t/ml) %(corsia 3.5 m)
%==== 6st col: 6o sdc 'carico militare isolato da 74,5t' in (t/ml). (corsia 3.5 m)

matrix_A=[1	16.00	24.00	28.00	12.16	36.00
    1.5	10.67	16.00	24.64	10.29	31.68
    2	8.00	12.00	20.86	9.79	26.82
    2.5	6.40	9.60	17.83	9.11	22.93
    3	5.33	8.00	15.49	8.77	19.92
    3.5	4.90	7.35	13.67	8.62	17.57
    4	4.50	6.75	12.22	8.40	15.71
    4.5	4.15	6.22	11.03	8.22	14.19
    5	3.84	5.76	10.06	7.94	12.93
    5.5	3.57	5.36	9.50	7.62	12.23
    6	3.33	5.00	9.15	7.29	11.80
    6.5	3.31	4.69	8.79	6.97	11.36
    7	3.27	4.41	8.43	6.66	10.92
    7.5	3.20	4.16	8.09	6.37	10.49
    8	2.13	3.94	7.77	6.13	10.08
    8.5	3.05	3.74	7.46	5.95	9.69
    9	2.96	3.56	7.26	5.81	9.40
    9.5	2.93	3.39	7.06	5.72	9.13
    10	2.88	3.24	6.88	5.64	8.87
    11	2.78	2.98	6.61	5.55	8.47
    12	2.67	2.75	6.33	5.51	8.08
    13	2.65	2.56	6.06	5.45	7.70
    14	2.61	2.39	5.81	5.36	7.36
    15	2.56	2.24	5.61	5.23	7.07
    16	2.53	2.11	5.41	5.12	6.80
    17	2.49	1.99	5.22	5.03	6.54
    18	2.44	1.84	5.10	4.98	6.29
    19	2.44	1.80	5.05	4.96	6.06
    20	2.42	1.71	5.01	4.94	5.84
    21	2.40	1.63	4.95	4.92	5.64
    22	2.38	1.56	4.88	4.88	5.44
    23	2.03	1.432	4.00	4.30	4.48
    24	2.03	1.378	3.94	4.29	4.36
    25	2.02	1.327	3.87	4.26	4.26
    26	2.04	1.28	3.82	4.22	4.15
    27	2.04	1.236	3.79	4.18	4.05
    28	2.04	1.196	3.77	4.15	3.95
    29	2.04	1.157	3.74	4.12	3.86
    30	2.03	1.121	3.70	4.11	3.77
    31	2.03	1.088	3.67	4.10	3.69
    32	2.03	1.056	3.63	4.11	3.00
    33	2.03	1.026	3.58	4.11	3.52
    34	2.02	0.997	3.57	4.13	3.45
    35	2.01	0.971	3.57	4.14	3.33
    36	2.01	0.945	3.58	4.16	3.30
    37	2.01	0.921	3.61	4.17	3.23
    38	2.02	0.898	3.63	4.17	3.17
    39	2.02	0.876	3.65	4.17	3.11
    40	2.02	0.856	3.66	4.17	3.04
    45	2.02	0.765	3.67	4.11	2.77
    50	2.01	0.691	3.66	4.12	2.54
    55	2.01	0.631	3.65	4.14	2.35
    60	2.00	0.58	3.62	4.11	2.18
    70	2.01	0.5	    3.57	4.12	1.90
    80	2.01	0.439	3.60	4.10	1.69
    90	2.00	0.391	3.60	4.11	1.52
    100	2.00	0.353	3.58	4.11	1.38
    120	2.00	0.295	3.58	4.11	1.16
    140	2.00	0.253	3.57	4.10	1.01
    160	2.00	0.222	3.58	4.10	0.89
    180	2.00	0.198	3.57	4.11	0.79
    200	2.00	0.178	3.57	4.09	0.72];


%La larghezza d'ingombro trasversale per lo schema 1 e 2 è pari a 3.11m,
%quella per gli schemi 4,5,6 è 3.50m.
%E' necessario specificare la categoria del ponte in esame per
%procedere all'analisi.Per quanto riguarda la disposizione trasversale dei carichi,
%le condizioni peggiori si hanno nelle  2 casi seguenti: max carico(per
%ottenere il max sforzo normale e la condizione di max eccentricità
%Se il ponte appartiene alla 1 categoria adotteremo
%uno schema militare: schema più gravoso tra 4/5/6 affiancato affiancato da una o piu colonne di
%autocarri(schema 1) e folla compatta sui marciapiedi(schema 3).Questo
%schema rappresenta la condizione a pieno carico per ottenere il max
%sforzo normale sulla pila. Se il ponte appartiene alla 2 categoria adotteremo la
%condizione più sfavorevole tra le seguenti:
%-una o più colonne indefinite di autocarri e folla compatta sui
%marciapiedi;
%-uno o più rulli compressori affiancati e folla compatta sui marciapiedi.

%Se l'impalcato è costituito da travi in semplice appoggio, si assume che
%il carico accidentale gravi sulla trave di luce maggiore e quindi la luce di calcolo 
%è la max tra l1 ed l2; diversamente sel'impalcato è costituito da uno schema di trave continua
%si assume come luce di calcolo la luce media di l1 ed l2.

% Definizione delle corsie di carico
if categoria_ponte==1
    n_corsie= 1+floor((l_sedestrad-3.5)/3.11);
elseif categoria_ponte==2
    n_corsie = floor(l_sedestrad/3.11) 
end

%elaboration
% n_corsie=floor((l_sedestrad-3.50)/3.11);%n.corsie di carico riferite allo schema 1
% n_corsie2=floor(l_sedestrad/3.11);%n.corsie di carico riferite allo schema 2

% %Equivalent load vector extraction
% matrixEQLOAD=interp1(matrix_A(:,1),matrix_A, sum(tributL));

if tributL<100
    fi_d=1+((100-sum(tributL))^2/(100*(250-sum(tributL))));%incremento carico accidentale
elseif tributL>=100
    fi_d=1;
end

%bridge structural typology definition
if length(dist_appoggi)==1
    bridgetype='cont';
elseif length(dist_appoggi)==2
    bridgetype='simplysupp';
else
    warning('check number of bearings')
end

%Other loads
q1F=4;% kN/m %carico folla
matrixEQLOAD=10*interp1(matrix_A(:,1),matrix_A(:,2:end), sum(tributL));

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
    q1A=max(matrixEQLOAD(3:5)); %carico militare più gravoso 1a cat
    q1B=matrixEQLOAD(1); %carico civile (autocarri) 
    QdistMC(:,2)=(q1A.*tributL.*fi_d);
    try
        QdistMC(:,3:end-1)= q1B.*tributL.*fi_d;
    end
    
elseif categoria_ponte==2
    
    q1B=max(matrixEQLOAD(1:2)); %carico civile più gravoso 2a cat
    QdistMC(:,2)=(q1B.*tributL.*fi_d);
    try
        QdistMC(:,3:end-1)= q1B.*tributL.*fi_d;
    end
end

QdistME=zeros(2,(n_corsie+2));
QdistME(:,1:2)=QdistMC(:,1:2);

%Definition of longitudinal load distribution 
if strcmp(bridgetype,'simplysupp')
    
    %Distributed loads
    MATRIXloadD={QdistMC,... %carico entrambe le campate MC
        [QdistMC(1,:);zeros(1,length(QdistMC))],...%carico solo la prima campata MC
        [zeros(1,length(QdistMC));QdistMC(1,:)],...
        QdistME,... %carico entrambe le campate ME
        [QdistME(1,:);zeros(1,length(QdistME))],...%carico solo la prima campata ME
        [zeros(1,length(QdistME));QdistME(1,:)]}; %carico solo la seconda campata ME;

elseif strcmp(bridgetype,'cont')
    
    %Distributed loads
    MATRIXloadD={sum(QdistMC),... %carico entrambe le campate MC
        sum(QdistME)}; %carico entrambe le campate ME};
end

%Definisco vettore eccentricità (dist risultante carichi da asse pila)
if categoria_ponte==1
    length_corsie=[0 3.5 3.11*ones(1,(n_corsie-1))];
elseif categoria_ponte==2
    length_corsie=[0 3.11 3.11*ones(1,(n_corsie-1))];
end
for m=1:n_corsie %vettore eccentricità del carico q1B
    tempECC(m)=ldeck/2-l_marciapiede-sum(length_corsie(1:m))-length_corsie(m+1)/2;
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
        sum(MATRIXloadD{i}(2,:))*-dist_appoggi(2)];
    else
            momtempy=0;
    end
    
    MsyACC(i)=abs(sum(momtempy));
end

clear momtempx momtempy

%% AZIONE FRENANTE
%Si calcola, la forza di frenatura e il momento longitudinale prodotto dalla stessa;
%si assume agente nella direzione dell’asse della strada ed al livello della
%superficie stradale, d'intensità pari  1/10 del sovraccarico costituito dallo schema 1.

q1fren=matrixEQLOAD(1);
if strcmp(bridgetype, 'simplysupp')
    for i=1:length(dist_appoggi)
        if strcmp(appoggio_testapila{i},'fixed')
            Fx_trafT(i)=(1/10*(q1fren*(tributL(i)*2)).*fi_d); %kN, forza di frenatura
        elseif strcmp(appoggio_testapila{i},'free')
            Fx_trafT(i)=0;
            My_traffic(i)=0;
        else
            warning('error in appoggio_testapila')
        end
    end
else
    Fx_trafT=(1/10*(q1fren*sum(tributL)).*fi_d); %kN, forza di frenatura
end
Fx_traffic=sum(Fx_trafT);
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
    Fy_windT(2)=sum(Fy_wind_deck)+Fy_wind_bent+Fy_wind_acc(1); %kN, carico accidentale su prima campata
    Fy_windT(3)=sum(Fy_wind_deck)+Fy_wind_bent+Fy_wind_acc(2); %kN  carico accidentale su prima campata
end
Fy_wind=[Fy_windT,...
    sum(Fy_wind_deck)+Fy_wind_bent]; %kN  carico accidentale non presente

%Le forze del vento che agiscono sull'impalcato producono momenti flettenti
%trasversali alla base della pila.
Mx_windT(1)=sum(Fy_wind_deck)*(h_bent+(h_slab+h_girder)/2)+sum(Fy_wind_acc)*(h_bent+h_slab+h_girder+1.5)+Fy_wind_bent*h_bent/2;%kNm
if strcmp(bridgetype,'simplysupp')
    Mx_windT(2)=sum(Fy_wind_deck)*(h_bent+(h_slab+h_girder)/2)+sum(Fy_wind_acc)*(h_bent+h_slab+h_girder+1.5)+Fy_wind_bent*h_bent/2;%kNm
    Mx_windT(3)=sum(Fy_wind_deck)*(h_bent+(h_slab+h_girder)/2)+(Fy_wind_acc(1))*(h_bent+h_slab+h_girder+1.5)+Fy_wind_bent*h_bent/2;%kNm
end
Mx_wind=[Mx_windT,...
    (sum(Fy_wind_deck)*(h_bent+(h_slab+h_girder)/2)+Fy_wind_bent*h_bent/2)];%kN  carico accidentale non presente

clear Fy_windT My_windT


end