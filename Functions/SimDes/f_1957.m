 function[nf_x_def,fi_x_def,fi_staffey_def,passo_staffe_def,nb_y_def,nf_y_def,fi_y_def,fi_staffex_def,nb_x_def]=f_1957(Lx_pier,Ly_pier,n,Nsd,Msd_x,Msd_y,Vsd_x,Vsd_y,sigmac_r,sigmas_adm)

%Progettazione a Pressoflessione e taglio della sezione rettangolare con MTA

%Convenzione assi: asse x è l'asse longitudinale dell'impalcato
%                  asse y è l'asse trasversale dell'impalcato

%% Elaborazione degli input
%Calcolo sigma ammissibile del cls
if sigmac_r>=250 %kg/cmq
    sigmac_adm= (75+(sigmac_r-225)/9); %kg/cmq,tensione ammissibile del conglomerato
elseif sigmac_r<250
    sigmac_adm=sigmac_r/3; %kg/cmq
end

a=Lx_pier*100;%cm
b=Ly_pier*100;%cm
copr= 7; %cm, copriferro
plotter='plot'; %write plot or char to plot results
%% Progetto dell'armatura longitudinale in direzione longitudinale
[As_y,a_tempy,nfy,Asy,nf_y_def,fi_y_def]=LreinfdesignMTA(a,b,n,Nsd,Msd_y,sigmac_adm,sigmas_adm,'l',plotter,copr);

%% Progetto dell'armatura longitudinale in direzione trasversale
[As_x,a_tempx,nfx,Asx,nf_x_def,fi_x_def]=LreinfdesignMTA(b,a,n,Nsd,Msd_x,sigmac_adm,sigmas_adm,'t',plotter,copr);

%% CHECK FINALE 
%calcolo tensione ammissibile ridotta del cls per sezioni soggette a
%compressione semplice
if min(a,b)>25
    sigmac_red= (0.7*sigmac_adm); %kg/cmq,tensione ridotta di compressione semplice
elseif min(a,b)<25
    sigmac_red = 0.7*(1-0.03*(25-min(a,b)))*sigmac_adm; %kg/cmq, tensione ridotta di compressione semplice
end

Ac_min=(max(Nsd))/sigmac_red; %cmq, area dicalcestruzzo strettamente necessaria per sezione 
%soggetta a N centrato.

%%%%Limiti normativi
%-L'armatura longitudinale deve essere maggiore del 0.8% di Ac_min se
%Ac_min<2000cmq
%-L'armatura longitudinale deve essere minore del 0.5% di Ac_min se
%Ac_min>8000cmq
%adottare per i casi intermedi l'interpolazione lineare
As_tot=As_x+As_y; %cmq,armatura longitudinale totale della sezione

if Ac_min<2000
    R=As_tot>=0.008*Ac_min; %Index for minimum reinforcement verification
elseif Ac_min>8000
    R=As_tot>=0.005*Ac_min;
elseif Ac_min>=2000 && Ac_min<=8000
    % Calcolo moltiplicatore Ac_min attraverso l'interpolazione lineare
    m=interp1([2000 8000],[0.008 0.005],Ac_min);
%     m=(0.008-0.005)/(8000-2000);%coefficinte angolare della retta
%     moltiplicatore=0.005+ m*(Ac_min-2000);
    R=As_tot>=m*Ac_min;
end

% OUTPUT As,long
%Se sono rispettati i minimi da normativa, allora trova gli output,
%altrimenti aumenta l'armatura in modo da verificare i minimi
if R==1
    % export data
    di_x_def=fi_x_def/10;%[cm]
    di_y_def=fi_y_def/10;%[cm]
    
elseif  R==0
    %Aumento prendendo la configurazione con As successiva, è una
    %condizione che si verifica molto difficilmente
    warning('limiti normativi non rispettati, aumento armatura')
    fi=[16 18 20 22 24 26 28];%[mm] tipolgia armatura della sezione
    
    [x2,~]=find(a_tempx==As_x);
    As_x2=a_tempx(x2+1);
    [x3,y3]=find(Asx==As_x2);
    nf_x_def=nfx(y3)*2;
    fi_x_def=fi(x3);
    di_x_def=fi(x3)/10;
    
    [x4,~]=find(a_tempy==As_y);
    As_y2=a_tempy(x4+1);
    [x5,y5]=find(Asy==As_y2);
    nf_y_def=nfy(y5)*2;
    fi_y_def=fi(x5);
    di_y_def=fi(x5)/10;    
end

%% Progetto a taglio, con forza di taglio 
% limiti min e max da norma
if sigmac_r<=120% valida  per  conglomerati di cemento idraulico normale(Portland), d'alto forno e pozzolanico.
    tau_adm=[4 14]; %kg/cmq, per conglomerai di cemento idraulico normale(Portland), d'alto forno e pozzolanico.
%     tau_adm2=14;%kg/cmq,massima tensione tangenziale ammissibile
elseif sigmac_r>=160 %valida per conglomerati di cemento ad alta resistenza ed alluminoso.
    tau_adm=[6 14]; %kg/cmq
%     tau_adm2=14;%kg/cmq,massima tensione tangenziale ammissibile.
end

%fisso tre passi di staffatura usuali dell'anno di costruzione, questa scelta verifica la condizione che ci siano almeno 3 staffe al metro
passo_staffe=[25 20 15];

% dir X
Ved_x=max(Vsd_x);%kg
%Calcolo tensione tangenziale max
tau_max=Ved_x/((0.9*a)*(b-copr));%kg/cmq

%Definisco i passi di tentativo
index=passo_staffe(:)>10*di_x_def | passo_staffe(:)>0.5*min(a,b);
passo_staffe(index)=0;

[fi_staffex_def,passo_staffex_def,nb_x_def]=TreinfdesignMTA(a,b,tau_max,tau_adm,Ved_x,sigmas_adm,passo_staffe);

%% Progetto a taglio,con forza di taglio lungo y

Ved_y=max(Vsd_y);%kg
tau_max=Ved_y/((0.9*b)*(a-copr));%kg/cmq

%Definisco i passi di tentativo
index=passo_staffe(:)>10*di_y_def | passo_staffe(:)>0.5*min(a,b);
passo_staffe(index)=0;

[fi_staffey_def,passo_staffey_def,nb_y_def]=TreinfdesignMTA(b,a,tau_max,tau_adm,Ved_y,sigmas_adm,passo_staffe);

passo_staffe_def=max(passo_staffey_def,passo_staffex_def);
 end
