function [trasvbar]=TreinfdesignMTAcirc(year,D,Ved,tau_adm,sigmas_adm,passo_staffe,varargin)

%Progetto a taglio una sezione rettangolare con MTA
%% Optional input

% Maximum number of optional inputs
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 1 optional inputs');
end

% set defaults for optional inputs
optargs = {0.07};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
cover = optargs{:};

%% INPUT
%se tau_max< tau_adm1,non è necessaria una progettazione a taglio
%basta adottare un'armatura a taglio minima.
h= D-cover; %m, altezza utile
fi_staffe=[6 8 10 12 14 16 18 20]*10^-3;% [m] tipologia staffe
area_staffe=(pi.*(fi_staffe/2).^2); %mq, area sezione della staffa
geomdim=[D, D;
    D, D];

%Calcolo tensione tangenziale max
%Uso 0.7 d per passare a sezione rettangolare equivalente
tau_max=Ved'./((0.7*(geomdim(:,1)-cover)).*(geomdim(:,2)));%kN/mq

%% Progetto a taglio

%%%%% BRACCI
%Definisco 1.5 m come limite empirico, oltre il quale non si possono usare
% solo due bracci, ma devo incrementare. Fisso il numero di bracci
nb=[2 2]; %Sezione rettangolare equivalente a due bracci. Il prof ha detto che se la tau max era maggiore
%della tau minima si incrementava il diametro della sezione

% for i=1:2
%     if geomdim(i,2)>1.5 %m
%         nb(i)=2+(ceil( geomdim(i,2)/1.5)-1);
%     elseif  geomdim(i,2)<=1.5 %m
%         nb(i)=2; %numero bracci minimo
%     end
% end

for i=1:2
    if tau_max(i)<=tau_adm(1)
        
%         Ast_min=0.1*geomdim(i,2)*10^(-2)*(1+0.15*geomdim(i,1)/geomdim(i,2)); %mq/m, armatura minima delle staffe per metro
%         
%         Ast_y=(area_staffe.*nb(i)); %mq, area staffa parallela asse y comprensiva di n bracci
%         
%         Ast_y1=Ast_y'.*(floor(1./passo_staffe));%mq/m armatura delle staffe per metro
%         f_temp=sort(Ast_y1(:));%  variabile temporanea che rappresenta in ordine crescente
%         %gli elementi del vettore Ast_y1
%         
%         %trovo l'area di staffe minima per rispettare i limiti normativi
%         index=find(f_temp(:)>Ast_min,1,'first');
%         
%         Ast_y_def=f_temp(index);% mq/m, armatura delle staffe per metro
%         [x,y]=find(Ast_y1==Ast_y_def);%è neccessario trovare gli indici dell'elemento
%         %della matrice Ast_x1 che verificano la condizione Ast_y1==Ast_y_def

        %%%%% 1) Uso i Limiti normativi   
        fi_staffe_def(i)=0.006; %m tipologia staffa di progetto
        passo_staffe_def(i)=0.25; %m,  passo staffe di progetto
        
    elseif tau_max(1)>tau_adm(1) && tau_max(1)<=tau_adm(2)
        
        %%%%% progettazione a taglio
        
        Ast=(Ved(i)*passo_staffe)/((geomdim(i,1)-cover)*sigmas_adm); %mq, area minima della singola staffa comprensiva di n bracci       
        Ast_y=(area_staffe.*nb(i)); %cmq, area staffa parallela asse y comprensiva di n bracci
        
        %h_temp è la variabile temporanea ottenuta dal rapporto delle due armature,
        %quando il rapporto è maggiore o uguale ad 1, vorrà dire che l'area staffa
        %è maggiore dell'area minima
        h_temp=Ast_y'./Ast;
        m_temp=sort(h_temp(:));
        ind=find(m_temp>=1,1,'first');
        [min_val]=m_temp(ind);
        
        %è neccessario trovare gli indici dell'elemento
        %della matrice h_temp che verificano la condizione rapp_minimo==h_temp)
        [x,y]=find(min_val==h_temp);
        
        fi_staffe_def(i)=fi_staffe(x);% tipologia staffa di progetto
        passo_staffe_def(i)=passo_staffe(y);%m, passo staffe di progetto
        
        
        
    elseif tau_max> tau_adm(2)
        tau_max=NaN;
        disp 'modificare le caratteristiche geometriche della sezione in modo da portare tau_max<tau_adm1'
        fi_staffe_def=0;% tipologia staffa di progetto
        passo_staffe_def=0;%cm, passo staffe di progetto
        nb=0;
    end
    
end

%Take the minumum needed Ast/s. (per sezioni circolari dovrebbe venire uguale)

%Calcolo il taglio in X con la config di As/s scelta per la direzione y
Vr_temp(1)=pi.*(fi_staffe_def(2)/2).^2./passo_staffe_def(2)*nb(1)*0.9*(geomdim(1)-cover)*sigmas_adm;
%Calcolo il taglio in Y con la config di As/s scelta per la direzione x
Vr_temp(2)=pi.*(fi_staffe_def(1)/2).^2./passo_staffe_def(1)*nb(2)*0.9*(geomdim(2)-cover)*sigmas_adm;
%Scelgo la configurazione che mi da rapporto Vr/Ve > 1 nelle due direzioni
ind=find(Vr_temp./Ved>1);
if isempty(ind)
    ind=1; %Avoid empty ind (when there are the minimum in both directions)
end
%% OUTPUT
trasvbar=[fi_staffe_def(ind), passo_staffe_def(ind), nb];
end
