function [F1, F2]=TLM(area, no_vibration)
        
        
%%%%%%%%%%%%%%%  Calcul de la fonction de transfert et des formantscc
%%%%%%%%%%%%%%%  Pierre BADIN
        
        ORAL.L = area(:,1)';
        ORAL.A = area(:,2)';
        
        nbfreq = 500; Fmax = 10000; Fmin = Fmax / (nbfreq - 1);  
        
        f = linspace(Fmin, Fmax, nbfreq);

        [H_eq, F_form] = vtn2frm_ftr_oral(ORAL, nbfreq, Fmax, Fmin, no_vibration);
        
        spec=abs(H_eq);

        %------------------------------------------------------------------------------
        % Nettoyage des eventuels poles/zeros redondants, mise dans l'ordre et nettoyage
        % des poles/zeros au dela de Fmax (AS le 28/02/2006)
        seuil2 = 10;
         % Poles
    if ~isempty(F_form)
        % Mise en ordre
        [bid, ind_F_ord] = sort(F_form(:,2)); F_form = F_form(ind_F_ord,:);
        % Nettoyage des poles superieurs a Fmax et negatifs
        F_form = F_form(find((F_form(:,2) <= Fmax) & (F_form(:,2) >= 0)),:);
        % Nettoyage des p?les redondants
        ind_F_diffOK = find(diff(F_form(:,2)) >= seuil2);
        if ~isempty(ind_F_diffOK)
            ind_F_Nredondant = [ind_F_diffOK; ind_F_diffOK(end) + 1];
        elseif isempty(ind_F_diffOK) & ~isempty(F_form) % if ~isempty(ind_F_diffOK)
            ind_F_Nredondant = 1;
        else  % if ~isempty(ind_F_diffOK)
            ind_F_Nredondant = [];
        end  % if ~isempty(ind_F_diffOK)
        F_form = F_form(ind_F_Nredondant,:);
    end  % if ~isempty(F_form)
    if isempty(F_form); F_form = NaN * ones(1,4); end 



    %---------------------------------------------------
    % Calcul des valeurs de la FT aux zeros et aux poles (AS le 27/02/2006)
    for i_pole = 1:size(F_form,1)
        F_form(i_pole,3) = 20*log10(abs(aire2spectre_oral(ORAL, 1, F_form(i_pole,2), F_form(i_pole,2), no_vibration)));
    end  
        
%       F1=round(F_form(1,2));
%       F2=round(F_form(2,2));
%       F3=round(F_form(3,2));
      F1=F_form(1,2);
      F2=F_form(2,2);
      % F3=F_form(3,2);
return
 
%---c----C----c----C----c----C----c----C----c----C----c----C----c----C
%
%                     INSTITUT DE LA COMMUNICATION PARLEE
%                               GRENOBLE, FRANCE
%
%                               ** nraph3.ftn **
%
%       Creation : 7 Juin 1987
%
%
% Hugo SANCHEZ
% Creation : 19 Juin 1984
% Mis a jour : 19 Juin 1984
% Derniere modification : Pierre BADIN (21 novembre 1985)
%     "         "       : 86-08-29 (Possibilite de recherche de zeros
%                       et de poles simples et complexes conjugues,
%                       aussi bien dans TP que dans TZ.
%     "         "       : 87-02-13 (Commentaires pour une version globale 0)
% Derniere Modification : 87-03-04 (Passage sur VAX)
%     "         "       : 87-06-08 (Conditions sur INAS)
%-----------------------------------------------------------------------------
%       Entree
%               FE, BNPE, Estimation de depart des parametres
%               ITERMX, Nombre d'iterations autorisees
%               MODE =1 recherche des poles de TP
%                    =2     "         poles de TZ
%                    =3     "     des zeros de TP
%                    =4     "         zeros de TZ
%               FMAX, Frequence maximum autorisee
%
%       Sorties
%               F, BP, Parametres estimes par le programme
%               NRAPH, Nombre d'iterations (=0 si pas de solution)


function [F, BP, ITER, FRQ, HHH] = nraph_oral(FE, BNPE, ITERMX, FMAX, ORAL, F_form, no_vibration);

DELTAS = complex(30, 30); SEUIL = 0.3; ITER = 1;


SI = complex(-BNPE*pi, FE*2*pi);

F = NaN; BP = NaN;

% Suivant que Z_form existe ou non, on cherche les pôles ou les zéros

% Recherche des poles ou des zéros
while ITER < ITERMX
	% Calcul du premier point
	FIcx = -j * SI / (2*pi); 
	Q = 1 / aire2spectre_cor_oral(ORAL, 1, FIcx, FIcx, F_form, no_vibration);

	% Calcul du deuxieme point
	SIP = SI + DELTAS; FIPcx = -j * SIP / (2*pi);
	QP = 1 / aire2spectre_cor_oral(ORAL, 1, FIPcx, FIPcx, F_form, no_vibration);
	% Point d'intersection avec 0
	Q1D = (QP - Q) / DELTAS; SISU = SI - (Q / Q1D);
	
	% Si on a trouve ...
	if(abs(SISU-SI) < SEUIL) F = imag(SISU) / (2*pi); BP = - real(SISU)/pi; return
	% ... sinon
	else SI = SISU; ITER = ITER + 1;	end
end  % while ITER < ITERMX	


return

function [H_eq, F_form] = vtn2frm_ftr_oral(ORAL, nbfreq, Fmax, Fmin, no_vibration);

% Auteur: Pierre Badin Gipsa-Lab

f = linspace(Fmin, Fmax, nbfreq);

FE = 150; BNPE = 50; % Valeurs initiales de recherche
FINC = 100; % Incrément pour initialiser le point suivant

ITERMX = 100;

seuil2 = 10; % Différence minimale entre deux formants

FMAX = Fmax; 


F_form = []; 
% Détermination des pics du module de fonction de transfert initiale
H_eq = aire2spectre_cor_oral(ORAL, nbfreq, Fmax, Fmin, F_form, no_vibration);
[maxi, ind_maxi, mini, ind_mini] = peakpick(20*log10(abs(H_eq)));
frq_max_first = f(ind_maxi); 

	
% -------------------------------------------------------------------------------------
NF = 0; F = 0;
while F <= FMAX
	if NF > 0 FE = F + FINC; end

	% Détermination des pics du module de fonction de transfert
	H_eq = aire2spectre_cor_oral(ORAL, nbfreq, Fmax, Fmin, F_form, no_vibration);
	[maxi, ind_maxi, mini, ind_mini] = peakpick(20*log10(abs(H_eq)));
	frq_min = f(ind_mini); frq_max = f(ind_maxi); 
	FEST = frq_max; 
	FE = FEST(1);
	
	% Recherche des poles de H_eq par nraph
  [F, BP] = nraph_oral(FE, BNPE, ITERMX, FMAX, ORAL, F_form, no_vibration);
	NF = NF + 1; F_form(NF, 2) = F; F_form(NF, 4) = BP;
% 	round(F_form)
end % while F <= FMAX

% On supprime dans frq_max_first les formants proches de ceux trouvés ci-dessus
seuil3 = 20;
for ifirst = 1:length(frq_max_first)
	ind = find(abs(frq_max_first(ifirst) - F_form(:, 2)) < seuil3);
	if ~isempty(ind) frq_max_first(ifirst) = NaN; end
end 
frq_max = frq_max_first(find(isfinite(frq_max_first)));

% Boucle sur tous les pics trouvés
for ifirst = 1:length(frq_max)
	FE = frq_max(ifirst);
	% Recherche des poles de H_eq par nraph
  [F, BP] = nraph_oral(FE, BNPE, ITERMX, FMAX, ORAL, F_form, no_vibration);
	NF = NF + 1; F_form(NF, 2) = F; F_form(NF, 4) = BP;
end % for NF = 1:length(frq_max)

% On nettoie les éléments qui débordent ...
ind_KO = find((F_form(:, 2) < 0) | F_form(:, 2) > FMAX); F_form(ind_KO, :) = [];

% Pôles
if ~isempty(F_form)
	% Mise en ordre
	[bid, ind_F_ord] = sort(F_form(:,2)); F_form = F_form(ind_F_ord,:);
	% Nettoyage des pôles supérieurs à Fmax et négatifs!!!
	F_form = F_form(find((F_form(:,2) <= Fmax) & (F_form(:,2) >= 0)),:);
	% Nettoyage des pôles redondants
	ind_F_diffOK = find(diff(F_form(:,2)) >= seuil2);
	if ~isempty(ind_F_diffOK)
		ind_F_Nredondant = [ind_F_diffOK; ind_F_diffOK(end) + 1];
	elseif isempty(ind_F_diffOK) & ~isempty(F_form) % if ~isempty(ind_F_diffOK)
		ind_F_Nredondant = 1;
	else  % if ~isempty(ind_F_diffOK)
		ind_F_Nredondant = [];
	end  % if ~isempty(ind_F_diffOK)
	F_form = F_form(ind_F_Nredondant,:);
end  % if ~isempty(F_form)
if isempty(F_form); F_form = NaN * ones(1,4); end

% Calcul de la fonction d'aire complète
H_eq = aire2spectre_oral(ORAL, nbfreq, Fmax, Fmin, no_vibration);

return

function [H,YE]=spectrelec(w,A,zr,l,no_vibr_paroi);

% Auteur: Pierre Badin GIPSA-Lab 

% Si rien n'est prŽcisŽ, on prend bien en compte
% les vibrations de paroi
if nargin < 5
    no_vibr_paroi = 0;
end  % if nargin < 5

% utilisation des variables globales
global CONST_DAT

% Par defaut, toutes les longueurs sont egales a 1 cm
if (nargin < 4 )
   l = ones(size(A));
end;
 
% Recuperation des variables
%---------------------------
c = CONST_DAT(1);
ro = CONST_DAT(2);
lambda = CONST_DAT(3);
eta = CONST_DAT(4);
mu = CONST_DAT(5);
cp = CONST_DAT(6);
bp = CONST_DAT(7);
mp = CONST_DAT(8);

% Determination de L,C,R,G,et YP
%-------------------------------

% PeŽrime?tre du tube ŽlŽmentaire
S = 2*sqrt(A*pi) ;
% Caracteristiques du tube
L = ro./A.*l ;
C = A.*l/ro/c/c; % (A*l)/(ro*c^2)

% Perte de Fant 1960
R_coef = sqrt(ro*mu/2*w) ;
G_coef = (eta-1)/(ro*c^2)*sqrt(lambda*w/(2*cp*ro));
R = S.*l./(A.^2) * R_coef ;
G = S.*l * G_coef ;
% pertes par vibration et viscosite le long des parois
YP_coef = 1./(bp^2+mp^2*w.^2) ;
YP = S.*l * ( (bp-j*mp*w).*YP_coef )  ;
if no_vibr_paroi;
    YP = 0 ;
end;

% Determination de Y, Z
%----------------------
Z = R+j*L*w ;
Y = G+j*C*w+YP ;

% Liste des matrices de quadripole
%---------------------------------
aa = (1 + (Z.*Y/2));
bb = - (Z + Z.^2 .* Y/4);
cc = - Y;
dd = aa;
% calcul de la fonction de transfert
%-----------------------------------
% H=[aa bb ;
%    cc dd]

%    | A  B |
%H = |      |
%    | C  D |
%
%    | a  b |    | A  B |
%H = |      | * |      |
%    | c  d |    | C  D |

aaa = aa(1,:) ;
bbb = bb(1,:) ;
ccc = cc(1,:) ;
ddd = dd(1,:) ;
for ind = (1:length(A)-1)
    proda = aa(ind+1,:).*aaa + bb(ind+1,:).*ccc ;
    prodb = aa(ind+1,:).*bbb + bb(ind+1,:).*ddd ;
    prodc = cc(ind+1,:).*aaa + dd(ind+1,:).*ccc ;
    prodd = cc(ind+1,:).*bbb + dd(ind+1,:).*ddd ;
    aaa = proda ;
    bbb = prodb ;
    ccc = prodc ;
    ddd = prodd ;
end

if (zr==inf)
    H=0;
    YE= -ccc./ddd ;
else
    % calcul de la fonction de transfert
    % pl = aaa * pg + bbb * ug
    % ul = ccc * pg + ddd * ug
    % or pl = zr * ul
    % on obtient donc, (aaa*ddd-bbb*ccc =1)
    % ul / ug = 1 / (aaa - ccc * zr)
    H = ones(size(aaa))./( aaa - ccc.*zr) ;
    % calcul de la conductance d'entree
    YE = -( aaa - ccc.*zr ) ./ ( bbb - ddd.*zr ) ;
end;

function [lg_out,aire_out] = lgvar2lgfix(lg_in, aire_in, nb_sections);

% Pierre Badin Gipsa-Lab
% permet de surechantillonner une fonction d'aire initiale en nb_sections
% Calcul de la longueur du conduit vocal d'origine
% la fonction d'aire est en colonne (aire, longueur);

DL2 = sum(lg_in)/nb_sections;
[nb_sect_in bid] = size(lg_in);


if bid > 1
    lg_in = lg_in';
    aire_in = aire_in';
    [nb_sect_in bid] = size(lg_in);
end


lg_out = DL2*ones(1, nb_sections);
aire_out = zeros(1, nb_sections);

% Calcul effectif
I1 = 1;
XC = 0;
XFIN1 = lg_in(I1);

% Boucle sur les trames de sortie
for I2 = 1 : nb_sections,
  XFIN2 = I2*DL2;
  Vol = 0;
  while (XFIN1 < XFIN2 - 1.0e-10) & (I1 <= nb_sect_in)
    Vol = Vol + aire_in(I1)*(XFIN1-XC);
    I1 = I1 + 1;
    XC = XFIN1;
    XFIN1 = XC + lg_in(I1);
  end
    Vol = Vol + aire_in(I1)*(XFIN2-XC);
    XC = XFIN2;
    aire_out(I2) = Vol/DL2;

end

% [maxi, ind_maxi, mini, ind_mini] = peakpick(signal);
% DŽtection des maximas et minimas locaux d'un signal
% Auteur: Pierre Badin Gipsa-Lab

% EntrŽe:
%   signal(1:nb_ech) : signal ˆ analyser
%
% Sorties:
%   maxi, ind_maxi   : valeurs et indices des maximas locaux
%   mini, ind_mini   : valeurs et indices des minimas locaux
%

function [maxi, ind_maxi, mini, ind_mini] = peakpick(signal);

% DeŽtection des minimas et maximas locaux d'un signal
mini = []; ind_mini = []; maxi = []; ind_maxi = [];

nb_ech = length(signal);
for ind = 2:nb_ech-1
  if (signal(ind-1) < signal(ind)) & (signal(ind+1) < signal(ind))
    ind_maxi = [ind_maxi, ind]; maxi = [maxi, signal(ind)];
  end
  if (signal(ind-1) > signal(ind)) & (signal(ind+1) > signal(ind))
    ind_mini = [ind_mini, ind]; mini = [mini, signal(ind)];
  end    
  
end

%  Si on n'a pas trouveŽ d'extremum local, il faut prendre l'extre?miteŽ du  signal
if isempty(maxi)
    [maxi, ind_maxi] = max(signal);
end
if isempty(mini)
    [mini, ind_mini] = min(signal);
end
return

function H_eq = aire2spectre_cor_oral(ORAL, nbfreq, Fmax, Fmin, F_form, no_vibration);

% Auteur: Pierre Badin Gipsa-Lab

H_eq = aire2spectre_oral(ORAL, nbfreq, Fmax, Fmin, no_vibration);

f = linspace(Fmin, Fmax, nbfreq); 
if isreal(f) SI = complex(0, 2 * pi * f); else SI = 2 * pi * j * f; end

NF = size(F_form, 1);
for I = 1 : NF
	SI1 = complex(F_form(I, 4) * pi, F_form(I, 2) * 2 * pi);
	H_eq = H_eq .* ((SI - SI1) .* (SI-conj(SI1))) ./ (SI1 .* conj(SI1));
end

return

function [H_oral, Y_oral] = aire2spectre_oral(ORAL, nbfreq, Fmax, Fmin, no_vibration);

% Auteur: Pierre Badin Gipsa-Lab
% function [H_oral, Y_oral] = aire2spectre_oral(ORAL{, nbfreq, Fmax, Fmin});
%
% Calcule:
%   - la FT totale pour les freq Fmin à Fmax à partir de la fonction d'aire ORAL et NASAL
%   - Les formants associés
%   - la FT orale
%
% Entrées
%      ORAL.A(nb_tubes_oral)	  : Aires des tubes oraux (glotte -> lèvres)
%      ORAL.L(nb_tubes_oral)	  : Longueurs des tubes oraux (glotte -> lèvres)
%
%      nbfreq ([256])           : Nombre de point d'échantillonnage fréquentiel entre 0 et Fmax
%      Fmax ([5000])            : Fréquence maximale de la fonction de transfert
%      Fmin ([Fmax/nbfreq])     : Fréquence minimale de la fonction de transfert
%
%
% Sorties
%      H_oral(nbfreq)           : Fonction de transfert orale (linéaire)
%      Y_oral(nbfreq)           : Impédance d'entrée du tuyau oral

global CONST_DAT

if nargin == 2 nbfreq = 256; Fmax = 5000; end;
if nargin == 3 Fmax = 5000; end;

% valeurs des constantes Mrayati
c = 35100; % cm/s
rho = 1.14e-3; % g/cm^3  ; Badin & Fant 

%rho = 1.18e-3; % g/cm^3  ; Helie Thomas
lambda = 5.5e-5;
eta = 1.4;
mu = 1.86e-4;
cp = 0.24;
bp = 1600; % dynes.s/cm
mp = 1.4; % g/cm²

CONST_DAT = [c ; rho ; lambda ; eta ; mu ; cp ; bp ; mp];

% Variables frequences et pulsations
%%% ATTENTION : Pas de Fréquence Nulle !!!
f = linspace(Fmax/nbfreq, Fmax, nbfreq);
if nargin == 5 f = linspace(Fmin, Fmax, nbfreq); end  % if nargin == 5
w=2*pi*f;


% Fonction de transfert Bucale (= de l'embranchement avec le conduit nasal aux lèvres)
%=====================================================================================
% Fonction de Transfert
% Impédance de rayonnement aux lèvres
Zr_oral =  rho/(2*pi*c)*(w.^2)+j*8*rho/(3*pi.*sqrt(pi*ORAL.A(end)))*w; 
% Zr_oral =  Inf * ones(size(f)); 
% Zr_oral =  0 * ones(size(f)); 
[H_oral, Y_oral] = spectrelec(w, ORAL.A', Zr_oral, ORAL.L',no_vibration); % 0 = vibration des parois

% %***************************************
% figure
% clf; hold on; grid on;
% plot(f, 20*log10(abs(H_oral)))
% %***************************************
return









