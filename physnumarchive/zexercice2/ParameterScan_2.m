d% Ce script Matlab automatise la production de résultats
% lorsqu'on doit faire une série de simulations en
% variant un des paramètres d'entrée.
% 
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un paramètre du fichier d'input
% par la valeur désirée.
%

%% Parametres %%
%%%%%%%%%%%%%%%%

% MODIFIER SELON VOS BESOINS LES NOMS ET LES VALEURS
 % Chemin d'accès au code compilé (NB: enlever le ./ sous Windows)
 repertoire = '';
executable = 'Exercice2_Comerre_Vion.exe'; % Nom de l'exécutable (NB: l'extension.exe est nécessaire sous Windows)
input = 'configuration.in.example'; % Nom du fichier d'entrée 

nsteps = [500 600 700 800 900 1000 1100 1200 1300 1400 1600 2000 2600];
nsimul = numel(nsteps); % Nombre de simulations a faire:
% numel est une fonction Matlab qui retourne le nombre d'elements d'un tableau.

% Autre exemple: 
%nsteps = round(logspace(2,4,nsimul)); % Nombre d'iterations entier de 10^2 a 10^4
% logspace est une fonction Matlab retournant un tableau de valeurs dont les logarithmes sont equidistants
% tapez 'help logspace' pour plus de details
% Voir aussi la fonction 'linspace'

tfin = 8.e-08; % TODO: Verifier que la valeur de tfin est EXACTEMENT la meme que dans le fichier input

dt = tfin ./ nsteps;

paramstr = 'nsteps'; % Nom du parametre a scanner  MODIFIER SELON VOS BESOINS
param = nsteps;  % Valeurs du parametre a scanner  MODIFIER SELON VOS BESOINS

%% Simulations %% 
%%%%%%%%%%%%%%%%%
% Lance une serie de simulations (= executions du code C++)
% Normalement, on ne devrait pas avoir besoin de modifier cette partie

output = cell(1, nsimul); % Tableau de cellules contenant les noms des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i)), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
    disp(cmd)
    system(cmd);
    disp('Done.')
    
end

%% Analyse %%
%%%%%%%%%%%%%
% Ici, on aimerait faire une etude de convergence: erreur fonction de dt, sur diagramme log-log.
% A MODIFIER ET COMPLETER SELON VOS BESOINS
q = 1.6022e-19; %charge of the proton [C]
m = 1.6726e-27; %mass of the proton [kg]
B = 4;          %intensity of the magnetic field [T]
omega = q*B/m;  %cyclotron frequency [rad/s]
v = 5e5;        %proton initial velocity modulus [m/s]
r=v/omega;      %radius of the trajectory [m]

error = zeros(1,nsimul);
norme_position_finale = [];
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data  = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    t     = data(:,1); 
    vx    = data(end,5);
    vy    = data(end,6);
    x = data(end,2);
    y = data(end,3);
    z = data(end,4);
% TODO:  inserer ici les expressions de la solution exacte
    x_th  = (1-cos(omega * tfin))*v/omega; 
    y_th  = sin(omega * tfin)*v/omega; 
    vx_th = v*sin (omega*tfin); 
    vy_th = v * cos(omega * tfin); 
    tampon = [sqrt( x^2 + y^2 + z^2)];
    norme_position_finale = [norme_position_finale tampon];
    
    error(i) = max(abs(vx-vx_th),abs(vy-vy_th));
end

pente_x = [ 100 1000 10000];
pente_2_y = [10000 100 1];
pente_1_y = [1000000 100000 10000];
lw=1.5; fs=16;
figure
loglog(nsteps, error, 'r+-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('N_{steps}')
ylabel('Err v finale [m/s]')
hold on
loglog(pente_x,pente_1_y,'k--','linewidth',lw)
loglog(pente_x,pente_2_y, 'k--','linewidth',lw)
axis( [ 1000 10000 100000 1000000])
grid on
legend('BorisBuneman')
% %legend ('Euler explicite')
% %legend ('Euler semi-implicite')
% legend ('Euler implicite')
% %legend ('Runge Kutta2')
hold off

% Si on n'a pas la solution analytique: on représente la quantite voulue
% (ci-dessous v_y, modifier selon vos besoins)
% en fonction de (Delta t)^norder, ou norder est un entier.


norder=2; % TODO: modifier si besoin est

% 
figure
plot(dt.^norder,norme_position_finale,'k+-','linewidth',lw)
set(gca,'fontsize',fs)

xlabel('(\Delta t)^2 [s^2]')
ylabel('norme position finale[m]')
grid on



