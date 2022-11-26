% Nom du fichier d'output a analyser (modifiez selon vos besoins)
filename = 'output.out'; 

% Chargement des donnees
data = load(filename);

% Extraction des quantites d'interet
% (Le code c++ ecrit t, x(t), v(t), E_mec(t)  ligne par ligne, 
%  une ligne par pas de temps)
t    = data(:,1);
mF   = data(:,2);
x    = data(:,3);
v    = data(:,4);
Emec = data(:,5);
P    = data(:,6); 

% nombre de pas de temps effectués:
nsteps = length(t);
% longueur du pas de temps:
dt = t(2)-t(1);

% Figures
% line width and font size (utile pour la lisibilité des figures dans le
% rapport)
lw=2; fs=16; 
figure('Name', [filename ': x(t)'])
plot(t, x, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t [s]')
ylabel('x [m]')
grid on

figure('Name', [filename ': v(t)'])
plot(t, v, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t [s]')
ylabel('v [m/s]')
grid on

figure('Name', [filename ': (x,v)'])
plot(x, v, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('x [m]')
ylabel('v [m/s]')
grid on

figure('Name', [filename ': Emec(t)'])
plot(t, Emec, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t [s]')
ylabel('Emec [J]')
grid on

%calcul de la dérivée de l'énergie mécanique
dEmec = (Emec(3:end) - Emec(1:end-2))*(1/(2*dt));
figure ()

plot(t(2:end-1),dEmec,'r--','LineWidth',3.5)
hold on
set(gca,'fontsize',fs)
xlabel('t [s]')
ylabel('P_{mec}')
plot(t(2:end-1),P(2:end-1),'Color','black','linewidth',1.25)
hold off


figure()
plot(t(2:end-1),abs(dEmec- P(2:end-1)),'LineWidth',2)
set(gca,'fontsize',fs)
xlabel('t [s]')
ylabel('erreur sur P_{mec}')


