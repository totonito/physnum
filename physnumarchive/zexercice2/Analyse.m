close all
% Nom du fichier d'output a analyser (TODO: modifier selon vos besoins)
filename = 'output_example.out';

% Chargement des donnees
output = load(filename);

% Extraction des quantites d'interet 
% TODO: verifier la consistance avec l'ecriture du fichier output par le code  C++
t  = output(:,1);
x  = output(:,2);   
y  = output(:,3);
z  = output(:,4);
vx = output(:,5);
vy = output(:,6);
vz = output(:,7);
mu = output(:,8);
E  = output(:,9);
vpn = output(:,10);
vpx = output(:,11);
vpy = output(:,12);
vpz = output(:,13);

% Figures
fs=16; % font size
lw=2;  % linewidth
% 
figure
% 
% % Si on veut faire une figure combinée de 4 sous-figures, 
% % enlever le commentaire des lignes %subplot(...)
% % et commenter les lignes 'figure' CI-APRES
% % subplot(2,2,1) % array 2x2 of subfigures, 1st subfigure
% %subplot(2,2,1)
plot(x,y,'linewidth',lw)
set(gca,'fontsize',fs)
axis equal
grid on
xlabel('x [m]')
ylabel('y [m]')

% 
% 
figure
plot3(x,y,z,'linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

figure
%subplot(2,2,2) % array 2x2 of subfigures, 2nd subfigure
plot(vx,vy,'linewidth',lw)
set(gca,'fontsize',fs)
axis equal
grid on
xlabel('v_x [m/s]')
ylabel('v_y [m/s]')
% 
figure
%subplot(2,2,3) % array 2x2 of subfigures, 3rd subfigure
plot(t,x,t,y,'linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('t [s]')
ylabel('x,y [m]')
legend('x','y')
% 
% 
figure
%subplot(2,2,3) % array 2x2 of subfigures, 3rd subfigure
plot(t,z,'linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('t [s]')
ylabel('z[m]')
% 
figure
%subplot(2,2,4) % array 2x2 of subfigures, 4th subfigure
plot(t,vx,t,vy,'linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('t [s]')
ylabel('v_x,v_y [m/s]')
legend('v_x','v_y')
% 
figure
%subplot(2,2,4) % array 2x2 of subfigures, 4th subfigure
plot(t,vz,'linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('t [s]')
ylabel('v_z [m/s]')

%figure

%figure
subplot(1,2,1)
plot(t,vpn,'linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('t [s]')
ylabel('v_{para} [m/s]')
%figure
subplot(1,2,2)
plot3(vpx,vpy,vpz,'linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('v_{para x} [m/s]')
ylabel('v_{para y} [m/s]')
zlabel('v_{para z} [m/s]')

Emec_init = 2.09075*10^-16; % A MODIFIER SELON LA SITUATION !!


subplot(1,2,1)
plot(t,abs(E-Emec_init),'linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('t [s]')
ylabel('Erreur sur E_{mec} [J]')

 tfin = 8.e-08;                 %A MODIFIER SI BESOIN EST !
 dt = tfin ./ 1000;%

%figure
dE = (E(3:end) - E(1:end-2))*(1/2*dt);
subplot(1,2,2)
plot(t(2:end-1),dE,'k-','LineWidth',2)
set(gca,'fontsize',fs)
xlabel('t [s]')
ylabel('dE_{mec}[J/s]')
grid on

figure

plot(t,mu,'linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('t [s]')
ylabel('Moment magnétique[A\cdot m^2]')

