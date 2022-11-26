% figures pour une simulation

% TODO: changer le nom du fichier de sortie pour qu'il corresponde à celui
% que vous voulez analyser
output = load('output.out');

% TODO: les variables et indices 1 à 10 doivent correspondre à ce que votre code C++
% écrit dans le fichier de sortie
t   = output(:,1);
xA  = output(:,2);
yA  = output(:,3);
vxA = output(:,4);
vyA = output(:,5);
Emec= output(:,6);
Pnc = output(:,7);
rho = output(:,8);
acc = output(:,9);
dt  = output(:,10);

clear output

RT=6.3781e6;
theta = linspace(0,2*pi,200);

lw=1; fs=16; % TODO: adapter selon ce que vous voulez

% trajectoire d'Apollo 13 et surface de la terre
figure
plot(xA,yA,'b-', 'linewidth',lw)
hold on
plot(RT*cos(theta),RT*sin(theta),'r-')
axis equal
set(gca,'fontsize',fs)
xlabel('x[m]')
ylabel('y[m]')
grid on
hold off

% TODO: faire d'autres figures selon vos besoins

