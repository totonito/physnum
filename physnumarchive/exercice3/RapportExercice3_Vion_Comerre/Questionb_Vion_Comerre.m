filename="output_example.out";
output=load(filename);
t = output(:,1);
dt = output(2,1);
theta = output(:,2);
thetadot = output(:,3);
Emec = output(:,4);
Pnc = output(:,5);
dEmec = (Emec(3:end) - Emec(1:end-2))*(1/(2*dt));

lw=2;
fs=16;

% figure
% plot(t,theta, 'r', 'linewidth', lw)
% set(gca,'fontsize',fs)
% xlabel('$t$ [s]', 'Interpreter','latex', 'FontSize', 24)
% ylabel ('$\theta$ [$\mathrm{rad}$]', 'Interpreter','latex', 'FontSize', 24)
% axis([0 18 -2 2])
% 
% figure
% plot(t,thetadot, 'r', 'linewidth', lw)
% set(gca,'fontsize',fs)
% xlabel('$t$ [s]', 'Interpreter','latex', 'FontSize', 24)
% ylabel ('$\dot{\theta}$ [$\mathrm{rad}\cdot\mathrm{s}^{-1}$]', 'interpreter','latex', 'FontSize', 24)
% axis([0 18 -12 12])

% figure
% plot(theta,thetadot, 'r', 'linewidth', lw)
% set(gca,'fontsize',fs)
% xlabel('$\theta$ [$\mathrm{rad}$]', 'interpreter','latex', 'FontSize', 24)
% ylabel ('$\dot{\theta}$ [$\mathrm{rad}\cdot\mathrm{s}^{-1}$]', 'Interpreter','latex', 'FontSize', 24)
% axis([-1.75 1.90 -11 12])
% 
% figure
% plot(t,Emec, 'linewidth', lw)
% set(gca,'fontsize',fs)
% xlabel('$t$ [s]', 'Interpreter','latex', 'FontSize', 24)
% ylabel ('$E_{mec}$ [$\mathrm{J}$]', 'Interpreter','latex', 'FontSize', 24)
% axis([0 18 0 0.3])
% 
% figure
% plot(t(2:end-1),dEmec-Pnc(2:end-1), 'linewidth', lw)
% set(gca,'fontsize',fs)
% xlabel('$t$ [s]', 'Interpreter','latex', 'FontSize', 24)
% ylabel ('$\dot{E}_{mec} - P_{nc}$ [W]', 'Interpreter','latex', 'FontSize', 24)
% axis([0 18 -0.66e-6 0.66e-6])
% 
% 
% figure
% plot(t(2:end-1),Pnc(2:end-1), 'linewidth', lw)
% hold on 
% plot(t(2:end-1),dEmec, 'r', 'linewidth', lw)
% hold off
% set(gca,'fontsize',fs)
% xlabel('$t$ [s]', 'Interpreter','latex', 'FontSize', 24)
% ylabel ('$P_{nc}$, $dE_{mec}$ [W]', 'Interpreter','latex', 'FontSize', 24)
% axis([0 18 -0.7 0.7])
% legend('$P_{nc}$', '$E_{mec}$', 'Interpreter', 'latex')
% 
% figure
% subplot(1,2,1)
% plot(t,theta, 'r', 'linewidth', lw)
% set(gca,'fontsize',fs)
% xlabel('$t$ [s]', 'Interpreter','latex', 'FontSize', 24)
% ylabel ('$\theta$ [$\mathrm{rad}$]', 'Interpreter','latex', 'FontSize', 24)
% axis([0 18 -2 2])
% subplot(1,2,2)
% plot(t,thetadot, 'r', 'linewidth', lw)
% set(gca,'fontsize',fs)
% xlabel('$t$ [s]', 'Interpreter','latex', 'FontSize', 24)
% ylabel ('$\dot{\theta}$ [$\mathrm{rad}\cdot\mathrm{s}^{-1}$]', 'interpreter','latex', 'FontSize', 24)
% axis([0 18 -12 12])
% 
% figure
% subplot(1,2,1)
% plot(t(2:end-1),Pnc(2:end-1), 'linewidth', lw)
% set(gca,'fontsize',fs)
% xlabel('$t$ [s]', 'Interpreter','latex', 'FontSize', 24)
% ylabel ('$P_{nc}$ [W]', 'Interpreter','latex', 'FontSize', 24)
% axis([0 18 -0.7 0.7])
% subplot(1,2,2)
% plot(t(2:end-1),dEmec, 'r', 'linewidth', lw)
% xlabel('$t$ [s]', 'Interpreter','latex', 'FontSize', 24)
% ylabel ('$\dot{E}_{mec}$ [W]', 'Interpreter','latex', 'FontSize', 24)
% axis([0 18 -0.7 0.7])

figure
scatter(mod(theta-pi,2*pi)-pi,thetadot,'blue', 'filled');
xlabel('$\theta$ [rad]', 'Interpreter', 'latex', 'FontSize', 24)
ylabel('$\dot{\theta}$ [rad$\cdot$s$^{-1}$]','Interpreter', 'latex', 'FontSize', 24)
axis([-3.2 3.2 -25 25])