data = load('output_example.out');
t    = data(:,1);
theta   = data(:,2);
thetadot = data(:,3);

theta_ana = 1e-6 * cos(omega * t);
thetadot_ana = -omega * 1e-6  * sin(omega * t);

figure(1)
plot (t,theta,'r-','linewidth',1)
hold on
plot(t,theta_ana,'k--','linewidth',lw)
set (gca,'fontsize',fs)
grid on
legend('$\theta_{num}$','$\theta_{ana}$','interpreter','latex')
xlabel('$t$ [s]','interpreter','latex')
ylabel('$\theta$ [rad]','interpreter','latex')

figure(2)
plot(t,thetadot,'r-','linewidth',1)
hold on
plot(t,thetadot_ana,'k--','linewidth',lw)
set (gca,'fontsize',fs)
grid on
legend ('$\dot\theta_{num}$','$\dot\theta_{ana}$','interpreter','latex')
xlabel('$t$ [s]','interpreter','latex')
ylabel('$\dot\theta$ [rad$\cdot$ s$^{-1}$]','interpreter','latex')