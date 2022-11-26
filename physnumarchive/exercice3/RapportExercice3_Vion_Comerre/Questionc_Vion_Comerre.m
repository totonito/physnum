data = load('output_example.out');
t    = data(:,1);
theta   = data(:,2);
thetadot = data(:,3);
fs = 16;
figure(1)
plot (t,theta,'r-','linewidth',1)
hold on
set (gca,'fontsize',fs)
grid on
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\theta$ [rad]','Interpreter','latex')

figure(2)
plot(t,thetadot,'r-','linewidth',1)
hold on
set (gca,'fontsize',fs)
grid on
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\dot\theta$ [rad $\cdot$ s$^{-1}$]','Interpreter','latex')