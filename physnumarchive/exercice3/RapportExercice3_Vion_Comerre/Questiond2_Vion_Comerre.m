data = load('output_example.out');
t    = data(:,1);
theta   = data(:,2);
thetadot = data(:,3);
fs = 16;
figure(1)
plot (t(200:end),theta(200:end),'r-','linewidth',0.5)
hold on
set (gca,'fontsize',fs)
grid on
xlabel('t[s]')
ylabel('\theta')
% 
% figure(2)
% plot(t(200:end),thetadot(200:end),'r-','linewidth',0.5)
% hold on
% set (gca,'fontsize',fs)
% grid on
% xlabel('t[s]')
% ylabel('thetadot')
    
figure()
plot(theta(200:end),thetadot(200:end),'r+','linewidth',0.5)
hold on
set (gca,'fontsize',fs)
grid on
xlabel('theta')
ylabel('thetadot')