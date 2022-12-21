fs = 16;
lw = 2;
close all
filename = 'output_example.out';
data = load(filename); 

%paramètres, DOIVENT ETRE IDENTIQUES QUE CEUX DANS LE CONFIGURATION
g = 9.81;
gamma = 1.4;
rho_0 = 1.2;
P_0 = 100000;
K = P_0 * rho_0^ (-gamma);

%Recupération des données
z           = data(:,1); 
rho       = data(:,2);
rho_ana = ( rho_0.^(gamma-1) - g.*z.* (gamma-1) ./(K*gamma) ).^(1/(gamma-1));
erreur = abs(rho - rho_ana);
%plots
figure(1)
plot(z,rho,'LineWidth',lw-0.5)
hold on
plot(z,rho_ana,'--','LineWidth',lw+1)
grid on 
set (gca,'fontsize',fs)
xlabel('altitude [m]','Interpreter','latex')
ylabel('$\rho$ [kg$\cdot$m$^{-3}$]','Interpreter','latex')
legend('$\rho_{num}$ [kg$\cdot$m$^{-3}$]','$\rho_{ana}$ [kg$\cdot$m$^{-3}$]','Interpreter','latex')

figure(2)
plot(z,erreur,'b.','LineWidth',lw)
xlabel ('Hauteur [m]','Interpreter','latex');
ylabel ('erreur [kg$\cdot$m$^{-3}$]','Interpreter','latex');
set (gca,'fontsize',fs)
grid on 