repertoire = ''; % Chemin d'accès au code compilé (NB: enlever le ./ sous Windows)
executable = 'Exercice4_Comerre_Vion_rho.exe';
input      = 'configuration.in.example';
fs = 16;
lw = 2;

%paramètres, DOIVENT ETRE IDENTIQUES QUE CEUX DANS LE CONFIGURATION
g= 9.81;
gamma = 1.4;
rho_0 = 1.2000000000000;
P_0 = 100000;
K = P_0 * rho_0^ (-gamma);

        nsteps = [100 200 400 800 1600 3200];
        paramstr_1 = 'nsteps'; 
        param_1 = nsteps; 
        

nsimul_1 = numel(param_1);
%% Simulations %% 

output = cell(1, nsimul_1);
for i = 1:nsimul_1
    output{i} = [paramstr_1, '=', num2str(param_1(i)), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr_1, param_1(i), output{i});
    disp(cmd)
    system(cmd);
    disp('Done.')
end
        
%% Analyse %%
erreur = zeros(1,nsimul_1);
dz    = zeros(1,nsimul_1);
rho   = zeros(1,nsimul_1);
z_init= zeros(1,nsimul_1);

for i = 1:nsimul_1 
    data        = load(output{i});  
    rho(i)       = data(end,2);
    z_init(i) = data(1,1);
    dz(i) = z_init(i)./nsteps(i);     
  erreur(i) = rho_0 - rho(i);



end

figure()
p = polyfit(log(nsteps),log(erreur),1);
z = polyval(p,log(nsteps));
z_ = exp(z);
loglog(nsteps,z_,'g--')
hold on
loglog (nsteps,erreur,'gx','LineWidth',lw)
grid on
set (gca,'fontsize',fs)
xlabel('$N_{steps}$','Interpreter','latex')
ylabel('erreur [kg$\cdot$m$^{-3}$]','Interpreter','latex')
legend('y = -3.92x + 7.07','$\epsilon = 1000$','Interpreter','latex')
