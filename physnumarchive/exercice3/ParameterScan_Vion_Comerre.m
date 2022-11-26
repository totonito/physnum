
%% Parametres %%
close all
repertoire = ''; % Chemin d'accès au code compilé (NB: enlever le ./ sous Windows)
executable = 'Exercice3.exe';
input      = 'configuration.in.example';
g= 9.81;
L= 0.2;
omega = sqrt(g/L);
fs = 16;
lw = 2;
scan_case = "scan_variation_theta"; % A MODIFIER; 

switch scan_case

    case 'scan_dt_sans_exc'
        nsteps = [200 300 400 500 600 800 1000 2000 4000 8000 16000];
        tfin = 10; %Verifier que la valeur de tfin est EXACTEMENT la meme que dans le fichier input
        dt = tfin ./ nsteps;
        paramstr = 'dt'; 
        param = dt;  
    
    case 'scan_dt_avec_exc'
       
        tfin = 10;
        nsteps = [1000 1100 1250 1400 1500 1750 2000 3000 5000 7500 10000];
        dt = tfin ./ nsteps;
        dt_carre = dt.^2;
        paramstr = 'dt'; 
        param = dt;
    
    case 'question_d1' %Poincaré
        thetadot0 = linspace(0,20,5);
        paramstr = 'thetadot0'; 
        param = thetadot0;  

    case 'scan_variation_theta' %thetadot= 2 pour non chaotique et 15 pour chaotique thetadot = 10 pour la question E
        theta0 = [0 0.00000001];
        paramstr = 'theta0'; 
        param = theta0;
    
end


nsimul = numel(param);

%% Simulations %% 

output = cell(1, nsimul);
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i)), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
    disp(cmd)
    system(cmd);
    disp('Done.')
end
        
%% Analyse %%
nsteps = 44858;                   %SI ERREUR : définir ici le nombre de pas de temps (taille de la matrice t) pour les questions d et e 
error            = zeros(1,nsimul);
theta_end_vector = zeros(1,nsimul);
max_Emec         = zeros(1,nsimul);
theta_           = zeros(nsteps,nsimul);
theta_poincare   = zeros(nsteps,nsimul);
thetadot_poincare= zeros(nsteps,nsimul);
error_           = zeros(nsteps,1);
for i = 1:nsimul 
    data        = load(output{i}); 
    t           = data(:,1); 
    theta       = data(:,2);
    thetadot    = data(:,3);
    theta_poincare(:,i) = mod(data(:,2)+pi,2*pi)-pi;%modulo 2pi
    theta_(:,i) = theta;                            %même chose sans le modulo
    thetadot_poincare(:,i) = data(:,3);


%     max_Emec(i) = max(data(:,4));
     theta_ana     = 1e-6 * cos(omega * t);
     theta_dot_ana = -omega * 1e-6  * sin(omega * t);
     %error(i)      =  sqrt((theta).^2+(thetadot).^2 ); % à modifier si la solution analytique est connue
     %theta_end_vector(i) = theta;
    error_ = sqrt( (thetadot_poincare(:,end)-thetadot_poincare(:,1)).^2 + omega^2*(theta_(:,end)-theta(:,1)).^2 );

end


switch scan_case

    case 'scan_dt_sans_exc'
   
loglog (nsteps,error,'k+','linewidth',lw)
hold on
p1 = polyfit(log(nsteps),log(error),1);
z1 = polyval(p1,log(nsteps));
z1_ = exp(z1);                       %Fit linéaire dans l'échelle loglog de l'erreur en fonction de Nsteps
loglog(nsteps,z1_, 'b--','linewidth',1)
set (gca, 'fontsize',fs)
xlabel ('N_{steps}')
ylabel ('erreur')
legend('data','y  = A MODIFIERx + C') %le coefficient est le premier élément de p, C est le deuxième élément de p
grid on

    case 'scan_dt_avec_exc'
    figure(1)
    n = 2;
    plot(dt.^n,theta_end_vector,'k+', 'linewidth', lw)
    hold on
    p2 = polyfit(dt.^n,theta_end_vector,1);
    z2 = polyval(p2,dt.^n);
    plot(dt.^n,z2,'b--','linewidth',1 )
    hold off
    set(gca,'fontsize',fs)
    grid on
    xlabel('$(\Delta t)^2$ [s$^{2}$]', 'interpreter', 'latex', 'FontSize', 20)
    ylabel ('$\theta_{final}$ $[\mathrm{rad}]$', 'interpreter', 'latex', 'FontSize', 20)
    legend('Data', '$y=-13964(\Delta t)^2+5,9167$', 'Interpreter', 'latex')


    case 'question_d1'
figure()
hold on 
map_ = [];
for i = 1:nsimul 
plot(theta_poincare(500:end,i),thetadot_poincare(500:end,i),'.','linewidth',0.2,'Color',[1-(i-1)/nsimul 0 (i-1)/nsimul]);
set (gca,'fontsize',fs)
map = [1-(i-1)/nsimul 0 (i-1)/nsimul];
map_ = [map_ ;map];
colormap(map_)%
axis([-pi pi -20 21])
grid on
xlabel('$\theta$ [rad]','Interpreter','latex')
ylabel('$\dot\theta$ [rad $\cdot$ s$^{-1}$]','Interpreter','latex')
end
    
    case 'scan_variation_theta'
figure(1)
hold on
plot (t,theta_(:,1),'r-','linewidth',1)
plot (t,theta_(:,2),'k--','linewidth',1)
set (gca,'fontsize',fs)
grid on
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\theta$ [rad]','Interpreter','latex')
legend('$\theta_0$  = 0 rad','$\theta_0$  = 10$^{-8}$ rad','Interpreter','latex')
%pbaspect([2,1,1])
hold off

figure(2)
%p3 = polyfit(t(2 : 23000,1),log(error_(2 : 23000,1)),1); % 1 a été exclu pour éviter un log(0)
%z3 = polyval(p3,t(2 : 23000,1));
%z3_= exp(z3);
semilogy(t,error_,'k-','linewidth',1)
hold on
%semilogy(t(2 : 23000,1),z3_,'b-','linewidth',2)
set (gca,'fontsize',fs)
grid on
%pbaspect([1,1,1])
%legend('erreur','y = C_1e^{0.8421t}+C_2') %A MODIFIER avec la première valeur de p3 et C_2 est la deuxième valeur de p_3 
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$d(t)$ [m $\cdot$ s $^{-1}$]','Interpreter','latex')
hold off
end

