% Parameter scan
% This script is used to scan the parameters of the model

repertoire = './';
executable = 'Exercice4_Comerre_Vion';
input = 'configuration.in.example';
% repertoire = ";%./Users/linleyvion/Desktop/Rapports PhysNum/Exercice 4";

%% Simulations %%
% ----------
nsteps = linspace(20000,259200,30);
nsimul = numel(nsteps);

tfin = 259200.0e0;
dt = tfin./nsteps;

paramstr = 'nsteps';
param = nsteps;

output = cell(1,nsimul);
for i = 1:nsimul
    output{i} = [paramstr,'=',num2str(param(i)),'.out']
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i})
    disp(cmd)
    system(cmd)
    disp('Done.')
end

%% Variables %%
% ---------

% ---------- Theorique ----------
% x_theo = theo(:,2);
% y_theo = theo(:,3);
% vX_theo = theo(:,4);
% vY_theo = theo(:,5);

RT = 6378137.0e0;
theta = linspace(0,2*pi,200);
hmin = zeros(1,nsimul);
vmax = zeros(1,nsimul);

for i = 1:nsimul % Parcours des resultats de toutes les simulations

    data        = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    %%Atmosphère
    %rho           = data(end,1); % the index (i.e. 1,2,4.. ) can change according to how you save data in the c++ code
    %z             = data(end,2);

    %Pas de temps fixe => param = nsteps
    t   = data(:,1);
    xA  = data(:,2);
    yA  = data(:,3);
    vxA = data(:,4);
    vyA = data(:,5);


    r=sqrt(xA.^2 + yA.^2);
    v=sqrt(vxA.^2 + vyA.^2);

    %r1=data(:,11);
    %v1=data(:,12);
    

    [alpha,beta]=min(abs(r));
    [gamma,delta]=max(abs(v));


    ALPHA=[r(beta-1) r(beta) r(beta+1)];
    GAMMA=[v(delta-1) v(delta) v(delta+1)];
    %GAMMA=[v(beta-1) v(beta) v(beta+1)];
    ta=[t(beta-1) t(beta) t(beta+1)];
    tg=[t(delta-1) t(delta) t(delta+1)];

    pa=polyfit(ta,ALPHA-RT,2);
    test=linspace(t(beta-1),t(beta+1),500);
    qa=polyval(pa,test);
    pg=polyfit(tg,GAMMA-11134.284,2);

    vmax(i)= -(pg(2)^2)/(4*pg(1)) + pg(3);
    hmin(i)=-(pa(2)^2)/(4*pa(1)) + pa(3);
    nsteps2(i)=length(t)-1;

    figure(i)
    plot(t,abs(r-RT),'b*','linewidth',3)
    hold on
    plot(test,qa,'r.-')
    xlim([t(beta-1) t(beta+1)])
    grid on
    hold off  
end

% for i = 1:nsimul
%     data = load(output{i});
%     t = data(:,1);
%     x = data(:,2);
%     y = data(:,3);
%     vX = data(:,4);
%     vY = data(:,5);
%     h = sqrt(x.^2+y.^2);
%     v = sqrt(vX.^2+vY.^2);
%     
%     [val_h,ind_h] = min(abs(h));
%     [val_v,ind_v] = max(abs(v));
%     
%     ALPHA = [h(ind_h-1) h(ind_h) h(ind_h+1)];
%     BETA = [v(ind_v-1) v(ind_v) v(ind_v+1)];
%     
%     ta = t(ind_h-1:ind_h+1);
%     tb = t(ind_v-1:ind_v+1);
%     
%     a = polyfit(ta,ALPHA-RT,2);
%     test_a = linspace(t(ind_h-1),t(ind_h+1),500);
%     q_a = polyval(a,test_a);
% 
%     b = polyfit(tb,BETA-11134.28343,2);  
% 
%     h_min(i) = -(a(2)^2)/(4*a(1)) + a(3);
%     v_max(i) = -b(2)^2/(4*b(1)) + b(3);
%     nsteps2(i) = length(t)-1;
% 
%     figure(i)
%     plot(t,abs(h-RT),'b*','linewidth',3)
%     hold on
%     plot(test_a,q_a,'r.-')
%     xlim([t(ind_h-1) t(ind_h+1)])
%     grid on
%     hold off
%     end

lw = 1; fs = 24;

% figure
% loglog(dt,error,'-o','LineWidth',lw)
% set(gca,'FontSize',fs)
% xlabel('$\Delta t$ [s]', 'Interpreter', 'latex')
% ylabel('error [m]', 'Interpreter', 'latex')
% grid on

%% PLOT %%

figure
loglog(dt, hmin, '+','LineWidth',lw)
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize',fs);
xlabel('$\Delta t$ [s]', 'Interpreter', 'latex')
ylabel('minimum altitude [m]', 'Interpreter', 'latex')
grid on

figure
loglog(dt, vmax,'+','LineWidth',lw)
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize',fs);
xlabel('$\Delta t$ [s]', 'Interpreter', 'latex')
ylabel('$v_{max}$ [m/s]', 'Interpreter', 'latex')
grid on

    

