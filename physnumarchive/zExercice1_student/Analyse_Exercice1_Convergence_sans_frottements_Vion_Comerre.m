 nsteps_num  = [5 10 20 40 80 160 320]; % vous complétez ici 'a la main' 
 vfin_num = [1263.89 1341.51 1381.68 1402.11 1412.41 1417.58 1420.17030453103]; % vous complétez ici 'a la main' 
 lw=2; fs=16;
 figure
 
 plot(1./nsteps_num, vfin_num, 'k+-','linewidth',lw)
 yline(1423.96,'r','linewidth',2)
 set(gca,'fontsize',fs)
 xlabel('1/N_{steps}')
 ylabel('v_{final}[m]')
 grid on
 

% si on a la solution analytique:
 vfin_ana = 1423.96;
 error_vfin = vfin_num-vfin_ana;
 figure
 plot(nsteps_num, abs(error_vfin),'k+-','linewidth',lw)
 set(gca,'fontsize',fs)
 set(gca,'xscale','log')
 set(gca,'yscale','log')
 
 xlabel('N_{steps}')
 ylabel('Erreur sur v_{fin}[m]')
 grid on



 xfin_num = [47156.1 57587 63302.7 66292.5 67821.3 68594.2 68982.8685014566]; % vous complétez ici 'a la main' 
 lw=2; fs=16;
 figure
 plot(1./nsteps_num, xfin_num, 'k+-','linewidth',lw)
 yline(69444.9,'r','linewidth',2)
 set(gca,'fontsize',fs)
 xlabel('1/N_{steps}')
 ylabel('x_{final}[m]')
 grid on

% si on a la solution analytique:
 xfin_ana = 69444.9;
 error_xfin = xfin_num-xfin_ana;
 figure
 plot(nsteps_num, abs(error_xfin),'k+-','linewidth',lw)
 set(gca,'fontsize',fs)
 set(gca,'xscale','log')
 set(gca,'yscale','log')
 xlabel('N_{steps}')
 ylabel('Erreur sur x_{fin}[m]')
 grid on



