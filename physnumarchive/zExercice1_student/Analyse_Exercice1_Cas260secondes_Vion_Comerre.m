%CAS 260 SECONDES
 vfin_num = [7039.4 9751.06666666667 12559.9065600035 15419.2807635328 18304.4129715712 21202.5471447469 24107.2130569549]; % vous complétez ici 'a la main' 
 lw=2; fs=16;
 figure
 plot(1./nsteps_num, vfin_num, 'k+-','linewidth',lw)
 set(gca,'fontsize',fs)
 xlabel('1/N_{steps}')
 ylabel('v_{final}')
 grid on


 xfin_num = [328057.6 473736.466666667 580564.314719955 651907.225037037 696787.932842394 723895.498389787 739798.708141224]; % vous complétez ici 'a la main' 
 lw=2; fs=16;
 figure
 plot( (260./nsteps_num).^0.7, xfin_num, 'k+-', 'linewidth' ,lw);%ordre 0.70?
 set(gca,'fontsize',fs)
 xlabel('(\Delta t)^{0.7}[s]')
 ylabel('x_{final}[m]')
 grid on