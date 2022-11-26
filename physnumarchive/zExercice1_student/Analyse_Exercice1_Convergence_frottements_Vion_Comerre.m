 nsteps_num  = [5 10 20 40 80 160 320]; % vous complétez ici 'a la main' 
 vfin_num = [1190.6 1284.1 1331.6 1355.5 1367.5 1373.45944322111 1376.45262447741]; % vous complétez ici 'a la main' 
 lw=2; fs=16;
 figure
 plot(1./nsteps_num, vfin_num, 'k+-','linewidth',lw)
 set(gca,'fontsize',fs)
 xlabel('1/N_{steps}')
 ylabel('v_{final}[m.s^{-2}]')
 grid on

 xfin_num = [45088.1 55130.7 60768.9 63745.7 65273.8 66047.7684902886 66437.2389532103]; % vous complétez ici 'a la main' 
 lw=2; fs=16;
 figure
 plot(1./nsteps_num, xfin_num, 'k+-','linewidth',lw)
 set(gca,'fontsize',fs)
 xlabel('1/N_{steps}')
 ylabel('x_{final}')
 grid on


