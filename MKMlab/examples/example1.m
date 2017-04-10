filename = 'CO_oxidation.mkm';
T = 600;
gases = {'CO','O2','CO2'};
P_i = [0.33 0.67 1];
TOF_species = 'CO2';
sol = run_mkm(filename,T,gases,P_i,'TOF_species',TOF_species);
plot_coverages(sol);
fprintf('TOF is %E s^(-1)\n',sol.TOF)