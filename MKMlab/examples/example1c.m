filename = 'CO_oxidation.mkm';
T = 600;
gases = {'CO','O2','CO2'};
P_i = [0.33 0.67 1];
sol = run_mkm(filename,T,gases,P_i,'TOF_species','CO2','tspan',[0 1E-6]);
plot_coverages(sol);
fprintf('TOF for CO2 is: %g s^(-1)\n',sol.TOF)