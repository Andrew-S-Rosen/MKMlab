filename = 'CO_oxidation.mkm';
T = 600;
gases = {'CO','O2','CO2'};
P_i = [0.33 0.67 1];
tspan = [0 1E-11];
sol = run_mkm(filename,T,gases,P_i,'tspan',tspan,'log',true);
plot_coverages(sol);