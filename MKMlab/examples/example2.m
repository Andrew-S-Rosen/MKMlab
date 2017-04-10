filename = 'NH3_synthesis.mkm';
T = 700;
gases = {'N2','H2','NH3'};
P_tot = 100;
P_i = [0.25*P_tot 0.75*P_tot,10];
TOF_species = 'NH3';
sol = run_mkm(filename,T,gases,P_i,'TOF_species',TOF_species);
plot_coverages(sol);
