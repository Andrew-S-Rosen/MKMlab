filename = 'NH3_synthesis.mkm';
T = linspace(400,700,100);
gases = {'N2','H2'};
P_i = [25 75];
TOF_species = 'NH3';

TOF = zeros(1,length(T));
for i = 1:length(T)
    sol = run_mkm(filename,T(i),gases,P_i,'TOF_species',TOF_species);
    TOF(i) = sol.TOF;
end

plot(T,TOF)
xlabel('Temperature (K)')
ylabel('TOF (site^{-1} s^{-1})')