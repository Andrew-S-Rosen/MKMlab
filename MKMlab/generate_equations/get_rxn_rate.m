function r = get_rxn_rate(gas_species,site_species,P,reactants,nu,k)
%gets string describing unidirectional reaction rate

%define string variables theta(1) theta(2) ... theta(N)
theta_vars = cell(1,length(site_species));
for i = 1:length(site_species)
    theta_vars{i} = strcat('theta(',num2str(i),')');
end

%define site species and gas species in given rxns
gases = intersect(reactants,gas_species);
sites = intersect(reactants,site_species);

%define stoichiometric numbers for site species and gas species in rxn
nu_gases = nu(contains(reactants,'*') == false);
nu_sites = nu(contains(reactants,'*') == true);

%numerical factor in rate expression
num_factor = k;

%multiply effect of each partial pressure in rxn
idx_P = zeros(1,length(gases));
for i = 1:length(gases)
    idx_P(i) = find(strcmp(gas_species,gases{i}) == true);
    num_factor = num_factor*P(idx_P(i))^(nu_gases(i));
end

%get effect of every coverage variable in rxn
coverage_terms = cell(1,length(sites));
for i = 1:length(sites)
    theta_i = theta_vars{strcmp(site_species,sites{i}) == true};
    coverage_terms{i} = strcat(theta_i,'^',num2str(nu_sites(i)));
end

%rxn equation is numerical factor times coverage terms
r = strcat('(',strjoin([{num2str(num_factor,16)},coverage_terms],'*'),')');

end