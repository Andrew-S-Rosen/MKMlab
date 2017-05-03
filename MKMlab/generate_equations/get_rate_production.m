function net_rate_production = get_rate_production(species_to_solve,reactants,products,nu_r,nu_p,R)
%gets equation for net rate of production of a species for a given reaction

%find index of species in reactant or product list for given rxn
loc_react = strcmp(species_to_solve,reactants);
loc_prod = strcmp(species_to_solve,products);

%if species is a reactant in rxn
if sum(loc_react) == true
    net_rate_production = strcat(num2str(-nu_r(loc_react)),'*(',R,')');
end

%if species is a product in rxn
if sum(loc_prod) == true
    net_rate_production = strcat(num2str(nu_p(loc_prod)),'*(',R,')');
end

end