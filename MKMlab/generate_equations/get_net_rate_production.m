function [eqn_handle_str,site_species] = get_net_rate_production(filename,T,gases_str,provided_P,eqn_type,log_val,TOF_species)
%generate string representation for net rates of production
%
%INPUTS:
%filename - string: filename of reaction mechanism file
%T - double: absolute temperature
%gases_str - cell array of strings: gas-phase species
%provided_P - vector of doubles: partial pressures for each gas species
%
%OPTIONAL INPUTS:
%log_val - logical: true if information should be printed to screen or false if not (default: true)
%eqn_type - double: (1) gets system of ODEs; (2) gets system of ODEs; (3) gets equation for calculating TOF
%TOF_species - string: species to obtain TOF equation for (only included if eqn_type = 3)
%
%OUTPUTS:
%eqn_handle_str - string: string representation of function handle for system of equations
%site_species - cell array of strings: list of all site species, including the bare catalyst site

%set defaults for DAE and log_val if not specified
if nargin == 5
    log_val = true;
end

if exist('TOF_species','var') == true && eqn_type ~=3
    error('TOF_species can only be specified if eqn_type is set to 3')
elseif exist('TOF_species','var') == false && eqn_type == 3
    error('TOF_species must be included if eqn_type is set to 3')
end

if eqn_type ~=1 && eqn_type ~= 2 && eqn_type ~= 3
    error('Input argument eqn_type can only take on values of 1, 2, or 3')
end

%parse input file for reaction details
[reacts,prods,nu_r,nu_p,A,Ea,flag_rev] = parse_mechanism_file(filename,log_val);

%combine reactants and products to get all species (includes duplicates)
rxn_species = [reacts prods];

%find all unique species in the provided input file
species = unique([rxn_species{:}]);

%define number of rxns
n_rxns = length(flag_rev);
n_rev = sum(flag_rev);
n_rxns_tot = n_rxns + n_rev;

%for every reaction, calculate k
k = zeros(1,n_rxns_tot);
for i = 1:n_rxns_tot
    k(i) = calculate_k(A(i),Ea(i),T);
end

%define gas species as those without a *
gas_species = species(contains(species,'*') == false);

%throw error if user provides invalids gas-phase species
if length(intersect(gases_str,gas_species)) < length(gases_str)
    error('Gas species in setup file is not in mechanism')
end

%for every gas-phase species, assign partial pressure
P = zeros(1,length(gas_species));
for i = 1:length(gas_species)
    
    gas_loc = find(strcmp(gases_str,gas_species{i}));
    if isempty(gas_loc) == true
        P(i) = 0; %set pressure to 0 as default
    else
        P(i) = provided_P(gas_loc); %assign pressure based on user input
    end
end

%define adsorbed species
ads_species = species(contains(species,'*') == true & strcmp(species,'*') == false);

%define all site species (i.e. adsorbed species and bare site)
site_species = horzcat(ads_species,'*');

%pre-allocate rate expressions
R = cell(1,n_rxns); %net rate for each rxn
r_f = cell(1,n_rxns); %forward rate for each rxn
r_r = cell(1,n_rxns); %reverse rate for each rxn
rev_counter = 0; %counter for number of reversible rxns

%for every reaction line in mechanism file, calculate net rxn rate
for i = 1:n_rxns
    
    %index for reaction number, counting reversible rxns twice
    idx = i + rev_counter;
    
    %calculate forward rxn rate
    r_f{i} = get_rxn_rate(gas_species,site_species,P,reacts{i},nu_r{i},k(idx));
    
    %if reversible, account for reverse rate
    if flag_rev(i) == true
        r_r{i} = get_rxn_rate(gas_species,site_species,P,prods{i},nu_p{i},k(idx+1));
        R{i} = strcat(r_f{i},'-',r_r{i});
        rev_counter = rev_counter + 1;
    else
        R{i} = r_f{i};
    end
end

%number of equations to generate
if eqn_type == 1
    n_eqns = length(site_species);
elseif eqn_type == 2
    n_eqns = length(ads_species);
elseif eqn_type == 3
    n_eqns = 1;
end

%get expressions for net rates of production
rate_production = cell(1,n_eqns);
for i = 1:n_eqns
    
    if eqn_type == 1
        species_to_solve = site_species{i};
    elseif eqn_type == 2
        species_to_solve = ads_species{i};
    elseif eqn_type == 3
        species_to_solve = TOF_species;
    end
    
    %for each reaction, get rate of production for species i
    iter = 1;
    for j = 1:n_rxns
        
        %if species not found in rxn j, go to next iteration
        if sum(strcmp(species_to_solve,reacts{j})) == 0 && sum(strcmp(species_to_solve,prods{j})) == 0
            continue
        end
        
        %get expression for rate of production of species i in rxn j
        net_rate_j = get_rate_production(species_to_solve,reacts{j},prods{j},nu_r{j},nu_p{j},R{j});
        
        %net rate of production of species i is sum of rate of its rate of production in each rxn
        if iter == 1
            rate_production{i} = strcat('(',net_rate_j,')');
        else
            rate_production{i} = strcat(rate_production{i},'+(',net_rate_j,')');
        end
        
        iter = iter + 1;
    end
end

eqn_handle_str = strcat('[',strjoin(rate_production,';'),']');

end