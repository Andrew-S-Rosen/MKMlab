function rxn_struct = get_net_rate_production(filename,conditions,options)
%generate string representation for net rates of production

%unpack conditions and options structures
T = conditions.T;
gases_str = conditions.gases;
provided_P = conditions.provided_P;
DAE = options.ODE_options.DAE;
TOF_species = options.TOF_species;

%parse input file for reaction details
rxn_struct = parse_mechanism_file(filename,conditions,options);

%unpack rxn structure
reacts = rxn_struct.reacts;
prods = rxn_struct.prods;
nu_r = rxn_struct.nu_r;
nu_p = rxn_struct.nu_p;
A = rxn_struct.A_vec;
Ea = rxn_struct.Ea_vec;
flag_rev = rxn_struct.flag_rev;

%combine reactants and products to get all species (includes duplicates)
rxn_species = [reacts prods];

%find all unique species in the provided input file
species = unique([rxn_species{:}]);

%define number of rxns
n_rxns = length(flag_rev);

%for every reaction, calculate k
k = calculate_k(A,Ea,T);

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
if DAE == false
    n_eqns = length(site_species);
else
    n_eqns = length(ads_species);
end

%get expressions for net rates of production
rate_production = cell(1,n_eqns);

for i = 1:n_eqns
    
    if DAE == false
        species_to_solve = site_species{i};
    else
        species_to_solve = ads_species{i};
    end
    
    %for each reaction, get rate of production for species i
    iter = 1;
    for j = 1:n_rxns
        
        %if species not found in rxn j, go to next iteration
        if sum(strcmp(species_to_solve,reacts{j})) == 0 && sum(strcmp(species_to_solve,prods{j})) == 0
            continue
        end
        
        %get expression for rate of production of species i in rxn j
        net_rate_j = get_rate_production(species_to_solve,reacts{j},...
            prods{j},nu_r{j},nu_p{j},R{j});
        
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

if isempty(TOF_species) == false

    species_to_solve = TOF_species;
    %for each reaction, get rate of production for species i
    iter = 1;
    for j = 1:n_rxns

        %if species not found in rxn j, go to next iteration
        if sum(strcmp(species_to_solve,reacts{j})) == 0 && sum(strcmp(species_to_solve,prods{j})) == 0
            continue
        end

        %get expression for rate of production of species i in rxn j
        net_rate_j = get_rate_production(species_to_solve,reacts{j},...
            prods{j},nu_r{j},nu_p{j},R{j});

        %net rate of production of species i is sum of rate of its rate of production in each rxn
        if iter == 1
            rate_production = {strcat('(',net_rate_j,')')};
        else
            rate_production = {strcat(rate_production,'+(',net_rate_j,')')};
        end

        iter = iter + 1;
    end

    TOF_eqn_handle_str = strcat('[',strjoin(rate_production,';'),']');
    [rxn_struct.TOF_eqn_handle_str] = TOF_eqn_handle_str;

end

%add fields to rxn_struct
[rxn_struct.k] = deal(k);
[rxn_struct.P] = deal(P);
[rxn_struct.site_species] = deal(site_species);
[rxn_struct.gas_species] = deal(gas_species);
[rxn_struct.n_rxns] = deal(n_rxns);
[rxn_struct.eqn_handle_str] = deal(eqn_handle_str);

end