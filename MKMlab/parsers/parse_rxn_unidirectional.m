function [species,nu] = parse_rxn_unidirectional(eq_uni)
%gets species and stoichiometric numbers for one side of a reaction equation

%pre-allocate
species_with_num = {};
i = 1;

%define location of + signs
loc_plus = strfind(eq_uni,'+');

%if there are no + signs in eq, do not modify string
if isempty(loc_plus) == true
    species_with_num{i} = eq_uni;
end

%if one + sign, split into two sets of species/coeffs
if length(loc_plus) == 1
    species_with_num{i} = eq_uni(1:loc_plus-1);
    species_with_num{i+1} = eq_uni(loc_plus+1:end);
end

%if two or more + signs, split into appropriate sets of species/coeffs
if length(loc_plus) >= 2
    species_with_num{i} = eq_uni(1:loc_plus(1)-1);
    for j = 2:length(loc_plus)
        i = i + 1;
        species_with_num{i} = eq_uni(loc_plus(j-1)+1:loc_plus(j)-1);
    end
    species_with_num{i+1} = eq_uni(loc_plus(end)+1:end);
end

%parse species and coeffs from string and put in cell arrays
species = cell(1,length(species_with_num));
nu = zeros(1,length(species_with_num));
for k = 1:length(species_with_num)
    
    %locate coefficient as number preceeding species
    loc_stoich_num = regexp(species_with_num{k},'(?<=^)\d(?i)([A-Z]|\*)');
    nu(k) = str2double(species_with_num{k}(1:loc_stoich_num));
    
    %if there is no coefficient written, assume it's 1
    if isempty(loc_stoich_num) == true
        nu(k) = 1;
    end
    
    species{k} = species_with_num{k}(regexp(species_with_num{k},'(?i)([A-Z]|\*)'):end);
    
end

end