function rxn_struct = parse_mechanism_file(filename,inputs,options)
%parses input file for species, reactions, and reaction parameters

%import each reaction line from the reaction text file
rxn_lines = import_rxn_lines(filename);

%pre-allocate
n_rxn = length(rxn_lines);
reacts = cell(1,n_rxn);
prods = cell(1,n_rxn);
nu_r = cell(1,n_rxn);
nu_p = cell(1,n_rxn);
flag_rev = false(1,n_rxn);

%total number of reactions (double-counting each reversible reaction)
n_rxn_tot = n_rxn + sum(flag_rev);

%pre-allocate
A_vec = zeros(1,n_rxn_tot);
Ea_vec = zeros(1,n_rxn_tot);
rev_counter = 0;

%for every listed reaction
for i = 1:n_rxn
    
    %split reaction string based on commas (first index is rxn)
    split_line = strsplit(rxn_lines{i},',');
    
    %get reaction details
    [reacts{:,i},prods{:,i},nu_r{:,i},nu_p{:,i},flag_rev(i)] = parse_rxn(split_line{1});
    
    %print reaction details
    if options.log == true
        fprintf('\nReaction %g\nReactants:',i)
        disp(reacts{:,i})
        fprintf('Products:')
        disp(prods{:,i})
        fprintf('Stoich (reactants):')
        disp(nu_r{:,i})
        fprintf('Stoich (products):')
        disp(nu_p{:,i})
        fprintf('Reversible (1 = true; 0 = false):')
        disp(flag_rev(i))
    end
    
    %if there is the same species on the reactant and product side, throw an error
    if isempty(intersect(reacts{:,i},prods{:,i})) == false
        error('Reactant or product present on both sides of reaction equation #%d. Simplify the expression.',i)
    end
    
    %get A and Ea for each reaction
    [A,Ea] = parse_arrhenius_parameters(rxn_lines{i},inputs);
    
    %index accounting for double-counting of reversible reactions
    idx = i+rev_counter;
    
    %if irreversible, A and Ea are scalars
    if flag_rev(i) == false
        if length(A) ~= 1 || length(Ea) ~= 1
            error('For irreversible reaction #%d, only 1 A and Ea should be supplied')
        end
        A_vec(idx) = A;
        Ea_vec(idx) = Ea;
        
        %if reversible, A and Ea are vectors of length 2
    else
        if length(A) ~= 2 || length(Ea) ~= 2
            error('For reversible reaction #%d, 2 A and Ea values should be supplied in the format A=[A1;A2], Ea=[Ea1;Ea2]',i)
        end
        A_vec(idx:idx+1) = A;
        Ea_vec(idx:idx+1) = Ea;
        rev_counter = rev_counter + 1;
        
    end
    
end

rxn_struct = struct('reacts',{reacts},'prods',{prods},'nu_r',{nu_r},...
    'nu_p',{nu_p},'A_vec',A_vec,'Ea_vec',Ea_vec,'flag_rev',flag_rev);

end