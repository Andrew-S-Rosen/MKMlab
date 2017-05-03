function [reacts,prods,nu_r,nu_p,flag_rev] = parse_rxn(eq)
%parses reaction string for species and stochiometric numbers

%Convert equation string into consistent format
eq = strrep(eq,'->','>');
eq = strrep(eq,'=>','>');
if contains(eq,'<-') == true || contains(eq,'<=') == true
    error('Only forward (->) and reversible (=) reactions are allowed')
end
eq = strrep(eq,'<->','=');
eq = strrep(eq,'<=>','=');

%Find location of rxn arrow and flag if rxn is reversible
if contains(eq,'>') == true
    flag_rev = false;
    loc_arrow = strfind(eq,'>');
elseif contains(eq,'=') == true
    flag_rev = true;
    loc_arrow = strfind(eq,'=');
else
    error('Missing -> or = operator')
end

%Define left (reactants) and right (products) sides of rxn
eq_left = eq(1:loc_arrow-1);
eq_right = eq(loc_arrow+1:end);

%Get reactants, products, and their stoichiometric numbers
[reacts,nu_r] = parse_rxn_unidirectional(eq_left);
if length(unique(reacts)) ~= length(reacts)
    error('Identical species not grouped together in reaction: %s',eq_left)
end
[prods,nu_p] = parse_rxn_unidirectional(eq_right);
if length(unique(prods)) ~= length(prods)
    error('Identical species not grouped together in reaction: %s',eq_right)
end
%Confirm number of * on reactants side is same as on products side
if sum(nu_r(contains(reacts,'*'))) ~= sum(nu_p(contains(prods,'*')))
    warning('Unabalanced number of sites in reaction equation: %s',eq)
end
end