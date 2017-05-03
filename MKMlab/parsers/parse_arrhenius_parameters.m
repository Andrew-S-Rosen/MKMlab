function [A,Ea] = parse_arrhenius_parameters(rxn_line,conditions)
%parse the Arrhenius parameters from a given reaction string

%unpack structure
k = const.k_B;
T = conditions.T;
h = const.h;

%Split components of reaction line based on commas
delim = ',';
split_line = strsplit(rxn_line,delim);

%Parse pre-factor based on 'A=' and activation energy based on 'Ea=' identifier
loc_A = contains(split_line,'A=');
loc_Ea = contains(split_line,'Ea=');
if sum(loc_A) == 0 || sum(loc_Ea) == 0
    error('Missing pre-exponential factor or activation energy for reaction: %s',rxn_line)
elseif sum(loc_A) > 1 || sum(loc_Ea) > 1
    error('Too many pre-exponential factors or activation energies for reaction: %s',rxn_line)
end
A_line = char(split_line(loc_A));
Ea_line = char(split_line(loc_Ea));

split_line = strsplit(A_line,'=');
split_fr = strsplit(split_line{2},';');
if length(split_fr) == 1
    if isnan(str2double(split_line{2})) == false
        A = str2double(split_line{2});
    else
        if strcmpi(split_line{2},'kt/h') == true
            A = k*T/h;
        end
    end
elseif length(split_fr) == 2
    if isnan(str2double(split_fr{1}(2:end))) == false
        A(1) = str2double(split_fr{1}(2:end));
    else
        if strcmpi(split_fr{1}(2:end),'kt/h') == true
            A(1) = k*T/h;
        end
    end
    if isnan(str2double(split_fr{2}(1:end-1))) == false
        A(2) = str2double(split_fr{2}(1:end-1));
    else
        if strcmpi(split_fr{2}(1:end-1),'kt/h') == true
            A(2) = const.k_B*conditions.T/const.h;
        end
    end
end
if isempty(A) == true
    error('Invalid syntax for prefactor in reaction: %s',rxn_line)
end

split_line = strsplit(Ea_line,'=');
split_fr = strsplit(split_line{2},';');
if length(split_fr) == 1
    Ea = str2double(split_line{2});
elseif length(split_fr) == 2
    Ea(1) = str2double(split_fr{1}(2:end));
    Ea(2) = str2double(split_fr{2}(1:end-1));
end

end