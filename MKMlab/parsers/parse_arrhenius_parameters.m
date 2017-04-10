function [A,Ea] = parse_arrhenius_parameters(rxn_line)
%parse the Arrhenius parameters from a given reaction string
%
%INPUTS:
%rxn_line - string: a reaction line from the input file
%
%OUTPUTS:
%A - double or vector: Arrhenius pre-factor(s)
%Ea - double or vector: activation energy/energies

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
A = str2num(split_line{2});
split_line = strsplit(Ea_line,'=');
Ea = str2num(split_line{2});
end