function rxn_lines = import_rxn_lines(filename)
%import each line containing a reaction in the input file

%open file
fid = fopen(filename);

%pre-allocate
rxn_lines = cell(1);
tline = blanks(1);
i = 1;

%read file until last line
while ischar(tline)
    
    %get current line
    tline = fgetl(fid);
    
    %if line is blank, go to next line
    if sum(tline) == -1 || isempty(tline) == true
        continue
    end
    
    %remove leading/trailing whitespace in string
    rxn_lines{i} = strtrim(tline);
    
    %keep only non-whitespace characters
    rxn_lines{i} = rxn_lines{i}(regexp(rxn_lines{i},'\S'));
    
    i = i + 1;
end

%close file
fclose(fid);
end