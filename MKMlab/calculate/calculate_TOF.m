function TOF = calculate_TOF(eqn_handle_str,theta)
%Calculate turnover frequency given equation and values of theta
%
%INPUTS:
%eqn_handle_str - string: a string representing the function handle for the system of ODEs
%theta - vector of doubles: coverages for each spcies in system of ODEs
%
%OUTPUTS:
%TOF - double: turnover frequency for desired species

%convert string to function handle
fun = str2func(strcat('@(theta)',eqn_handle_str));

%substitute in coverages and solve for TOF of desired species
TOF = fun(theta);

end