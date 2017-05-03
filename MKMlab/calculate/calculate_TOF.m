function TOF = calculate_TOF(eqn_handle_str,theta)
%calculate turnover frequency given equation and values of theta

%convert string to function handle
fun = str2func(strcat('@(theta)',eqn_handle_str));

%substitute in coverages and solve for TOF of desired species
TOF = fun(theta);

end