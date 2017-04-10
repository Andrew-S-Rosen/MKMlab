function sol = root_solver(eqn_handle_str,n_vars,guess,log_val,DAE,optim_options)
%solves system of equations for roots using fsolve
%INPUTS:
%eqn_handle_str - string: string representation of function handle for system of equations
%n_vars - double: number of variables (i.e. number of unique site species)
%guess - vector of doubles: initial guess for each variable in system of equations
%
%OPTIONAL INPUTS:
%log_val - logical: true if information should be printed to screen or false if not (default: true)
%DAE - logical: true if system of DAEs is to be solved or false if ODEs are to be solved (default: false)
%optim_options - structure: structure containing fsolve options (see optimoptions optional arguments in run_mkm for details)
%
%OUTPUTS:
%sol - structure: solution to the microkinetic model, containing the following field:
%sol.theta - array of doubles: fractional coverages for each site species

if log_val == true
    fprintf('*******************************************\nAlgebraic Root-Finding:\n\n')
end

%generate function handle from string (if DAE is enabled, include site balance)
if DAE == true
    theta_vars = cell(1,n_vars);
    for i = 1:n_vars
        theta_vars{i} = strcat('theta(',num2str(i),')');
    end
    site_balance = strjoin(theta_vars,'+');
    site_balance = strcat(site_balance,'-1');
    fun = str2func(strcat('@(theta)',eqn_handle_str(1:end-1),';',site_balance,']'));
else
    fun = str2func(strcat('@(theta)',eqn_handle_str));
end

%use fsolve to find roots
f_sol = fsolve(fun,guess,optim_options);

%compute difference in integration and root-finding procedures
err_algebraic = f_sol-guess;
if max(abs(err_algebraic)) >= 0.01
    warning('Solution found via root-finding may differ from the steady-state ODE solution. Try setting the steady-state ODE tolerance to be lower')
end

%store fsolve solution to output structure
sol = struct('theta',f_sol);

%print root-finding details
if log_val == true
    fprintf('d(theta)/dt:\n')
    disp(fun(f_sol))
    fprintf('\nError in site balance: %g\n',1-sum(f_sol))
end

if abs(1-sum(f_sol)) >= 0.01
    error('Site balance has too much error after root-finding procedure.')
end
end