function sol = root_solver(guess,rxn_struct,options)

eqn_handle_str = rxn_struct.eqn_handle_str;
n_vars = length(rxn_struct.site_species);
log = options.log;
DAE = options.ODE_options.DAE;
optim_options = options.optim_options;

%solves system of equations for roots using fsolve

if log == true
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
if log == true
    fprintf('d(theta)/dt:\n')
    disp(fun(f_sol))
    fprintf('\nError in site balance: %g\n',1-sum(f_sol))
end

if abs(1-sum(f_sol)) >= 0.01
    error('Site balance has too much error after root-finding procedure.')
end
end