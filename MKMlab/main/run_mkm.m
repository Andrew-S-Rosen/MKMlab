function sol = run_mkm(filename,T,gases_str,provided_P,varargin)
%given input file and experimental conditions, construct and solve a microkinetic model
%
%INPUTS:
%filename - string: filename of reaction text file
%T - double: absolute temperature
%gases_str - cell array of strings: gas-phase species
%provided_P - vector of doubles: partial pressures for each gas species
%
%OPTIONAL INPUTS (positional values):
%sites_str - cell array of strings: site species, including bare sites, used to define initial conditions (default: {'*'})
%provided_y0 - vector of doubles: initial fractional coverages for each site species (default: 1 for bare site)
%
%OPTIONAL INPUTS (name-value pairs):
%ODE_algorithm - string: ODE solver to use. options include ode15s, ode23s, ode23t, ode23tb, ode45, ode23, and ode 113 (default: 'ode15s')
%tspan - column vector of length 2: time span for integration. if not specified, QSS is assumed (default: [])
%DAE - logical: true if conservation law is explicitly included in system of equations (i.e. system of DAEs) or false if not (i.e. system of ODEs) (default: false)
%TOF_species - string: species to evaluate TOF for
%ss_deriv_tol - double: steady-state criterion based on d(theta)/dt for ODE/DAE solver (default: 1E-6)
%ss_diff_tol - double: steady-state criterion based on relative change in theta for ODE/DAE solver (default: 1E-10)
%ODE_rel_tol - double: relative tolerance specified in odset option (default: 1E-3)
%ODE_abs_tol - double: absolute tolerance specified in odset option (default: 1E-6)
%root_algorithm - string: root-finding algorithm for fsolve specified in optimoptions (default: 'levenberg-marquardt')
%root_max_fun_eval - double: maximum function evaluations for fsolve specified in optimoptions (default: 10000)
%root_fun_tol - double: function tolerance for fsolve specified in optimoptions (default: 1E-6)
%root_max_iter - double: maximum iterations for fsolve specified in optimoptions (default: 10000)
%root_step_tol - double: step tolerance for fsolve specified in optimoptions (default: 1E-6)
%log - logical: true if information should be printed to screen or false if not (default: false)

%OUTPUTS:
%sol - structure: solution to the microkinetic model, containing the following fields:
%sol.theta - N x M matrix of doubles: fractional coverages for M species over N time points (N = 1 if QSS is implied)
%sol.t - vector of doubles of length N: time points for integration
%sol. site_species - cell array of strings of length N: adsorbed species, including the bare catalyst site, corresponding to sol.theta
%sol.TOF - double: TOF for specified TOF_species

%ensure user provides partial pressure for each gas species listed
if length(gases_str) ~= length(provided_P)
    error('Array of partial pressures is not same length as array of gas species')
end

%ensure temperature is positive
if T <= 0
    error('Temperature is in Kelvins (must be > 0)')
end

%ensure pressures are positive
if sum(provided_P < 0) > 0
    error('Partial pressures must be >= 0')
end

%parse input
input_results = parse_inputs(filename,T,gases_str,provided_P,varargin);

%define optional parameters based on user input
sites_str = input_results.sites_str;
provided_y0 = input_results.provided_y0;
tspan = input_results.tspan;
DAE = input_results.DAE;
ODE_algorithm = input_results.ODE_algorithm;
ss_deriv_tol = input_results.ss_deriv_tol;
ss_diff_tol = input_results.ss_diff_tol;
ODE_rel_tol = input_results.ODE_rel_tol;
ODE_abs_tol = input_results.ODE_abs_tol;
log_val = input_results.log;
root_algorithm = input_results.root_algorithm;
root_max_fun_eval = input_results.root_max_fun_eval;
root_fun_tol = input_results.root_fun_tol;
root_max_iter = input_results.root_max_iter;
root_step_tol = input_results.root_step_tol;
TOF_species = input_results.TOF_species;

if sum(provided_y0) ~= 1
    error('Sum of initial coverages does not equal 1')
end

if DAE == false
    coverage_eqn_type = 1;
else
    coverage_eqn_type = 2;
end

%set ODE and fsolve options
ODE_options = struct('ODE_algorithm',ODE_algorithm,'DAE',DAE,'ss_deriv_tol',ss_deriv_tol,'ss_diff_tol',ss_diff_tol,'ODE_rel_tol',ODE_rel_tol,'ODE_abs_tol',ODE_abs_tol);
optim_options = optimoptions(@fsolve,'Display','off','FiniteDifferenceType','central','Algorithm',root_algorithm,'MaxFunctionEvaluations',root_max_fun_eval,'FunctionTolerance',root_fun_tol,'MaxIterations',root_max_iter,'StepTolerance',root_step_tol);

%print experimental conditions
if log_val == true
    fprintf('*******************************************\nExperimental Conditions:\n\n')
    fprintf('Temperature: %g K\n',T)
    disp(table(gases_str',provided_P','VariableNames',{'Gases','Partial_Pressures'}))
    fprintf('\n')
    disp(table(sites_str',provided_y0,'VariableNames',{'Species','Initial_Coverage'}))
end

%get rate of production for each site species
[eqn_handle_str,site_species] = get_net_rate_production(filename,T,gases_str,provided_P,coverage_eqn_type,log_val);

%Solve the system of equations for coverages
if length(sites_str) ~= length(provided_y0)
    error('Array of initial fractional coverages is not same length as array of site species')
end

%ensure provided site species exist in rxn mechanism
if length(intersect(sites_str,site_species)) < length(sites_str)
    error('Site species in setup file is not in mechanism')
end

%definte the number of variables in system of equations
n_vars = length(site_species);

%for every site species, assign initial coverage based on user-input
y0 = zeros(1,n_vars);
for i = 1:n_vars
    
    %find location of user-provided site species in the sites_str list
    site_loc = find(strcmp(sites_str,site_species{i}));
    
    if isempty(site_loc) == true
        y0(i) = 0; %default is bare (except for '*', which is 1)
    else
        y0(i) = provided_y0(site_loc); %assign based on user-input
    end
end

%solve for coverages using ODE/DAE solver
ODE_sol = ODE_solver(eqn_handle_str,y0,tspan,log_val,ODE_options);

%solve for roots using fsolve if QSS is implied using ODE_sol as initial guess
if isempty(tspan) == true && license('test','Optimization_Toolbox') == true
    guess = ODE_sol.theta';
    root_sol = root_solver(eqn_handle_str,n_vars,guess,log_val,DAE,optim_options);
    sol = root_sol;
elseif isempty(tspan) == true && license('test','Optimization_Toolbox') == false
    warning('Optimization toolbox not installed. Steady-state coverages determined from integration only')
    sol = ODE_sol;
else
    sol = ODE_sol;
end

%print coverage results
if log_val == true
    if isempty(tspan) == true
        fprintf('*******************************************\nResults at steady-state:\n\n')
        disp(table(site_species',sol.theta','VariableNames',{'Species','Coverage'}))
    else
        fprintf('*******************************************\nResults at end-time:\n\n')
        disp(table(site_species',sol.theta(:,end),'VariableNames',{'Species','Coverage'}))
    end
end

%include site species in output structure
[sol.site_species] = deal(site_species);

%solve for TOF and add to output structure
if isempty(TOF_species) == false
    TOF_eqn_type = 3;
    TOF_eqn_handle_str = get_net_rate_production(filename,T,gases_str,provided_P,TOF_eqn_type,log_val,TOF_species);
    TOF = calculate_TOF(TOF_eqn_handle_str,sol.theta');
    [sol.TOF] = deal(TOF);
end

end