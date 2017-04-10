function sol = ODE_solver(eqn_handle_str,y0,tspan,log_val,ODE_options)
%solve microkinetic model for coverages via a system of ODEs or DAEs
%
%INPUTS:
%eqn_handle_str - string: expression for system of equations
%y0 - row vector: initial coverage for each site species
%tspan - two-element column vevtor: time span for integration (disables QSS assumption)
%ODE_solver - stirng: ODE solver to use (ode15s, ode23s, ode23t, ode23tb, ode45, ode23, or ode113)
%
%OPTIONAL INPUTS:
%log_val - logical: true if information should be printed to screen or false if not (default: false)
%ODE_options - structure: structure containing ODE options (see odeset optional arguments in run_mkm for details)
%
%OUTPUTS:
%sol - structure: solution to the microkinetic model, containing the following fields:
%sol.theta - N x M matrix of doubles: fractional coverages for M species over N time points (N = 1 if QSS is implied)
%sol.t - vector of doubles of length N: time points for integration

%set integration options
if nargin == 3
    ODE_algorithm = default_params.ODE_algorithm;
    DAE = default_params.DAE;
    ss_deriv_tol = default_params.ss_deriv_tol;
    ss_diff_tol = default_params.ss_diff_tol;
    ODE_rel_tol = default_params.ODE_rel_tol;
    ODE_abs_tol = default_params.ODE_abs_tol;
elseif nargin == 4 || nargin == 5
    ODE_algorithm = ODE_options.ODE_algorithm;
    DAE = ODE_options.DAE;
    ss_deriv_tol = ODE_options.ss_deriv_tol;
    ss_diff_tol = ODE_options.ss_diff_tol;
    ODE_rel_tol = ODE_options.ODE_rel_tol;
    ODE_abs_tol = ODE_options.ODE_abs_tol;
else
    error('Invalid number of input arguments')
end

odeset_options = odeset('RelTol',ODE_rel_tol,'AbsTol',ODE_abs_tol);

if nargin == 4
    log_val = true;
end

if isempty(tspan) == true
    QSS = true;
    tspan_start = [0 1E-10];
    tspan = tspan_start;
else
    QSS = false;
end

%for every theta variable, define a new one that is a function of time
n_vars = length(y0);
theta_vars = cell(1,length(y0));
for i = 1:n_vars
    theta_vars{i} = strcat('theta(',num2str(i),')');
end

%if DAE is enabled, define site balance (solved for zero)
if DAE == true
    site_balance = strjoin(theta_vars,'+');
    site_balance = strcat(site_balance,'-1');
end

%Set up system of equations
if DAE == true
    system = str2func(strcat('@(t,theta)',eqn_handle_str(1:end-1),';',site_balance,']'));
    M = eye(n_vars);
    M(end) = 0;
    odeset_options.Mass = M;
else
    system = str2func(strcat('@(t,theta)',eqn_handle_str));
    odeset_options = odeset(odeset_options,'NonNegative',1:1:n_vars);
end

if log_val == true
    fprintf('\n*******************************************\nIntegration:\n')
end

%integrate system of ODEs/DAEs subject to initial conditions
if DAE == true
    if strcmp(ODE_algorithm,'ode15s') == true
        ODE_sol = ode15s(system,tspan,y0,odeset_options);
    elseif strcmp(ODE_algorithm,'ode23t') == true
        ODE_sol = ode23t(system,tspan,y0,odeset_options);
    else
        error('Only ode15s or ode23t can be used to solve system of DAEs')
    end
else
    if strcmp(ODE_algorithm,'ode15s') == true
        ODE_sol = ode15s(system,tspan,y0,odeset_options);
    elseif strcmp(ODE_algorithm,'ode23s') == true
        ODE_sol = ode23s(system,tspan,y0,odeset_options);
    elseif strcmp(ODE_algorithm,'ode23t') == true
        ODE_sol = ode23t(system,tspan,y0,odeset_options);
    elseif strcmp(ODE_algorithm,'ode23tb') == true
        ODE_sol = ode23tb(system,tspan,y0,odeset_options);
    elseif strcmp(ODE_algorithm,'ode45') == true
        ODE_sol = ode45(system,tspan,y0,odeset_options);
    elseif strcmp(ODE_algorithm,'ode23') == true
        ODE_sol = ode23(system,tspan,y0,odeset_options);
    elseif strcmp(ODE_algorithm,'ode113') == true
        ODE_sol = ode113(system,tspan,y0,odeset_options);
    else
        error('Only ode15s, ode23s, ode23t, ode23tb, ode45, ode23, or ode113 can be used')
    end
end

%calculate relative difference and derivative
theta_diff = ODE_sol.y(:,end)-ODE_sol.y(:,end-1);
deriv = system(1,ODE_sol.y(:,end)');

%print integration details
if log_val == true
    fprintf('\nTime: %.1e\nRelative change in coverages:\n',tspan(2))
    disp(theta_diff)
    fprintf('d(theta)/dt:\n')
    disp(deriv)
end

%if QSS is enabled, integrate (increasing end-time by 10x) until steady-state is achieved
if QSS == true
    max_deriv = Inf;
    theta_diff = Inf;
    j = 0;
    t_final = tspan(2)*10;
    
    while max_deriv > ss_deriv_tol || max(abs(theta_diff)) > ss_diff_tol
        
        if j >= 40
            warning('System of ODEs could not be solved to given steady-state tolerance of %g\n',ss_deriv_tol)
        end
        
        %extend ODE solution from prior integration
        ODE_sol = odextend(ODE_sol,[],t_final);
        
        %calculate relative difference and derivative
        theta_diff = ODE_sol.y(:,end)-ODE_sol.y(:,end-1);
        deriv = system(1,ODE_sol.y(:,end)');
        max_deriv = max(abs(deriv));
        
        %print integration details
        if log_val == true
            fprintf('\nTime: %.1e\nRelative change in coverages:\n',t_final)
            disp(theta_diff)
            fprintf('d(theta)/dt:\n')
            disp(deriv)
        end
        
        j = j + 1;
        t_final = t_final*10;
    end
    
    %print coverages at end of integration procedure
    if log_val == true
        fprintf('\nSteady-state reached!\nCoverages at end of integration procedure:\n')
        disp(ODE_sol.y(:,end))
    end
    
else
    %set warning if derivative is not close to 0
    deriv = system(1,ODE_sol.y(:,end)');
    max_deriv = max(abs(deriv));
    if max_deriv > 0.1
        warning('Maximum derivative at end timepoint is %g. System may not be integrated to steady-state\n',max_deriv)
    end
end

%Compute error in site balance
err_site_bal = 1-sum(ODE_sol.y(:,end));

if log_val == true
    fprintf('\nError in site balance: %g\n',err_site_bal)
end

if abs(err_site_bal) >= 0.01
    error('Site balance has too much error after integration. Try solving system of DAEs or using smaller ODE tolerances')
end

%output relevant data to sol structure
if QSS == true
    sol = struct('theta',ODE_sol.y(:,end));
else
    sol = struct('t',ODE_sol.x,'theta',ODE_sol.y);
end

end