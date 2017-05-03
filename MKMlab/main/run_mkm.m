function [soln_struct,rxn_struct] = run_mkm(filename,T,gases_str,provided_P,varargin)
%given input file and experimental inputs, construct and solve a microkinetic model

%parse inputs
[inputs, options] = parse_inputs(filename,T,gases_str,provided_P,varargin);

%unpack structure
log = options.log;
sites_str = inputs.sites_str;
provided_y0 = inputs.provided_y0;
tspan = options.ODE_options.tspan;
DAE = options.ODE_options.DAE;
optim_options = options.optim_options;
TOF_species = options.TOF_species;

%print experimental inputs
if log == true
    fprintf('*******************************************\nExperimental inputs:\n\n')
    fprintf('Temperature: %g K\n',T)
    disp(table(gases_str',provided_P','VariableNames',{'Gases','Partial_Pressures'}))
    fprintf('\n')
    disp(table(sites_str',provided_y0,'VariableNames',{'Species','Initial_Coverage'}))
end

%get rate of production for each site species
rxn_struct = get_net_rate_production(filename,inputs,options);

%unpack structure
site_species = rxn_struct.site_species;
eqn_handle_str = rxn_struct.eqn_handle_str;

%solve the system of equations for coverages
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

%solve for coverages using ODE/inputs.DAE solver
soln_struct = ODE_solver(eqn_handle_str,y0,options);

%solve for roots using fsolve if QSS is implied using ODE_sol as initial guess
if isempty(tspan) == true && license('test','Optimization_Toolbox') == true
    soln_struct = root_solver(soln_struct.theta',rxn_struct,options);
elseif isempty(tspan) == true && license('test','Optimization_Toolbox') == false
    warning('Optimization toolbox not installed. Steady-state coverages determined from integration only')
end

theta = soln_struct.theta;

%print coverage results
if log == true
    if isempty(tspan) == true
        fprintf('*******************************************\nResults at steady-state:\n\n')
        disp(table(site_species',theta','VariableNames',{'Species','Coverage'}))
    else
        fprintf('*******************************************\nResults at end-time:\n\n')
        disp(table(site_species',theta(:,end),'VariableNames',{'Species','Coverage'}))
    end
end

%add entries to structures
[soln_struct.site_species] = deal(site_species);
[rxn_struct.y0] = deal(y0);

%solve for TOF and add to soln_struct structure
if isempty(TOF_species) == false
    TOF_eqn_handle_str = rxn_struct.TOF_eqn_handle_str;
    if isempty(tspan) == true
        TOF = calculate_TOF(TOF_eqn_handle_str,theta');
    else
        TOF = calculate_TOF(TOF_eqn_handle_str,theta(:,end)');
    end
    [soln_struct.TOF] = deal(TOF);
end

end