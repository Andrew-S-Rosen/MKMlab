function [inputs,options] = parse_inputs(filename,T,gases_str,provided_P,optionals)
%parses input parameters for run_mkm

p = inputParser;
addRequired(p,'filename',@ischar);
addRequired(p,'T',@isnumeric);
addRequired(p,'gases_str',@iscellstr);
addRequired(p,'provided_P',@isnumeric);
addOptional(p,'sites_str',default_params.sites_str,@iscellstr);
addOptional(p,'provided_y0',default_params.provided_y0,@isnumeric);
addParameter(p,'tspan',default_params.tspan,@isnumeric);
addParameter(p,'ODE_algorithm',default_params.ODE_algorithm,@ischar);
addParameter(p,'DAE',default_params.DAE,@islogical);
addParameter(p,'ss_deriv_tol',default_params.ss_deriv_tol,@isnumeric);
addParameter(p,'ss_diff_tol',default_params.ss_diff_tol,@isnumeric);
addParameter(p,'ODE_rel_tol',default_params.ODE_rel_tol,@isnumeric);
addParameter(p,'ODE_abs_tol',default_params.ODE_abs_tol,@isnumeric);
addParameter(p,'root_algorithm',default_params.root_algorithm,@ischar);
addParameter(p,'root_max_fun_eval',default_params.root_max_fun_eval,@isnumeric);
addParameter(p,'root_fun_tol',default_params.root_fun_tol,@isnumeric);
addParameter(p,'root_max_iter',default_params.root_max_iter,@isnumeric);
addParameter(p,'root_step_tol',default_params.root_step_tol,@isnumeric);
addParameter(p,'log',default_params.log,@islogical);
addParameter(p,'TOF_species',default_params.TOF_species,@ischar);
parse(p,filename,T,gases_str,provided_P,optionals{:});
in = p.Results;

%basic error checking
if length(in.gases_str) ~= length(in.provided_P)
    error('Array of partial pressures is not same length as array of gas species')
end
if in.T <= 0
    error('Temperature is in Kelvins (must be > 0)')
end
if sum(in.provided_P < 0) > 0
    error('Partial pressures must be >= 0')
end
if sum(in.provided_y0) ~= 1
    error('Sum of initial coverages does not equal 1')
end

%define structures of input arguments
inputs = struct('T',T,'gases',{gases_str},'provided_P',provided_P,...
    'sites_str',{in.sites_str},'provided_y0',in.provided_y0);

ODE_options = struct('tspan',in.tspan,'ODE_algorithm',...
    in.ODE_algorithm,'DAE',in.DAE,'ss_deriv_tol',...
    in.ss_deriv_tol,'ss_diff_tol',in.ss_diff_tol,...
    'ODE_rel_tol',in.ODE_rel_tol,'ODE_abs_tol',in.ODE_abs_tol);

optim_options = optimoptions(@fsolve,'Display','off',...
    'FiniteDifferenceType','central','Algorithm',in.root_algorithm,...
    'MaxFunctionEvaluations',in.root_max_fun_eval,...
    'FunctionTolerance',in.root_fun_tol,'MaxIterations',...
    in.root_max_iter,'StepTolerance',in.root_step_tol);

options = struct('log',in.log,'TOF_species',in.TOF_species,...
    'ODE_options',ODE_options,'optim_options',optim_options);

end

