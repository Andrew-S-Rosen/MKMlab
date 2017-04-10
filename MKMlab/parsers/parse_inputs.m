function input_results = parse_inputs(filename,T,gases_str,provided_P,optionals)
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
input_results = p.Results;
end