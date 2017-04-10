classdef default_params
    %default parameters
    
    properties (Constant)
        ODE_algorithm = 'ode15s';
        DAE = false;
        ss_deriv_tol = 1E-6;
        ss_diff_tol = 1E-10;
        ODE_rel_tol = 1E-3;
        ODE_abs_tol = 1E-6;
        root_algorithm = 'levenberg-marquardt';
        root_max_fun_eval = 10000;
        root_fun_tol = 1E-6;
        root_max_iter = 10000;
        root_step_tol = 1E-6;
        log = false;
        TOF_species = [];
        tspan = [];
        provided_y0 = 1;
        sites_str = {'*'};
    end
    
end

