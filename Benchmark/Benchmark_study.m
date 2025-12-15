clear all;

%%% Adding 'Manopt' to the path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(pwd());
% Recursively add Manopt directories to the Matlab path.
cd('manopt');
addpath(genpath(pwd()));
cd('..');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
skip_BLOCparallel = 1;
skip_activeset = 1;
skip_interiorpoint = 0;
skip_sqp = 1;

fun_choices = {'ackley', 'griewank', 'rosenbrock'};  % index 1, 2, 3
input_vals = [3];  % example
M = 5;         % M = 5, 10, 20, 50, 100
Num_exp = 3;  % Num_exp = 10 for M = 5, 10, 20, 50, 100.
maxtime_each = 3600; % dont change
N = M*(M-1)/2;
All_comp_times = nan(Num_exp, 9);
All_funvals = nan(Num_exp, 9);

for i = 1:length(input_vals)
    which_fun = fun_choices{input_vals(i)};
    
    if strcmp(which_fun, 'ackley')
        objFun = @(C) modified_ackley(C);
    elseif strcmp(which_fun, 'griewank')
        objFun = @(C) modified_griewank(C);
    elseif strcmp(which_fun, 'rosenbrock')
        objFun = @(C) modified_rosenbrock(C);
    else
        error('Unknown function specified in which_fun.');
    end
    
    
    for ii = 1:Num_exp
        fprintf('Performing experiment number: %d using function: %s\n', ii, which_fun);
        rand_seed = ii;
        C0 = randCorrMatrix(M, rand_seed);
        Theta0 = Corr2Theta(C0);
        
        %%% BLOC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [C_BLOC, fval_BLOC, comp_time_BLOC] = BLOC(objFun, C0, 'MaxTime', maxtime_each);
        
        %%% BLOC: parallel
        [C_BLOCparallel, fval_BLOCparallel, comp_time_BLOCparallel] = deal(nan);
        if(skip_BLOCparallel == 0)
            if(isempty(gcp('nocreate')))
                parpool;
            end
            [C_BLOCparallel, fval_BLOCparallel, comp_time_BLOCparallel] = BLOCparallel(objFun, C0, 'MaxTime', maxtime_each);
        end
        
        %%% fmincon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x0 = vecFromCorrMatrix(C0); % Initial values for off-diagonal correlations
        lb = -ones(size(x0));
        ub = ones(size(x0));
        obj = @(x) objFun(corrMatrixFromVec(x, M));

        %%% fmincon: active-set
        [x_opt_activeset, fval_activeset, comp_time_activeset] = deal(nan);
        if (skip_activeset == 0)
            options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'active-set');
            tic;
            [x_opt_activeset, fval_activeset] = fmincon(obj, x0, [], [], [], [], lb, ub, @correlationConstraint, options);
            comp_time_activeset = toc;
        end
        
        %%% fmincon: interior-point
        [x_opt_interiorpoint, fval_interiorpoint, comp_time_interiorpoint] = deal(nan);
        if (skip_interiorpoint == 0)
            options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point');
            tic;
            [x_opt_interiorpoint, fval_interiorpoint] = fmincon(obj, x0, [], [], [], [], lb, ub, @correlationConstraint, options);
            comp_time_interiorpoint = toc;
        end
        
        %%% fmincon: sqp
        [x_opt_sqp, fval_sqp, comp_time_sqp] = deal(nan);
        if (skip_sqp == 0)
            options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
            tic;
            [x_opt_sqp, fval_sqp] = fmincon(obj, x0, [], [], [], [], lb, ub, @correlationConstraint, options);
            comp_time_sqp = toc;
        end
        %%% Manopt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        rng(rand_seed)
        manifold = sympositivedefinitefactory(M);
        problem.M = manifold;
        problem.cost = @(X) objFun(X);
        
        %%% Manopt: barzilai-borwein
        tic;
        options = struct();
        options.maxtime = maxtime_each;
        [Xopt_barzilaiborwein, xcost_barzilaiborwein, info_barzilaiborwein] = barzilaiborwein(problem, C0, options);
        comp_time_barzilaiborwein = toc;
        
        %%% Manopt: conjugate-gradient
        tic;
        options = struct();
        options.maxtime = maxtime_each;
        [Xopt_conjugategradient, xcost_conjugategradient, info_conjugategradient] = conjugategradient(problem, C0, options);
        comp_time_conjugategradient = toc;
        
        %%% Manopt: steepest-descent
        tic;
        options = struct();
        options.maxtime = maxtime_each;
        [Xopt_steepestdescent, xcost_steepestdescent, info_steepestdescent] = steepestdescent(problem, C0, options);
        comp_time_steepestdescent = toc;
        
        %%% Manopt: trust-region
        tic;
        options = struct();
        options.maxtime = maxtime_each;
        [Xopt_trustregion, xcost_trustregion, info_trustregion] = trustregions(problem, C0, options);
        comp_time_trustregion = toc;
        
        %%% Summary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        All_comp_times(ii,:) = [comp_time_BLOC, comp_time_BLOCparallel, comp_time_activeset,...
            comp_time_interiorpoint, comp_time_sqp, comp_time_barzilaiborwein,...
            comp_time_conjugategradient, comp_time_steepestdescent, comp_time_trustregion];
        
        All_funvals(ii,:) = [fval_BLOC, fval_BLOCparallel, fval_activeset, fval_interiorpoint,...
            fval_sqp, xcost_barzilaiborwein, xcost_conjugategradient, xcost_steepestdescent,...
            xcost_trustregion];
        
    end
    
    
    timestamp = datestr(now, 'yyyymmdd_HHMM');
    
    filename = ['Outputs/All_comp_times_' which_fun '_M_' num2str(M) '_reps_' num2str(Num_exp) '_' timestamp '.csv'];
    writematrix(All_comp_times, filename);
    
    filename = ['Outputs/All_funvals_' which_fun '_M_' num2str(M) '_reps_' num2str(Num_exp) '_' timestamp '.csv'];
    writematrix(All_funvals, filename);
    
    
end

