clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_BLOCparallel_asWell = 1;
fun_choices = {'ackley', 'griewank', 'rosenbrock', 'rastrigin'};  % index 1, 2, 3, 4
input_vals = [1];  % example: use 1 = Ackley only
M = 100;         % M = 100
Num_exp = 10;  
maxtime_each = 5*3600; % dont change
N = M*(M-1)/2;



for i = 1:length(input_vals)
    All_funvalNcomptimes = nan(Num_exp, 2);
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
        
        %%% BLOC: parallel
        
        if(isempty(gcp('nocreate')))
            parpool;
        end
        [C_BLOCparallel, fval_BLOCparallel, comp_time_BLOCparallel] = BLOCparallel(objFun, C0, 'MaxTime', maxtime_each);
        
        All_funvalNcomptimes(ii,:) = [fval_BLOCparallel, comp_time_BLOCparallel];
    end
    
    
    % timestamp = datestr(now, 'yyyymmdd_HHMM');
    %filename = ['Outputs/All_highdim_funvalNcomptimes_' which_fun '_M_' num2str(M) '_reps_' num2str(Num_exp) '_' timestamp '.csv'];
    
    filename = ['Outputs/All_highdim_funvalNcomptimes_' which_fun '_M_' num2str(M) '_reps_' num2str(Num_exp) '.csv'];
    
    writematrix(All_funvalNcomptimes, filename);

end
