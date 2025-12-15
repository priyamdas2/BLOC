function [C_opt, fval, comp_time] = BLOC(objFun, C0, varargin)
% BLOC: Blackbox optimization over correlation matrices
%
% Inputs:
%   objFun   - Handle to the objective function, takes M x M matrix p.d.
%   correlation matrix as input
%   C0       - Initial guess for the optimization, M x M matrix p.d.
%   correlation matrix
%   MaxTime  - maximum allowed execution time in seconds
%   MaxRuns  - number of maximum runs allowed
%   MaxIter  - number of maximum iterations allowed within each run
%   sInitial - initial step-size
%   rho      - step decay rate
%   TolFun1  - tolerance value 1 (controls step decay rate)
%   TolFun2  - tolerance value 1 (controls run restart decision)
%   phi      - minimum allowed step-size
%   DisplayUpdate - 1 = shows iteration update, 0 = doesn't show.
%   DisplayEvery  - Every that many seconds update is displayed only if
%                   DisplayUpdate is set to 1
%   PrintStepSize - 1 = Displays current step-size while updating, 0 = doesn't show. 
%   PrintSolution - 1 = Displays final solution where function is minimized, 0 = doesn't show. 

% Outputs:
%   C_opt   - The optimal solution (minimizer)
%   fval    - Function value at the optimal solution

%%%%%%%% parameters %%%%%%%%%%%%


% Parse optional inputs using inputParser or name-value style
p = inputParser;
p.addParameter('MaxTime', 3600);
p.addParameter('MaxRuns', 1000);
p.addParameter('MaxIter', 10000);
p.addParameter('sInitial', 1);
p.addParameter('rho', 2);
p.addParameter('TolFun1', 10^(-6));
p.addParameter('TolFun2', 10^(-20));
p.addParameter('phi', 10^(-20));
p.addParameter('DisplayUpdate', 1);
p.addParameter('DisplayEvery', 2);
p.addParameter('PrintStepSize', 1);
p.addParameter('PrintSolution', 0);
p.parse(varargin{:});
params = p.Results;

MaxTime        = params.MaxTime;
MaxRuns        = params.MaxRuns;
MaxIter        = params.MaxIter;
sInitial       = params.sInitial;
rho            = params.rho;
TolFun1        = params.TolFun1;
TolFun2        = params.TolFun2;
phi            = params.phi;
DisplayUpdate  = params.DisplayUpdate;
DisplayEvery   = params.DisplayEvery;
PrintStepSize  = params.PrintStepSize;
PrintSolution  = params.PrintSolution;



fprintf('========================= BLOC Starts =======================\n')
tic;

if norm(C0 - C0') > 10^-10 || any(eig(C0) < 0)
    error('C0 must be a symmetric positive semi-definite matrix');
end
if any(abs(diag(C0) - 1) > 1e-10)
    error('C0 must have 1s on its diagonal (within tolerance)');
end


Theta0 = Corr2Theta(C0);
N = length(Theta0);
M = (1 + sqrt(1 + 8*N)) / 2;
start_element_locs = ThetaGroupStartIndicesRow2onward(M);
end_element_locs = ThetaGroupEndIndicesRow2onward(M);

[domain_lb, domain_ub] = Domain_ub_lb(N);

ValAtInitialPoint  = objFun(C0);

fprintf('\n');
fprintf('=> Hmm... obj. fun. value at provided initial point is: %d. \n', ValAtInitialPoint);
fprintf('\n');
fprintf('=> Lets BLOC it up!!! \n');
fprintf('\n');

Theta_array = zeros(MaxRuns, N);
Loop_solution = zeros(MaxRuns, 1);
array_of_values = zeros(MaxIter,1);
last_toc = 0;
break_now = 0;

for iii = 1:MaxRuns
    epsilon = sInitial;
    if(iii == 1)
        theta = Theta0;
    else
        theta = transpose(Theta_array((iii-1),:));
    end
    
    for i = 1:MaxIter
        %disp([iii,i]);
        if(toc > MaxTime)
            break_now = 1;
            fprintf('=> As requested, BLOC has been terminated after %.2f seconds :( \n', MaxTime);
            fprintf('\n')
            break;
        end

        Corr = Theta2Corr(theta);
        current_lh = objFun(Corr);
        
        
        %%%% Time display %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        toc_now = toc;
        if(DisplayUpdate == 1)
            if(toc_now - last_toc > DisplayEvery)
                if(PrintStepSize == 1)
                    fprintf('=> Executing Run: %d, iter: %d, current obj. fun. value: %d, current log10(step-size): %.2f. \n', iii, i, current_lh, log10(epsilon));
                else
                     fprintf('=> Executing Run: %d, iter: %d, current obj. fun. value: %d. \n', iii, i, current_lh);
                end
                last_toc = toc_now;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Theta_in_domain = all(theta >= domain_lb & theta <= domain_ub);
        if(Theta_in_domain == 0)
            error('Theta is out of allowed domain');
        end
        
        total_lh = zeros(2*N,1);
        matrix_update_at_h = zeros(N,2*N);
         
       for location_number = 1:(2*N)         
            change_loc  = ceil(location_number/2);
            epsilon_temp = ((-1)^location_number)*epsilon;
            possibility = theta;
            value_at_position = possibility(change_loc);
            possibility_temp0 = value_at_position + epsilon_temp;
            
            if (change_loc == 1)
                possibility_temp = smooth_reflect_pi2(possibility_temp0);
            else
                if ismember(change_loc, end_element_locs)
                    possibility_temp = mod(possibility_temp0, 2*pi);
                else
                    if ismember(change_loc, start_element_locs)
                        possibility_temp = map_start_elements(possibility_temp0);
                    else
                        possibility_temp = map_middle_elements(possibility_temp0);
                    end
                end
            end

            possibility(change_loc) = possibility_temp;
            possibility_L = Theta2L(possibility);
            
            if all(diag(possibility_L) > 1e-3)
                possibility_Corr = possibility_L*possibility_L';
                total_lh(location_number) = objFun(possibility_Corr);
                matrix_update_at_h(:,location_number) = possibility;
            else
                total_lh(location_number) = current_lh;
                matrix_update_at_h(:,location_number) = theta;
            end  
       end
       
       [candidate,I] = min(total_lh);
       
       if(candidate < current_lh)
           theta = matrix_update_at_h(:,I);
       end
       array_of_values(i) =  min(candidate,current_lh);
       
       if(i > 1)
           if(abs(array_of_values(i) - array_of_values(i-1)) < TolFun1)
               if(epsilon > phi)
                   epsilon = epsilon/rho;
               else
                   break
               end
            end
        end
        
    end
    Theta_array(iii,:) = transpose(theta);
    corr = Theta2Corr(theta);
    Loop_solution(iii) = objFun(corr);
    if(iii > 1)
        old_soln = Theta_array(iii-1,:);
        new_soln = Theta_array(iii,:);
        old_soln_corr = Theta2Corr(old_soln);
        new_soln_corr = Theta2Corr(new_soln);
        if(norm(objFun(old_soln_corr) - objFun(new_soln_corr)) <TolFun2)
            break
        end
    end
end

Theta_opt = theta;
C_opt = Theta2Corr(Theta_opt);
fval = objFun(C_opt);
comp_time = toc;

% Final solution
if(PrintSolution == 1)
    fprintf('\n')
    fprintf('=> Final BLOC solution is: \n');
    disp(C_opt)
end
fprintf('\n')
fprintf('=> Obj. fun. value at BLOC minima: %d \n',fval);
fprintf('\n')
fprintf('=> Total time taken: %.4f secs.\n',comp_time);

fprintf('xxxxxxxxxxxxxxxxxxxxxx BLOC ends xxxxxxxxxxxxxxxxxxxxxxxxxx\n')

end


