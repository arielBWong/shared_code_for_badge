function bilevel_LP_test()
clc;
clear;
% This function is used to generate PF/PS for UL problems
% its LL problem is solved by using Linear programming
seed = 1;
rng(seed,"twister");
addpath(fullfile(pwd, 'methods', 'src','Utility'));
add_project_path();

% prob = real_world1();
prob = real_world2();

Params = load_parameters(prob, 0); % internal comparison


saved_result = cell(1, 6);
saved_result{1} = [];  % ULcount
saved_result{2} = [];  % LLcount;
saved_result{3} = [];  % generation igd (w.r.t. archive);  One more data point to calculate normalized igd
saved_result{4} = [];  % generation hv (w.r.t. archive);   One more data point to calculate normalized hv
saved_result{5} = {};  % generation data
saved_result{6} = {};  % all evaluated data (instead of generation, they are children)
saved_result{7} = {};  % archive of generation
saved_result{8} = {};  % last archive with full XL
saved_result{9} = 0;   % time recored
saved_result{10} = []; % LL idg of final last archive


tmp_Trgsolutions = solutions(); % temporary variable for training data collection
trgData_solutions = solutions();
Active_pop = init_active_population();

% Init UL population
XU = Initialize_UL_Pop(prob, Params);

gen_llcount = 0;
gen_ulcount = 0;

Number_OneLLsolutions = 0;
trgsize_control = 5000;

f1 = figure();

for k = 1: Params.UL_popsize 
    fprintf("[INFO] process individual %k \n", k);

    [XL, FL, FLC, LLcount, feasible_exist, xl_one, xfeasible] = LL_LinProg(XU(k,:), prob);

    UL_visualization(f1, XU(k,:), XL, xl_one, xfeasible,prob);

    if size(XL, 1)==1
        Number_OneLLsolutions = Number_OneLLsolutions + 1;
    end

    % fprintf('[INFO] LL return solutions size %d in ind %d \n', size(XL,1), k);
    gen_llcount = gen_llcount + LLcount;
   %  [XL, FL, FLC] = LL_postprocess(XL, FL, FLC); % XL sorting with f1 value

    % special case handling: LL returned solution is infeasible
    if size(FLC, 1) == 1 &&  FLC(1)>0 
        FU = ones(1, prob.ul_nobj) * Inf; % no evaluation,  use position holder
        FC = ones(1, prob.ul_ncon) * Inf;
    else 
        % Evaluating the UL objective
        [FU, FC] = UL_Evaluate(XU(k,:), XL, prob);
        gen_ulcount = gen_ulcount + size(XL, 1);
    end

    fprintf('[INFO] processind the %d th individual \n', k);
    if feasible_exist
        % accumulate training data
        trgData_solutions.add(XU(k,:), XL, FU, FL, FC, FLC);
    end

    % save active population
    Active_pop = struct_population_expand(repmat(XU(k,:),size(XL,1),1), FU, FC, FLC, XL, FL, Active_pop);
    tmp_Trgsolutions.add(XU(k,:), XL, FU, FL, FC, FLC);
end

saved_result{2} = [gen_llcount]; 
saved_result{1} = [gen_ulcount];

% post process after one generation
saved_result{5} = [saved_result{5}, tmp_Trgsolutions]; % generation pop data
saved_result{6} = [saved_result{6}, tmp_Trgsolutions]; % generation child data

% the following line should be before clear
tmp_Trgsolutions2 = solutions();
tmp_Trgsolutions2.copy(tmp_Trgsolutions);
tmp_Trgsolutions2.nd_sort;
saved_result{7} = [saved_result{7}, tmp_Trgsolutions2];

tmp_Trgsolutions3 = solutions();
tmp_Trgsolutions3.copy(tmp_Trgsolutions);
tmp_Trgsolutions3.nd_sort('keep all FL');
saved_result{8} = [saved_result{8}, tmp_Trgsolutions3];

clear tmp_Trgsolutions
clear tmp_Trgsolutions2
clear tmp_Trgsolutions3

% Determine what strategy to run, and what  type of experiments to run
funchn = [];
strategy = "linprog";

%---creating net
hiddenSizes = max(length(prob.ul_bu)+1, length(prob.ll_bu))*2;
net = net_ffn(hiddenSizes, trgsize_control);
net.data_accumulation(cat(1, trgData_solutions.XUE{:}), cat(1, trgData_solutions.XL{:}), cat(1, trgData_solutions.FLC{:}));

% tic;
net.train([prob.ul_bu, 1],[prob.ul_bl, 0], prob.ll_bu, prob.ll_bl, false);

[saved_result, net] = UL_Evolution_Linprog_for_LL(Params, prob, Active_pop, strategy, funchn, 0, net, saved_result, 10, seed, Number_OneLLsolutions);


remove_project_path();
end

function pop = init_active_population()
pop.XU=[]; pop.FU=[]; pop.FC=[];
pop.XL=[]; pop.FL=[]; pop.FLC=[];
pop.LLcount = [];
end


function UL_visualization(f1, xu, XLs, xl, xfeasible, prob)

num = size(XLs, 1);
XU = repmat(xu, num, 1);

[FU, FC] = prob.evaluate_u(XU, XLs);
scatter(FU(:, 1), FU(:, 2), 10, 'b'); hold on;

[fu, fc] = prob.evaluate_u(xu, xl);
scatter(fu(1), fu(2), 20, 'r', 'filled');

num = size(xfeasible, 1);
XU = repmat(xu, num, 1);
[FU, FC] = prob.evaluate_u(XU, xfeasible);

scatter(FU(:,1), FU(:, 2), 40, 'green');


end

