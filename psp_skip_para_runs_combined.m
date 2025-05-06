clear all
close all;
clc;

addpath(fullfile(pwd, 'methods', 'src','Utility'));
add_project_path()

warning("off");
batch_id = 20000;    % 20000, compare to state of the art experiment    2000 parameter
trg_size = 5000;     % [500, 1000, 2000, 5000];         % parameters experiments remember to change stopping condition to IGD                          
gen_gap = 10;        % [5, 10, 15, 20, inf];  

strategy.name = "net2skipLLsearch";
strategy.handle = '@(x, net,skip_flag)net_for_periodLL_search(x, prob, Params, net, 2, skip_flag)';% 1-IGD 2-HV

UL_termination = 2; % 

prob_str = {'DS1(10, 1)', 'DS2(10, 1)', 'DS3(10, 1)', 'DS1(10, -1)', 'DS2(10, -1)', 'DS3(10, -1)', 'TP1()', 'TP2(14)','DS4(5,4)', 'DS5(5, 4)','DS4D(5,4)', 'DS5D(5, 4)'}; 

[onerun_parameters, num_runs] = generate_pararun_parameter_list(prob_str, 21, trg_size, gen_gap, strategy, batch_id, UL_termination);

tic;
parfor i = 1:num_runs
    % arguments: prob_str, strategy, seed, UL termination_criterion(1-IGD, 2-HV), internal_comparison
    generation_assisted_BLMO_pspr3('prob_str', onerun_parameters(i).problem_str, ...
        'strategy', onerun_parameters(i).strategy, ...
        'seed', onerun_parameters(i).seed,...
        'termination_criterion', onerun_parameters(i).UL_termination,...
        'batch_ID', onerun_parameters(i).batch_ID,...
        're_evaluation', 0,...
        'trgsize_control', onerun_parameters(i).trg_size, ...
        'gen_gap', onerun_parameters(i).gen_gap);
end
toc;


remove_project_path();