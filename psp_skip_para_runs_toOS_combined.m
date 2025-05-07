clear all
close all;
clc;

addpath(fullfile(pwd, 'methods', 'src','Utility'));
add_project_path()

warning("off");
batch_id = 30000;    % 3000 (20popsize), 3(10D) 
trg_size = 5000;                              
gen_gap =  10;        

strategy1.name = "net2skipLLsearch";
strategy1.handle = '@(x, net,skip_flag)net_for_periodLL_search(x, prob, Params, net, 1, skip_flag)';% 1-IGD 2-HV

strategy2.name = "net4XLPS";
strategy2.handle = '@(x)netLL_search(x, prob, Params, net)';  % 1-IGD 2-HV

strategy3.name = "baseline";
strategy3.handle = '@(x)baselineLL_search(x, Params, prob, 1)';  % 1-IGD 2-HV

UL_termination = 1; % 

prob_str = {'DS1(10, 1)', 'DS2(10, 1)', 'DS3(10, 1)', 'DS1(10, -1)', 'DS2(10, -1)', 'DS3(10, -1)', 'TP1()', 'TP2(14)','DS4(5,4)', 'DS5(5, 4)',...
    'DS1(4, 1)', 'DS2(4, 1)', 'DS3(4, 1)', 'DS1(4, -1)', 'DS2(4, -1)', 'DS3(4, -1)', 'TP2(4)','DS4(3,2)', 'DS5(3, 2)',...
    'DS1(2, 1)', 'DS2(2, 1)', 'DS3(3, 1)', 'DS1(2, -1)', 'DS2(2, -1)', 'DS3(3, -1)',  'TP2(2)','DS4(2,1)', 'DS5(2, 1)',
    };  % Experiments with OS has DS4D, DS5D


[onerun_parameters, num_runs] = generate_internalComparison_parameter1_list(prob_str, 21, trg_size, gen_gap, strategy1, strategy2, strategy3, batch_id, UL_termination);

tic;
parfor i = 1:num_runs
    % arguments: prob_str, strategy, seed, UL termination_criterion(1-IGD, 2-HV), internal_comparison
    generation_assisted_BLMO_pspr3('prob_str', onerun_parameters(i).problem_str, ...
        'strategy', onerun_parameters(i).strategy, ...
        'seed', onerun_parameters(i).seed,...
        'termination_criterion', onerun_parameters(i).UL_termination,...
        'batch_ID', onerun_parameters(i).batch_ID,...
        're_evaluation', 1,...
        'trgsize_control', onerun_parameters(i).trg_size, ...
        'gen_gap', onerun_parameters(i).gen_gap);
end
toc;


%---------batch_id=3 -----
batch_id = 3;    % 3000 (20popsize), 3(10D) 
trg_size = 5000;    % [500, 1000, 2000, 5000];         % parameters experiments remember to change stopping condition to IGD                          
gen_gap =  10;      % [5, 10, 15, 20, inf];  

strategy1.name = "net2skipLLsearch";
strategy1.handle = '@(x, net,skip_flag)net_for_periodLL_search(x, prob, Params, net, 1, skip_flag)';% 1-IGD 2-HV

strategy2.name = "net4XLPS";
strategy2.handle = '@(x)netLL_search(x, prob, Params, net)';  % 1-IGD 2-HV

strategy3.name = "baseline";
strategy3.handle = '@(x)baselineLL_search(x, Params, prob, 1)';  % 1-IGD 2-HV

UL_termination = 1; % 


prob_str = {'DS1(10, 1)', 'DS2(10, 1)', 'DS3(10, 1)', 'DS1(10, -1)', 'DS2(10, -1)', 'DS3(10, -1)', 'TP1()', 'TP2(14)','DS4(5,4)', 'DS5(5, 4)',...
    'DS1(4, 1)', 'DS2(4, 1)', 'DS3(4, 1)', 'DS1(4, -1)', 'DS2(4, -1)', 'DS3(4, -1)', 'TP2(4)','DS4(3,2)', 'DS5(3, 2)',...
    'DS1(2, 1)', 'DS2(2, 1)', 'DS3(3, 1)', 'DS1(2, -1)', 'DS2(2, -1)', 'DS3(3, -1)',  'TP2(2)','DS4(2,1)', 'DS5(2, 1)',
    };  % Experiments with OS has DS4D, DS5D



[onerun_parameters, num_runs] = generate_internalComparison_parameter1_list(prob_str, 21, trg_size, gen_gap, strategy1, strategy2, strategy3, batch_id, UL_termination);

tic;
parfor i = 1:num_runs
    % arguments: prob_str, strategy, seed, UL termination_criterion(1-IGD, 2-HV), internal_comparison
    generation_assisted_BLMO_pspr3('prob_str', onerun_parameters(i).problem_str, ...
        'strategy', onerun_parameters(i).strategy, ...
        'seed', onerun_parameters(i).seed,...
        'termination_criterion', onerun_parameters(i).UL_termination,...
        'batch_ID', onerun_parameters(i).batch_ID,...
        're_evaluation', 1,...
        'trgsize_control', onerun_parameters(i).trg_size, ...
        'gen_gap', onerun_parameters(i).gen_gap);
end
toc;


remove_project_path();
