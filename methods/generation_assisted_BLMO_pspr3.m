function generation_assisted_BLMO_pspr3(varargin)
%---
% This version is used to combine 2 versions (non-causal origin + causal VAA ) together
% create argument parser
p = inputParser;
addParameter(p, 'prob_str', 'DS1(10, 1)');
addParameter(p, 'strategy', 'net4XLNDinit');
addParameter(p, 'seed', 1);
addParameter(p, 'termination_criterion', 2);
addParameter(p, 'batch_ID', 99);
addParameter(p, 're_evaluation', 1);
addParameter(p, 'trgsize_control', 5000);
addParameter(p, 'gen_gap', 10);
parse(p, varargin{:});
% parse input parameter 
prob_str = p.Results.prob_str;
strategy = p.Results.strategy;
seed = p.Results.seed;
termination_criterion = p.Results.termination_criterion;
batch_ID = p.Results.batch_ID;
re_evaluation = p.Results.re_evaluation;
trgsize_control = p.Results.trgsize_control;
gen_gap = p.Results.gen_gap;
%----------------------


warning("off");
prob = eval(prob_str);
rng(seed+1, 'twister');
% ------------
% checking whether exits
% batch_folder = strcat('Batch_', num2str(batch_ID));
% batch_folder = fullfile(pwd, batch_folder);
% filename = sprintf("%s_%s_trgsize_%d_frequence_%d_seed_%d_saved_result.mat",prob.name, strategy.name, trgsize_control, gen_gap, seed);
% savename = fullfile(batch_folder, filename);
% if exist(savename, "file" )
%     return
% end
% --------
pf_file = fullfile(pwd, 'problems', strcat(prob.name, '_ULPF1025.mat'));
load(pf_file);
pf_ref = 1.1 * max(pf, [], 1);
if strcmp(prob.name, 'TP1') pf_ref = [-0.9, 0.1];  end % for TP1 only
%-----------
[r, xl_causal_fl] = LL_variable_test(prob);

%r = false;
%-- check for the reproductivity
isSubstring = contains(strategy.name, "net2skipLLsearch"); % if strategy is PSP then possible switch
if isSubstring
    if ~r
        % rng(seed+1, 'twister');
        strategy.name = "net2skipLLsearch";
        strategy.handle = '@(x, net,skip_flag)net_for_periodLL_search(x, prob, Params, net, 2, skip_flag)';
    else
        strategy.name = "net2skipLLsearch_causal";
        strategy.handle = '@(x, net,skip_flag, xl_causal_fl, r)net_for_periodLL_search(x, prob, Params, net,2, skip_flag, xl_causal_fl, r)'; % 1-IGD, 2-HV
    end
end

if batch_ID < 4
    internal_comparison = true;
else
    internal_comparison = false;
end

Params = load_parameters(prob, internal_comparison);
fprintf('[INFO] Problem: %s, batch %d, seed %d,  starts  \n', prob.name, batch_ID, seed);

%------------------
nx = size(prob.ll_bl, 2); % for saving results purpose
if strcmp(prob.name, 'DS3') || strcmp(prob.name, 'DS3D')
    if nx == 3
        nx = 2;
    end
end

%======= check existing result ========
% batch_folder = strcat('Batch_', num2str(batch_ID));
% batch_folder = fullfile(pwd, batch_folder);
% filename = sprintf("%s_%s_trgsize_%d_frequence_%d_seed_%d_Nl%d_saved_result.mat",prob.name, strategy.name, trgsize_control, gen_gap, seed, nx);
% % filename = sprintf("%s_%s_trgsize_%d_frequence_%d_seed_%d_saved_result.mat",prob.name, strategy.name, trgsize_control, gen_gap, seed);
% savename = fullfile(batch_folder, filename);
% if exist(savename, 'file')
%     fprintf('[INFO] Problem: %s, batch %d, seed %d,  existing... return  \n', prob.name, batch_ID, seed);
%     return;
% end
% %=============================================


saved_result = cell(1, 6);
saved_result{1} = []; % ULcount
saved_result{2} = []; % LLcount;
saved_result{3} = []; % generation igd (w.r.t. archive);  One more data point to calculate normalized igd
saved_result{4} = []; % generation hv (w.r.t. archive);   One more data point to calculate normalized hv
saved_result{5} = {}; % generation data
saved_result{6} = {}; % all evaluated data (instead of generation, they are children)
saved_result{7} = {}; % archive of generation
saved_result{8} = {}; % last archive with full XL
saved_result{9} = 0;  % time recored
saved_result{10} = []; % LL idg of final last archive
saved_result{11} = -1;
tic;

% -----------------------
tmp_Trgsolutions = solutions(); % temporary variable for training data collection
Active_pop = init_active_population();

% Init UL population
XU = Initialize_UL_Pop(prob, Params);

gen_llcount = 0;
gen_ulcount = 0;


% === only for add one more point in convergence plot of DS4-5
FU0 = []; % FU0 is from FU values before UL additional search for VAA
% ===

for k = 1: Params.UL_popsize
    if ~r
        [XL, FL, FLC, LLcount] = LL_Evolution(XU(k,:), Params, prob, [], false, 0);
    else
        [XL, FL, FLC, LLcount, ULcount, fu0] = LL_Evolution_Causal(XU(k,:), Params, prob, [], false, 0, xl_causal_fl, r);
        gen_ulcount = gen_ulcount + ULcount;
        FU0 = [FU0; fu0];
    end
        
    % fprintf('[INFO] LL return solutions size %d in ind %d \n', size(XL,1), k);
    gen_llcount = gen_llcount + LLcount;
    [XL, FL, FLC] = LL_postprocess(XL, FL, FLC); % XL sorting with f1 value

    % special case handling: LL returned solution is infeasible
    if size(FLC, 1) == 1 &&  FLC(1)>0 
        FU = ones(1, prob.ul_nobj) * Inf; % no evaluation,  use position holder
        FC = Inf;
    else 
        % Evaluating the UL objective
        [FU, FC] = UL_Evaluate(XU(k,:), XL, prob);
        gen_ulcount = gen_ulcount + size(XL, 1);
    end
    % save active population
    Active_pop = struct_population_expand(repmat(XU(k,:),size(XL,1),1), FU, FC, FLC, XL, FL, Active_pop);
    tmp_Trgsolutions.add(XU(k,:), XL, FU, FL, FC, FLC);
end


% === only for add one more point in convergence plot of DS4-5
if size(FU0, 1) == 0
    saved_result{11} = inf;
else
    [FU0] = keep_nd(FU0);
    saved_result{11} = mean(min(pdist2(pf,FU0),[],2));
end
%===============================================================

saved_result{2} = [gen_llcount]; 
saved_result{1} = [gen_ulcount];

% post process after one generation
saved_result{5} = [saved_result{5}, tmp_Trgsolutions]; % generation pop data
saved_result{6} = [saved_result{6}, tmp_Trgsolutions]; % generation child data


% nested strategy does not accumulate training data
if ~strcmp(strategy.name, 'baseline') || ~strcmp(strategy.name, 'set_value')
    hiddenSizes = max(length(prob.ul_bu)+1, length(prob.ll_bu))*2;
    net = net_ffn(hiddenSizes, trgsize_control);
    net.data_accumulation(cat(1, tmp_Trgsolutions.XUE{:}), cat(1, tmp_Trgsolutions.XL{:}),  cat(1, tmp_Trgsolutions.FLC{:}));
    % tic;
    net.train([prob.ul_bu, 1],[prob.ul_bl, 0], prob.ll_bu, prob.ll_bl, false);
    % toc;
end

if strcmp(strategy.name, 'set_value')
    net = fiber_models_create(Active_pop);
end

% % %----figure  LL prediction example
% LL_prediction_visualization(1, prob, net, Params);
% return;
%---------------------

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
funchn = eval(strategy.handle);

% if strcmp(strategy.name, 'net4XLPS')
%     % run accuracy experiments
%     accuracy_experiment(batch_ID, prob, strategy, net, seed, 0);
%     return
% end
% 
% if strategy.name == "net2skipLLsearch" % if 5000 points are not accumulated, especially when D=2
%     if size(net.Trg_LHS, 1)>4999
%         accuracy_experiment(batch_ID, prob, strategy, net, seed, 0);
%         return
%     end
% end
%----------------- -arguments: Params, prob, Active_pop, strategy, funchn, termination_criterion, additional_termination, net


if ~r
    [saved_result, net] = UL_Evolution_v2(Params, prob, Active_pop, strategy, funchn, termination_criterion, net, saved_result, gen_gap, seed);
else
    [saved_result, net] = UL_Evolution_v2causal(Params, prob, Active_pop, strategy, funchn, termination_criterion, net, saved_result, gen_gap, xl_causal_fl, r);
end
toc;

% if isempty(saved_result)
%     fprintf('return from accuracy experiment \n');
%     return
% end

saved_result{9} = saved_result{9} + toc;

%----re-evaluation 
if re_evaluation && ~strcmp(strategy.name, 'baseline')
    tmp_Trgsolutions = solutions();
    gen_llcount = 0;
    gen_ulcount = 0;
    nd = saved_result{8};
    nd_xu = nd.xus;
    for ii = 1:size(nd_xu, 1)
        xu = nd_xu(ii, :);
        if ~r
            termination_criterion = 1;
            [XL, FL, FLC, LLcount] = net_as_initNDLL_search(xu, prob, Params, net, termination_criterion);
            ULcount = 0;
        else
            termination_criterion = 1;
            [XL, FL, FLC, LLcount, ULcount] = net_as_initNDLL_search_causal(xu, prob, Params, net, termination_criterion, xl_causal_fl, r);

        end
        
        gen_llcount = gen_llcount + LLcount;

        if size(FLC, 1) == 1 &&  FLC(1) > 0   % no feasible solution
            FU = ones(1, prob.ul_nobj) * Inf;    % no evaluation,  no count, use position holder
            FC = Inf;
        else
            [FU, FC] = UL_Evaluate(xu, XL, prob);
            gen_ulcount = gen_ulcount + size(XL, 1) + ULcount;
        end
        tmp_Trgsolutions.add(xu, XL, FU, FL, FC, FLC);
    end
    saved_result{2} = [saved_result{2}, saved_result{2}(end) + gen_llcount];
    saved_result{1} = [saved_result{1}, saved_result{1}(end) + gen_ulcount];

    % save ND front of this generation
    last_nd = saved_result{7}(end);
    tmp_solutions = solutions();
    tmp_solutions.merge(last_nd);
    tmp_solutions.merge(tmp_Trgsolutions);
    tmp_solutions.nd_sort;
    saved_result{7} = [saved_result{7}, tmp_solutions]; %nd front of each generation
    clear tmp_solutions

    % save ND front of this generation
    last_nd = saved_result{8}(end);
    tmp_solutions = solutions();
    tmp_solutions.merge(last_nd);
    tmp_solutions.merge(tmp_Trgsolutions);
    tmp_solutions.nd_sort('keep all FL');
    saved_result{8} = tmp_solutions; %nd front of each generation
    clear tmp_Trgsolutions
    clear tmp_solutions  

    % Record IGD and HV for re-evaluation
    current_nd_solutions = saved_result{7}(end);
    ndFU = current_nd_solutions.FUs;
    igd = mean(min(pdist2(pf,ndFU),[],2));
    hv = Hypervolume(ndFU, pf_ref);
    saved_result{3} = [saved_result{3}, igd];
    saved_result{4} = [saved_result{4}, hv];
end


% calculate normalized UL IGD
current_nd_solutions = saved_result{7}(end);
ndFU = current_nd_solutions.FUs;
pf_file = fullfile(pwd, 'problems', strcat(prob.name, '_ULPF1025.mat'));
load(pf_file);
ul_ideal = min(pf, [], 1);
ul_nadir = max(pf, [], 1);
norm_pf = (pf - ul_ideal) ./ (ul_nadir - ul_ideal);
normpf_ref = ones(1, 2) .* 1.1;
norm_nd = (ndFU - ul_ideal) ./ ((ul_nadir - ul_ideal));
norm_igd = mean(min(pdist2(norm_pf, norm_nd),[],2));
norm_hv = Hypervolume(norm_nd, normpf_ref);
saved_result{3} = [saved_result{3}, norm_igd];
saved_result{4} = [saved_result{4}, norm_hv];


% calculate LL IGD
nd = saved_result{8};
nd_xu = nd.xus;
fprintf("[INFO] Calculate LL igd for %d UL solutions \n", size(nd_xu, 1));
tic
ll_igdarray = [];
for ii = 1:size(nd_xu, 1)
    xu = nd_xu(ii, :);
    FL = nd.FL{ii};
    pf_FL = prob.PF_LL(1025, xu);

    ll_ideal = min(pf_FL, [], 1);
    ll_nadir = max(pf_FL, [], 1);

    norm_FL = (FL - ll_ideal) ./(ll_nadir - ll_ideal);
    norm_pf = (pf_FL - ll_ideal) ./ (ll_nadir - ll_ideal);
    ll_igd = mean(min(pdist2(norm_pf, norm_FL), [], 2));
    ll_igdarray = [ll_igdarray,  ll_igd];
end
toc;
LL_igd_time = toc;
fprintf("[INFO] Calculate LL igd for %d UL solutions takes %0.4f seconds\n", size(nd_xu, 1), LL_igd_time);

saved_result{10} = mean(ll_igdarray);

% end

% ======= parameter experiments save ==============
batch_folder = strcat('Batch_', num2str(batch_ID));
batch_folder = fullfile(pwd, batch_folder);
n = exist(batch_folder);
if n ~= 7
    mkdir(batch_folder)
end

strategy_term = strategy.name;

if contains(strategy_term, "causal") % for unify names
    strategy_term = 'net2skipLLsearch';
end

filename = sprintf("%s_%s_trgsize_%d_frequence_%d_seed_%d_Nl%d_saved_result.mat",prob.name, strategy_term, trgsize_control, gen_gap, seed, nx); % using Nl, is for internal comparison where variable size is changing
savename = fullfile(batch_folder, filename);
save(savename, 'saved_result');
% ======= parameter experiments save ==============

% record_results(strategy, seed, prob, saved_result, batch_ID, OS_lastEval, Params);
fprintf('[INFO] Problem: %s, batch %d, seed %d,  UL finish  \n', prob.name, batch_ID, seed);
end

function pop = init_active_population()
pop.XU=[]; pop.FU=[]; pop.FC=[];
pop.XL=[]; pop.FL=[]; pop.FLC=[];
pop.LLcount = [];
end


function[F] = keep_nd(F)
[fronts, ids, ~] = nd_sort(F, (1:size(F, 1))');
front_idx = [];
for j = 1:size(fronts,2)
    front_idx = [front_idx; ones(length(fronts(j).f), 1)*j];
end
F = F(ids, :);
id1 = front_idx == 1;
F = F(id1, :);

end

% 
% 
% function skip_flag = skip_baseline(strategy, batch_ID, prob, seed, nx)
% skip_flag = false;
% batch_folder = strcat('Batch_', num2str(batch_ID));
% batch_folder = fullfile(pwd, batch_folder);
% n = exist(batch_folder);
% if n ~= 7
%     mkdir(batch_folder)
% end
% 
% outfoldername = strcat('result_folder_', num2str(nx));
% savefolder_path = fullfile(batch_folder, outfoldername);
% n = exist(savefolder_path);
% if n ~= 7
%     mkdir(savefolder_path)
% end
% 
% strategy_path = fullfile(savefolder_path, strategy.name);
% n = exist(strategy_path);
% if n ~= 7
%     mkdir(strategy_path)
% end
% 
% problem_path = fullfile(strategy_path, prob.name);
% n = exist(problem_path);
% if n ~= 7
%     mkdir(problem_path)
% end
% 
% %-----------------
% filename = strcat('UL_archive_seed_', num2str(seed),'.mat');
% savename = fullfile(problem_path, filename);
% 
% if exist(savename)
%     skip_flag = true;
% end
% 
% end