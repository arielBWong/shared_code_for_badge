function [saved_result, net] = UL_Evolution_v2causal(Params, prob, generation_pop, strategy, funchn, termination_criterion, net, saved_result, gen_gap, xl_causal_fl, r)
%-
if termination_criterion == 1
    termination_length = 5;
else
    termination_length = 10;
end
termination_vector = ones(3, termination_length) .* Inf;    % IGD use 3 rows, HV use 1 row
pf_file = fullfile(pwd, 'problems', strcat(prob.name, '_ULPF1025.mat'));
load(pf_file);
pf_ref = 1.1 * max(pf, [], 1);
if strcmp(prob.name, 'TP1') pf_ref = [-0.9, 0.1];  end % for TP1 only

% adding saved result
current_nd_solutions = saved_result{7}(1);
ndFU = current_nd_solutions.FUs;
igd = mean(min(pdist2(pf,ndFU),[],2));
hv = Hypervolume(ndFU, pf_ref);
saved_result{3} = [saved_result{3}, igd];
saved_result{4} = [saved_result{4}, hv];

vis = false;
if vis
    f1 = figure(1);
    UL_visualization(prob,saved_result, f1, 1);
end

% UL evolution
param.popsize = Params.UL_popsize; % wrapper for using gsolver method
childXU = generate_child_DE(prob.ul_bl, prob.ul_bu, unique(generation_pop.XU, 'rows', 'stable'), param);
% Remove repeated ones
childXU = unique(childXU, 'rows', 'stable');
childXU = remove_repeated_solution(unique(generation_pop.XU, 'rows', 'stable'), childXU);

for gen = 1 : (Params.N_gen_UL - 1) % initialization caused one generation used

    fprintf(['[INFO] Problem: %s, size %d,  UL generation %d  \n'], prob.name, length(prob.ll_bl), gen+1);
    child_pop.XU = []; child_pop.FU = []; child_pop.FC = []; child_pop.XL = []; child_pop.FL = []; child_pop.FLC = []; % intermediate variable for the process of selecting children
    tmp_Trgsolutions = solutions(); % temporary variable for training data collection

    genCount_UL = 0;
    genCount_LL = 0;
    for i = 1:size(childXU , 1)
        % fprintf('[INFO] Child evaluation %d  \n', i);
        [XL, FL, FLC, LLcount, ULcount] = evaluate_ULChild(strategy, funchn, childXU(i, :), net, gen+1, gen_gap, xl_causal_fl, r);

        genCount_LL = genCount_LL + LLcount;
        genCount_UL = genCount_UL + ULcount;
        if size(FLC, 1) == 1 &&  FLC(1) > 0   % no feasible solution
            FU = ones(1, prob.ul_nobj) * Inf;    % no evaluation,  no count, use position holder
            FC = Inf;
        else
            [FU, FC] = UL_Evaluate(childXU(i,:), XL, prob);
            genCount_UL = genCount_UL + size(XL, 1);
        end
       
        tmp_Trgsolutions.add(childXU(i, :), XL, FU, FL, FC, FLC);
        child_pop = struct_population_expand(repmat(childXU(i, :), size(XL, 1), 1), FU, FC, FLC, XL, FL, child_pop);
    end
    saved_result{1} = [saved_result{1},saved_result{1}(end)+genCount_UL];
    saved_result{2} = [saved_result{2}, saved_result{2}(end)+genCount_LL];
    
    saved_result{6} = [saved_result{6}, tmp_Trgsolutions];      % generation child data

    % decide whether child should be used for update network training
    if update_net_decision(strategy.name, gen+1, net, gen_gap)
        net.data_accumulation(cat(1, tmp_Trgsolutions.XUE{:}), cat(1, tmp_Trgsolutions.XL{:}),  cat(1, tmp_Trgsolutions.FLC{:}));
        fprintf('[INFO] update training data and retrain network %d  \n',size(net.Trg_LHS, 1));
        net.train([prob.ul_bu, 1],[prob.ul_bl, 0], prob.ll_bu, prob.ll_bl, false);
    else
        fprintf('[INFO] network size: %d \n', size(net.Trg_LHS, 1));
    end

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
    saved_result{8} = tmp_solutions; % nd front of each generation
    clear tmp_Trgsolutions
    clear tmp_solutions

    % combine two population
    generation_pop = struct_population_expand(child_pop.XU, child_pop.FU, child_pop.FC,...
        child_pop.FLC, child_pop.XL, child_pop.FL, generation_pop);
    % Sort 2N population
    [generation_pop, front_idx] = population_sort(generation_pop);
    nd_idx_binary = front_idx == 1;

    % Extract nd xu solutions
    xu_candidate = generation_pop.XU(nd_idx_binary, :);
    % Unique nd xu for DSS selection
    unique_xuCandidate = unique(xu_candidate, 'rows', 'stable');
    select_candidateID = find(nd_idx_binary == 1);

    % SSD selection
    if size(unique_xuCandidate, 1) > Params.UL_popsize
        sparseXID = Sparse_Selection(select_candidateID, generation_pop.XU, prob.ul_bl, prob.ul_bu, Params.UL_popsize);
        select_from_compactID = select_candidateID(sparseXID);
    else
        tmp_uniqueXU = unique(generation_pop.XU, 'rows',  'stable');
        [~, ~, ib] = intersect(tmp_uniqueXU, generation_pop.XU, 'stable', 'rows'); % unique xu in order in generation_pop.XU
        select_from_compactID = ib(1:Params.UL_popsize);
    end

    % 2N->N
    XU = generation_pop.XU(select_from_compactID, :);
    % create new Active populuation for XU
    [generation_pop, solution_pop] = new_active_pop(generation_pop, XU, Params);

    % save progress information,
    saved_result{5} = [saved_result{5}, solution_pop]; % generation pop data
    current_nd_solutions = saved_result{7}(end);
    ndFU = current_nd_solutions.FUs;
    igd = mean(min(pdist2(pf,ndFU),[],2));
    hv = Hypervolume(ndFU, pf_ref);
    saved_result{3} = [saved_result{3}, igd];
    saved_result{4} = [saved_result{4}, hv];

    if vis UL_visualization(prob,saved_result, f1, gen+1); end

    % after each generation, check whether to terminate % termination_criterion,  gen, termination_vector, verbose
    [termination_flag, termination_vector] = Termination_criterion_activate(saved_result{7}, termination_criterion, gen+1, termination_vector, false);
    if termination_flag
        break
    else
        if Params.additional_stopping % special condition only for batch 5 experiments
            if igd < prob.target_igd
                break
            end
        end
    end

    if strcmp(strategy.name, 'set_value')
        every_period = mod(gen+1, gen_gap);
        if every_period == 0
            % fprintf('[INFO] Set valued scheme re train \n')
            try
                net = fiber_models_create(generation_pop);
            catch e
                fprintf('[INFO] The identifier was: %s \n',e.identifier);
                fprintf('[INFO] There was an error! The message was: %s \n',e.message);
                net = net;
            end
        end
    end

    % offspring generation
    childXU = generate_child_DE(prob.ul_bl, prob.ul_bu, XU, param);
    childXU = unique(childXU, 'rows', 'stable');

    all_evaluated_xu = [saved_result{6}.xu];
    all_evaluated_xu = cat(1, all_evaluated_xu{:});
    childXU = remove_repeated_solution(unique(all_evaluated_xu, 'rows', 'stable'), childXU);
end

if vis close(f1); end

end


function UL_visualization(prob,saved_result, f1, gen)

pf_filename = strcat(prob.name, '_ULPF1025.mat');
pf_path = fullfile(pwd, 'problems');
pf_path = fullfile(pf_path, pf_filename);
load(pf_path);

nd = saved_result{7}(end).FUs;
igd = saved_result{3}(end);
hv = saved_result{4}(end);
pop = saved_result{5}(end).FUs;
fprintf('[INFO] IGD is %0.4f \n', igd);
fprintf('[INFO] HV is %0.4f \n', hv);
clf(f1);
plot(pf(:, 1), pf(:, 2), 'k.'); hold on;
scatter(nd(:, 1), nd(:, 2), 20, 'filled');
scatter(pop(:, 1), pop(:, 2));
t = strcat('gen ', num2str(gen), ', ', num2str(saved_result{1}), ', ', num2str(saved_result{2}));
title(t);
s1 = strcat('IGD  ', num2str(igd));
s2 = strcat('HV  ' , num2str(hv));
legend(s1, s2);

pause(0.1);

end




function [Active_pop, front_idx] = population_sort(Active_pop)
% sorting, LL infeasible is assigned Inf F and Inf C, so they will
% always be in lowest rank
% data structure convertion
% sort normally according to F and C
popwrapper.X = Active_pop.XU;
popwrapper.F = Active_pop.FU;
popwrapper.C = Active_pop.FC;

popwrapper.XL = Active_pop.XL;
popwrapper.FL = Active_pop.FL;
popwrapper.FLC = Active_pop.FLC;
% call sort function
[popwrapper, front_idx] = pop_sort(popwrapper);

% Convert wrapper back to Archive
Active_pop.XU = popwrapper.X;
Active_pop.FU = popwrapper.F;
Active_pop.FC = popwrapper.C;
Active_pop.XL = popwrapper.XL;
Active_pop.FL = popwrapper.FL;
Active_pop.FLC = popwrapper.FLC;

end


function [Active_pop, solution_pop] = new_active_pop(Active_compactpop, XU, Params)
Active_pop.XU = []; Active_pop.FU = []; Active_pop.FC = [];
Active_pop.XL = []; Active_pop.FL = []; Active_pop.FLC = [];
solution_pop = solutions();

for i = 1: Params.UL_popsize
    ia = ismember(Active_compactpop.XU, XU(i, :), 'rows');
    Active_pop.XU = [Active_pop.XU; Active_compactpop.XU(ia, :)];
    Active_pop.FU = [Active_pop.FU; Active_compactpop.FU(ia, :)];
    Active_pop.FC = [Active_pop.FC; Active_compactpop.FC(ia, :)];
    Active_pop.XL = [Active_pop.XL; Active_compactpop.XL(ia, :)];
    Active_pop.FL = [Active_pop.FL; Active_compactpop.FL(ia, :)];
    Active_pop.FLC = [Active_pop.FLC; Active_compactpop.FLC(ia, :)];

    solution_pop.add(XU(i, :), Active_compactpop.XL(ia, :),  Active_compactpop.FU(ia, :), ...
        Active_compactpop.FL(ia, :),Active_compactpop.FC(ia, :), Active_compactpop.FLC(ia, :));
end
end




