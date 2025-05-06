function [onerun_parameters, num_runs] = generate_internalComparison_parameter1_list(prob_list, max_seed, trg_size_list, gen_gap_list, strategy_structure1,  strategy_structure2, strategy_structure3,save_batch_id, UL_termination_id)
% construct pararun parameter instances for parfor
onerun_parameters  = struct([]);
strc_id = 1;

%-------------
for dd = 1:length(trg_size_list)
    for  gg = 1:length(gen_gap_list)

        seeds = 1:max_seed;
        ns = length(seeds);
        np = length(prob_list);

        for i = 1 : np
            for j = 1 : ns
                onerun_parameters(strc_id). problem_str = prob_list{i};
                onerun_parameters(strc_id). seed = j;
                onerun_parameters(strc_id). strategy = strategy_structure1;
                onerun_parameters(strc_id). batch_ID = save_batch_id;
                onerun_parameters(strc_id). UL_termination = UL_termination_id;
                onerun_parameters(strc_id). gen_gap = gen_gap_list(gg);
                onerun_parameters(strc_id). trg_size = trg_size_list(dd);
                strc_id = strc_id + 1;
            end
        end

    end
end

% % temp
% num_runs = strc_id-1;
% return 
% %================


for dd = 1:length(trg_size_list)
    for  gg = 1:length(gen_gap_list)

        seeds = 1:max_seed;
        ns = length(seeds);
        np = length(prob_list);

        for i = 1 : np
            for j = 1 : ns
                onerun_parameters(strc_id). problem_str = prob_list{i};
                onerun_parameters(strc_id). seed = j;
                onerun_parameters(strc_id). strategy = strategy_structure2;
                onerun_parameters(strc_id). batch_ID = save_batch_id;
                onerun_parameters(strc_id). UL_termination = UL_termination_id;
                onerun_parameters(strc_id). gen_gap = gen_gap_list(gg);
                onerun_parameters(strc_id). trg_size = trg_size_list(dd);
                strc_id = strc_id + 1;
            end
        end

    end
end


for dd = 1:length(trg_size_list)
    for  gg = 1:length(gen_gap_list)

        seeds = 1:max_seed;
        ns = length(seeds);
        np = length(prob_list);

        for i = 1 : np
            for j = 1 : ns
                onerun_parameters(strc_id). problem_str = prob_list{i};
                onerun_parameters(strc_id). seed = j;
                onerun_parameters(strc_id). strategy = strategy_structure3;
                onerun_parameters(strc_id). batch_ID = save_batch_id;
                onerun_parameters(strc_id). UL_termination = UL_termination_id;
                onerun_parameters(strc_id). gen_gap = gen_gap_list(gg);
                onerun_parameters(strc_id). trg_size = trg_size_list(dd);
                strc_id = strc_id + 1;
            end
        end

    end
end

num_runs = strc_id-1;

end