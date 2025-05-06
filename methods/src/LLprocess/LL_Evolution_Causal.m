function[XL, FL, FLC, LLcount, ULcount, FU0] = LL_Evolution_Causal(xu, parameter, prob, init_pop, vis, termination_condition, xl_causal_fl, r)

LLcount = 0;
ulsearch_type = 3;
FU0 = [];

if r == 1 && ulsearch_type == 1
    related_to_fl = true;  % direct search to LL 
    [XL, LLcount] = partial_search(xu, xl_causal_fl, prob, parameter, related_to_fl, [], init_pop, termination_condition, ulsearch_type);
    ULcount = 0;
    refined_XL = [];
    for jj = 1:size(XL, 1)
        xl = XL(jj, :);
        related_to_fl = false;  % direct search to UL
        termination_condition_UL = 0;
        [XL_update, ul_count] = partial_search(xu, xl_causal_fl, prob, parameter, related_to_fl, xl, [], termination_condition_UL, ulsearch_type);
        ULcount = ULcount + ul_count;
        refined_XL = [refined_XL; XL_update];
    end

    XL = refined_XL;
    [FL, FLC] = prob.evaluate_l(repmat(xu, size(XL, 1), 1), XL);
    LLcount = LLcount + size(XL, 1);

elseif r == 1 && ulsearch_type == 2
    related_to_fl = true;  % direct search to LL
    [XL, LLcount] = partial_search(xu, xl_causal_fl, prob, parameter, related_to_fl, [], init_pop, termination_condition, ulsearch_type);
    ULcount = 0;
    refined_XL = [];

    % for partial UL search
    v = xl_causal_fl == 1;
    related_ids = repmat(v, size(XL, 1), 1);
    XL_LLpartial = XL(related_ids);
    XL_LLpartial = reshape(XL_LLpartial, size(XL, 1), sum(v));

    funh_obj = @(x)ULsearch_objective(xu, x, prob);
    funh_con = @(x)ULsearch_constraint(xu, x, prob);
    v = xl_causal_fl == 0;
    LL_nvar = sum(v);
    sub_ub = prob.ll_bu(v);
    sub_lb = prob.ll_bl(v);
    p.gen = 20;
    p.popsize = 60;

    visualization_ULadd = false;
    pf = [];
    if visualization_ULadd
        LL_ps = prob.PS_LL(129, xu);
        pf = prob.evaluate_u(repmat(xu, 129, 1), LL_ps);
    end

    [XL, ~, ~, archive, ~] = ea_solver(funh_obj, LL_nvar, sub_lb, sub_ub, [], funh_con, p, 'visualize', visualization_ULadd, 'pf', pf, 'termination_criterion', 0,  'init_addon', XL_LLpartial, 'xposition_indicator', v);
    ULcount = size(archive, 1);
    [FL, FLC] = prob.evaluate_l(repmat(xu, size(XL, 1), 1), XL);
    LLcount = LLcount + size(XL, 1);

elseif r == 1 && ulsearch_type == 3
    related_to_fl = true;  % direct search to LL
    [XL, LLcount] = partial_search(xu, xl_causal_fl, prob, parameter, related_to_fl, [], init_pop, termination_condition, ulsearch_type);
    
    replaceXL = unifrnd(repmat(prob.ll_bl,size(XL, 1),1),repmat(prob.ll_bu,size(XL, 1), 1));
    replaceXL = replaceXL(repmat(~xl_causal_fl, size(XL, 1), 1));
    replaceXL = reshape(replaceXL, size(XL, 1), sum(~xl_causal_fl));
    tmp_XL = XL;
    tmp_XL(repmat(~xl_causal_fl, size(XL, 1), 1)) = replaceXL;
    
    [FU0, FC0] = prob.evaluate_u(repmat(xu, size(XL, 1), 1), tmp_XL);
    FC0_id = sum(FC0, 2) <= 0;
    FU0 = FU0(FC0_id, :);

    
    ULcount = 0;
    refined_XL = [];
    refined_FL = [];
    refined_FLC = [];
    for jj = 1:size(XL, 1)
        xl = XL(jj, :);
        related_to_fl = false;  % direct search to UL
        termination_condition_UL = 0;
        [XL_update, ul_count] = UL_partial_search_transfer(xu, xl_causal_fl, prob, parameter, related_to_fl, xl, [], termination_condition_UL, jj, refined_XL, refined_FL, refined_FLC);
        ULcount = ULcount + ul_count;
        refined_XL = [refined_XL; XL_update];

        % prepare reference short cut library
        % assert(size(XL_update, 1) == 1, 'UL search ends up more than one solution ');
        if size(XL_update, 1) > 1
            fprintf('[INFO] UL additional search returned  more than one solution \n');
        end
        [fl, flc] = prob.evaluate_u(repmat(xu, size(XL_update, 1), 1), XL_update); % this is right because borrowing is to borrow the solution with UL feasible objective


        refined_FL = [refined_FL; fl];
        refined_FLC = [refined_FLC; flc];
    end

    XL = refined_XL;
    FL = refined_FL;
    FLC = refined_FLC;
    [FL, FLC] = prob.evaluate_l(repmat(xu, size(XL, 1), 1), XL);
    LLcount = LLcount + size(XL, 1);
else
     [XL, FL, FLC, LLcount] = LL_Evolution(xu, parameter, prob, init_pop, vis, termination_condition);
     ULcount = 0;
end

end


function plotMO(fighn, pop, pf, gen)
clf(fighn);
plot(pf(:, 1), pf(:, 2), 'go'); hold on;
plot(pop(:, 1), pop(:, 2), 'ro');
title(num2str(gen));
grid on;
pause(0.1);
end

%------wrapper----
%----not perfect data structure can cause objective evaluated twice, but it is fine
function  f = objective_func(prob, xu, xl)
xu = repmat(xu, size(xl, 1), 1);
[f, ~] = prob.evaluate_l(xu, xl);
end

function c = constraint_func(prob, xu, xl)
xu = repmat(xu, size(xl, 1), 1);
[~, c] = prob.evaluate_l(xu, xl);
end

function[f] = ULsearch_objective(xu, xl, prob)
xu = repmat(xu, size(xl, 1), 1);
[f, ~] = prob.evaluate_u(xu, xl);
end

function[c] = ULsearch_constraint(xu, xl, prob)
xu = repmat(xu, size(xl, 1), 1);
[~, c] = prob.evaluate_u(xu, xl);
end


function [f] = partial_search_objective_LL(x, xu, v, prob)
xl = zeros(1, length(prob.ll_bu));
[tmp_xu, tmp_xl] = assign_partial_xl(x, xu, xl, v);
[f, ~] = prob.evaluate_l(tmp_xu, tmp_xl);
end

function [c] = partial_search_constraint_LL(x, xu, v, prob)
xl = zeros(1, length(prob.ll_bu));
[tmp_xu, tmp_xl] = assign_partial_xl(x, xu, xl, v);
[~, c] = prob.evaluate_l(tmp_xu, tmp_xl);
end




% Partial search evaluation objective function
function [f] = partial_search_objective_UL(x, xu, xl, v, prob)
% v indicate LL variables relating to UL objective
[tmp_xu, tmp_xl] = assign_partial_xl(x, xu, xl, v);
[f, ~] = prob.evaluate_u(tmp_xu, tmp_xl);
end

function [c] = partial_search_constraint_UL(x, xu, xl, v, prob)
[tmp_xu, tmp_xl] = assign_partial_xl(x, xu, xl, v);
[~, c] = prob.evaluate_u(tmp_xu, tmp_xl);
end


function[tmp_xu, tmp_xl] = assign_partial_xl(x, xu, xl, v)
num_x = size(x, 1);
tmp_v = repmat(v, num_x, 1);

tmp_xu = repmat(xu, num_x, 1);
tmp_xl = repmat(xl, num_x, 1);
tmp_xl(tmp_v) = x;
end

function [XL, count] = UL_partial_search_transfer(xu, xl_causal_fl, prob, parameter, related_to_fl, xl, init_pop, termination_condition, ULid, refined_XL, refined_FL, refined_FLC)
v = xl_causal_fl == related_to_fl;
sub_ub = prob.ll_bu(v);
sub_lb = prob.ll_bl(v);

% p.gen = 40;
% p.popsize = 10; % parameter.LL_popsize;


p.gen = 80;
p.popsize = 5; % parameter.LL_popsize;

funh_obj = @(x)partial_search_objective_UL(x, xu, xl, v, prob);
funh_con = @(x)partial_search_constraint_UL(x, xu, xl, v, prob);

num = 200;
llps = prob.PS_LL(num, xu);
[upf, uc] = prob.evaluate_u(repmat(xu, num, 1), llps);
id = uc <= 0;
pf = upf(id, :);
vis = false;

if ULid > 1
    existingXL = refined_XL;    % existing neighbours
    existingFL = refined_FL;    
    existingFLC = refined_FLC;

    assert(size(existingFLC, 2) == 1, 'following steps only works for one constraint');
    binary_id = existingFLC <= 0;            % only consider feasible neighbours
    if sum(binary_id) > 0
        existingXL = existingXL(binary_id, :);
        existingFL = existingFL(binary_id, :);
        existingFLC = existingFLC(binary_id, :);

        % find out fl related xl variables
        fl_v = xl_causal_fl;
        row_num = size(existingXL, 1);
        col_num = sum(xl_causal_fl);
        
        fl_v = repmat(fl_v, row_num, 1);
        existingXL_short = existingXL(logical(fl_v));
        existingXL_short = reshape(existingXL_short, row_num, col_num);

        new_XL_short = xl(logical(xl_causal_fl));

        % existingXL_short is LL related
        [~, pos_id] = pdist2(existingXL_short, new_XL_short, 'euclidean', 'Smallest', 1); % for each (only one) in new_XL, the distance of existingXL to new_XL (results shown in column vector)
        init_pop = existingXL(pos_id, :);
        % v is picking up UL related
        init_pop = init_pop(logical(v));
        termination_condition = 2;
    end
end

LL_nvar = sum(v);
[XL_partial, ~, ~, archive, ~] = gsolver(funh_obj, LL_nvar, sub_lb, sub_ub, init_pop, funh_con, p, 'visualize', vis, 'pf', pf, 'termination_criterion', termination_condition);

count = size(archive.sols, 1);
if related_to_fl
    tmp_xl = zeros(size(XL_partial, 1), length(prob.ll_bu));
    tmp_v = repmat(v, size(XL_partial, 1), 1);
    tmp_xl(tmp_v) = XL_partial;
    XL = tmp_xl;
else
    additional_num = size(XL_partial, 1);
    tmp_v = repmat(v, additional_num, 1);
    tmp_xl = repmat(xl, additional_num, 1);
    tmp_xl(tmp_v) = XL_partial;
    XL = tmp_xl;
end

end

function [XL, count] = partial_search(xu, xl_causal_fl, prob, parameter, related_to_fl, xl, init_pop, termination_condition, ulsearch_type)
v = xl_causal_fl == related_to_fl;
sub_ub = prob.ll_bu(v);
sub_lb = prob.ll_bl(v);

% p.gen = parameter.LL_gensize;

if related_to_fl

    if ~isempty(init_pop)
        num_cols = sum(v);
        num_rows = size(init_pop, 1);
        ve = repmat(v, num_rows, 1);
        init_pop = init_pop(ve);
        init_pop = reshape(init_pop, num_rows, num_cols);
    end

    p.gen =  parameter.N_gen_LL;
    p.popsize = parameter.LL_popsize;
    funh_obj = @(x)partial_search_objective_LL(x, xu, v, prob);
    funh_con = @(x)partial_search_constraint_LL(x, xu, v, prob);
    pf = prob.PF_LL(129, xu);
    vis = false;

    
else
    p.gen = 40;
    p.popsize = 10; % parameter.LL_popsize;
    funh_obj = @(x)partial_search_objective_UL(x, xu, xl, v, prob);
    funh_con = @(x)partial_search_constraint_UL(x, xu, xl, v, prob);
    
    num = 200;
    llps = prob.PS_LL(num, xu);
    [upf, uc] = prob.evaluate_u(repmat(xu, num, 1), llps);
    id = uc<=0;
    pf = upf(id, :);
    vis = false;
    
    if ulsearch_type == 2


        return
    end
end

LL_nvar = sum(v);
[XL_partial, ~, ~, archive, ~] = gsolver(funh_obj, LL_nvar, sub_lb, sub_ub, init_pop, funh_con, p, 'visualize', vis, 'pf', pf, 'termination_criterion', termination_condition);

count = size(archive.sols, 1);
if related_to_fl
    tmp_xl = zeros(size(XL_partial, 1), length(prob.ll_bu));
    tmp_v = repmat(v, size(XL_partial, 1), 1);
    tmp_xl(tmp_v) = XL_partial;
    XL = tmp_xl;
else
    additional_num = size(XL_partial, 1);
    tmp_v = repmat(v, additional_num, 1);
    tmp_xl = repmat(xl, additional_num, 1);
    tmp_xl(tmp_v) = XL_partial;
    XL = tmp_xl;
end

end 