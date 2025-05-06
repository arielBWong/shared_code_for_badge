function [XL, FL, FLC, LLcount, feasible_exist, xlone, xfeasible] = LL_LinProg(xu, prob)

[xl, B, xlone, xfeasible] = generated_lp_code2(xu);

if ~isempty(xl)
    XL = unique(xl, "rows", "stable");
    nxl = size(XL, 1);
    [FL, FLC] = prob.evaluate_l(repmat(xu, nxl, 1), XL);
    feasible_exist = true;

    % visualize_LLfront(XL, FL, FLC);
else
    FL = ones(1, prob.ll_nobj) .* inf;
    FLC = ones(1, prob.ll_ncon) .* inf;
    XL = ones(1, length(prob.ll_bu));
    feasible_exist = false;
end

LLcount = 0;


end