function[XL, FL, FLC, LLcount] = baselineLL_search(xu, Params, prob, tc)
[XL, FL, FLC, LLcount] = LL_Evolution(xu, Params, prob, [], false, tc);

infeasible_id = any(FLC>1e-6, 2);
XL_feasible = XL(~infeasible_id,:);
FL_feasible = FL(~infeasible_id, :);
FLC_feasible = FLC(~infeasible_id, :);

if isempty(XL_feasible)
    fprintf('[INFO] LL return (NO*) feasible solution \n');
else
    fprintf('[INFO] LL return (FEASIBLE*) solution \n');
end

[XL, FL, FLC] = LL_postprocess(XL, FL, FLC);
end