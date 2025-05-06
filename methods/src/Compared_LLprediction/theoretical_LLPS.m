function[XL, FL, FLC] = theoretical_LLPS(xu, prob, Params)
n = Params.LL_popsize;
XL = prob.PS_LL(n, xu);

extend_xu = repmat(xu, size(XL, 1), 1);
[FL, FLC] = prob.evaluate_l(extend_xu, XL);  

global ll_num
ll_num = ll_num + size(XL, 1);

% Only return ND front
popwrapper.X = XL;
popwrapper.F = FL;
popwrapper.C = FLC;
[popwrapper, front_idx] = pop_sort(popwrapper);
select_binary = front_idx == 1;
XL = popwrapper.X(select_binary, :);
FL = popwrapper.F(select_binary, :);
if ~isempty(FLC)
    FLC = popwrapper.C(select_binary, :);
end

[XL, FL, FLC] = LL_postprocess(XL,FL, FLC); % remove repeated XL
end