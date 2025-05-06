function[XL, FL, FLC, LLcount] = netLL_search(xu, prob, Params, net)
%
% input order: xu,num_output bu_LHS, bl_LHS, bu_RHS, bl_RHS, num_ouput
XL = net.predict(xu, Params.LL_popsize, prob.ul_bu, prob.ul_bl, prob.ll_bu, prob.ll_bl);
extend_xu = repmat(xu, size(XL, 1), 1);
[FL, FLC] = prob.evaluate_l(extend_xu, XL);  

visualize = false;
if visualize
    f2 = figure(2);
    pf = prob.PF_LL(100, xu);
    scatter(pf(:, 1), pf(:, 2), 10, 'filled'); hold on;
    scatter(FL(:, 1), FL(:, 2));
    pause(1);
    close(f2);
end
LLcount = size(XL, 1);

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