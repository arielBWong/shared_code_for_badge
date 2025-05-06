function[XL, FL, FLC, LLcount, ULcount] = net_as_initNDLL_search_causal(xu, prob, Params, net, tc, xl_causal_fl, r)

% input order: xu,num_output bu_LHS, bl_LHS, bu_RHS, bl_RHS, num_ouput
XL = net.predict(xu, Params.LL_popsize, prob.ul_bu, prob.ul_bl, prob.ll_bu, prob.ll_bl);

% Only give ND front
LLcount = size(XL, 1);

XU = repmat(xu, size(XL, 1), 1);
[initFL, initFC] = prob.evaluate_l(XU, XL);

wrapper.X = XL;
wrapper.F = initFL;
wrapper.C = initFC;
[wrapper, nd_idx] = pop_sort(wrapper);
nd_binary = nd_idx == 1;
XL = wrapper.X(nd_binary, :);
ninit = size(XL, 1);

vis = false;
if vis
    f1 = figure(1);
    XU = repmat(xu, size(XL, 1), 1);
    [initFL, initFC] = prob.evaluate_l(XU, XL);
    pf = prob.PF_LL(1025, xu(1, :));
    plot(pf(:, 1), pf(:, 2), 'k.'); hold on;
    scatter(initFL(:, 1), initFL(:, 2),  30, 'r', 'filled');
    % title(num2str(gen));
    grid on;
    xlabel('f1','FontSize', 22);
    ylabel('f2', 'FontSize', 22);
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',0,'VerticalAlignment','middle');
    % title('LL initialization with PS generator','FontSize', 22);
    legend( 'LL PF', 'LL initialization', 'Location', 'northeast');
    lgd = legend;
    lgd.FontSize = 22;
    lgd.Location = 'northeast';
    pause(0.1);
    close(f1);
end

vis = false;
[XL, FL, FLC, LLcount, ULcount] = LL_Evolution_Causal(xu, Params, prob, XL, false, tc, xl_causal_fl, r);    
[XL, FL, FLC] = LL_postprocess(XL, FL, FLC);     % sorting according to fl

LLcount = LLcount + LLcount - ninit; % there is an re-evaluation step in gsolver

end