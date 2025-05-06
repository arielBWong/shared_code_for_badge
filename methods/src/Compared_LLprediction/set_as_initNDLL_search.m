function[XL, FL, FLC, LLcount] = set_as_initNDLL_search(xu, prob, Params, mdls, tc, skipflag)
%-----use set valued prediction to generate XL---
num_xl = size(prob.ll_bu, 2);
size_limit = size(mdls, 1);
XL = [];
for ii = 1:size_limit
    for jj = 1: num_xl
        xl(jj) = mdls{ii, jj}.constant + xu*mdls{ii, jj}.linear + xu*mdls{ii, jj}.sqmatrix*xu';
    end
    XL = [XL; xl];
end

%---------------------
maskLB = XL >= prob.ll_bl;
maskUB = XL <= prob.ll_bu;
XL =  XL .* maskLB .* maskUB + prob.ll_bl .* (~maskLB) + prob.ll_bu .* (~maskUB);

%------------------------------------------------------------------------------
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
FL = wrapper.F(nd_binary, :);
if isempty(initFC)
    FLC = [];
else
    FLC = wrapper.C(nd_binary, :);
end

ninit = size(XL, 1);

if skipflag
    % fprintf("[INFO] LL search skiped use predicted XL \n");
    [XL, FL, FLC] = LL_postprocess(XL, FL, FLC);     % remove repeated XL
    return;
end

% fprintf("[INFO] LL search with EA and use predicted XL \n");

vis = false;

if vis
    XU = repmat(xu, size(XL, 1), 1);
    [initFL, initFC] = prob.evaluate_l(XU, XL);
    f1 = figure(1);
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

[XL, FL, FLC, llnum] = LL_Evolution(xu,Params,prob, XL, vis, tc);
[XL, FL, FLC] = LL_postprocess(XL, FL, FLC);     % remove repeated XL

LLcount = LLcount + llnum - ninit; % there is an re-evaluation step in gsolver

end