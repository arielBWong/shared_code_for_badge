function [r, xl_causal_fl] = LL_variable_test(prob)
% causal relation determination
r = false;

xu = unifrnd(prob.ul_bl, prob.ul_bu);
xl = unifrnd(prob.ll_bl, prob.ll_bu);


nl = length(xl);
delta = 0.1;

base_fl = prob.evaluate_l(xu, xl);
xl_causal_fl = ones(1, nl);

% test LL variables on fl
for ii = 1:nl
    delta_xl = xl(ii) + delta;

    if delta_xl >= prob.ll_bu(ii)
        delta_xl = xl(ii) - delta;
    end

    tmp_xl = xl;
    tmp_xl(ii) = delta_xl;

    tmp_fl = prob.evaluate_l(xu, tmp_xl);

    if abs(tmp_fl - base_fl) < 1e-7 % no causal relation
        xl_causal_fl(ii) = 0;
    end
end

if sum(xl_causal_fl) < nl
    r = true;
end

end