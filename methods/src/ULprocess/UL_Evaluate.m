function [FU, FC]=UL_Evaluate(xu,XL, prob)
% function only for MO 

XU = repmat(xu, size(XL, 1), 1);
[FU, FC] = prob.evaluate_u(XU, XL);

if isempty(FC) % assume <=0 is feasible
    FC = zeros(size(XL, 1), 1);
end

end