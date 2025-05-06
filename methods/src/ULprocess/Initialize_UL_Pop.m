function XU = Initialize_UL_Pop(prob,Params)
% Create an initial population using LHS of XU for the problem
ub = prob.ul_bu;
lb = prob.ul_bl;



% XU = lhsdesign(Params.UL_popsize,size(ub, 2), 'criterion','maximin');
% XU = lb + (ub - lb) .* XU;

XU = unifrnd(repmat(lb,Params.UL_popsize,1),repmat(ub,Params.UL_popsize,1));

end
