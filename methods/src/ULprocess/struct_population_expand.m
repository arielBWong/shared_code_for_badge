function pop = struct_population_expand(XU, FU, FC, FLC, XL, FL, pop, varargin)
% this expansion assume pop has structure of xu, fu, fc, xl, fl, flc
pop.XU = [pop.XU; XU];
pop.FU = [pop.FU; FU];
pop.FC = [pop.FC; FC];
pop.XL = [pop.XL; XL];
pop.FL = [pop.FL; FL];
pop.FLC = [pop.FLC; FLC];
if ~isempty(varargin)
    pop.LLcount = [pop.LLcount; varargin{1}];
end
end