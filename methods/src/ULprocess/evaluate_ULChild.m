function [XL, FL, FLC, LLcount, varargout] = evaluate_ULChild(strategy, funchn, child_xu, LLnet, gen, gen_gap, varargin)
if strcmp(strategy.name, 'baseline') || strcmp(strategy.name, 'theorectical_LLPS')
    varargout{1} = 0;
    [XL, FL, FLC, LLcount] = funchn(child_xu);
end

if strcmp(strategy.name, 'nets4completeform') || strcmp(strategy.name, 'net4XLNDinit') ...
    [XL, FL, FLC, LLcount] = funchn(child_xu, LLnet);
end

if strcmp(strategy.name, 'set_value') 
   skipflag = true;
   every_period = mod(gen, gen_gap);
   if every_period == 0
       skipflag = false;
       
   end
   [XL, FL, FLC, LLcount] = funchn(child_xu, LLnet, skipflag);
end


if strcmp(strategy.name, 'net2skipLLsearch')    
   skipflag = determine_skipLLsearch(LLnet, gen, gen_gap);  
   [XL, FL, FLC, LLcount] = funchn(child_xu, LLnet, skipflag);
end


if strcmp(strategy.name, 'net2skipLLsearch_causal')
    skipflag = determine_skipLLsearch(LLnet, gen, gen_gap);
    xl_causal_fl = varargin{1};
    r = varargin{2};
    [XL, FL, FLC, LLcount, ULcount] = funchn(child_xu, LLnet, skipflag, xl_causal_fl, r);
    %  '@(x, net,skip_flag, xl_causal_fl, r)net_for_periodLL_search(x, prob, Params, net, 2, xl_causal_fl, r)';
    varargout{1} = ULcount;
end



if strcmp(strategy.name, 'net4XLPS')
    varargout{1} = 0;
    [XL, FL, FLC, LLcount] = funchn(child_xu);
end


if strcmp(strategy.name, 'linprog')
end

end