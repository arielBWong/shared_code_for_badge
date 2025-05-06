function flag = update_net_decision(strategy, gen, net, varargin)
flag = false;


if strcmp(strategy, 'net2skipLLsearch')
    gen_gap = varargin{1};
    skipflag = determine_skipLLsearch(net, gen, gen_gap);
    flag = ~skipflag;
elseif strcmp(strategy, 'net2skipLLsearch_causal')
    gen_gap = varargin{1};
    skipflag = determine_skipLLsearch(net, gen, gen_gap);
    flag = ~skipflag;

elseif strcmp(strategy, 'net4XLNDinit')
    flag = true;
elseif strcmp(strategy, 'set_value') % set value skip this checking
    flag = false;
end
end