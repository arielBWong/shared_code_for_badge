function[XL, FL, FLC, LLcount, varargout] = net_for_periodLL_search(xu, prob, Params, net, tc, skip_flag, varargin)
% skip_flag is for skip lower level search entirely
if skip_flag
    % fprintf('[INFO] evaluation through prediction \n');
    [XL, FL, FLC, LLcount] = netLL_search(xu, prob, Params, net);
    ULcount = 0;
    varargout{1} = ULcount;
else
    % fprintf('[INFO] evaluation through init and search \n');
    if isempty(varargin)
        [XL, FL, FLC, LLcount] = net_as_initNDLL_search(xu, prob, Params, net, tc);
    else
        xl_causal_fl = varargin{1};
        r = varargin{2};
        [XL, FL, FLC, LLcount, ULcount] = net_as_initNDLL_search_causal(xu, prob, Params, net, tc, xl_causal_fl, r);
        varargout{1} = ULcount;
    end
end

end