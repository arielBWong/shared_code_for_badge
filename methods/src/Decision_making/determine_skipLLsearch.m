function skipflag = determine_skipLLsearch(net, gen, gen_gap)
% if net train size is not big enough, no skip
% else gen ~= 5, 10, 15,... can skip
skipflag = false;
if size(net.Trg_LHS, 1) < net.trg_control
    skipflag = false;
else    
    every_period = mod(gen, gen_gap);
    if  every_period ~= 0 % means not in target generation
        skipflag = true;
    end
end
end