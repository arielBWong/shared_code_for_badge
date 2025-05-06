function Params = load_parameters(prob, internal_comparison)

if internal_comparison
    Params.N_gen_UL = 10 * length(prob.ul_bu);    % 10 * length(prob.ul_bu);  
    Params.N_gen_LL = 10 * length(prob.ll_bu);
    Params.UL_popsize = 10 * length(prob.ul_bu);  % pop size being 5D
    Params.LL_popsize = 10 * length(prob.ll_bu);
    Params.additional_stopping = false;
else

    Params.UL_popsize = 20;
    Params.LL_popsize = 20;
    Params.N_gen_UL = 1000;
    Params.N_gen_LL = 300; %300;
    Params.additional_stopping = false;
    
end
end
