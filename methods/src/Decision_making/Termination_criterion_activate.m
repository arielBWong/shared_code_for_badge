function [termination_flag, termination_vector] = Termination_criterion_activate(NDarchive, termination_criterion,  gen, termination_vector, verbose)

if termination_criterion == 1       % IGD termination
    % archive, sliding_window, threshold, gen, termination_vector, verbose
    [termination_flag, termination_vector] = Termination_criterion_IGD(NDarchive, 5, 0.01, gen, termination_vector, verbose);
    return

end

if termination_criterion == 2       % HV termination
    % archive, sliding_window, threshold, gen, termination_vector, verbose
    [termination_flag, termination_vector] = Termination_criterion_HV(NDarchive, 10, 0.001, gen, termination_vector, verbose);
    return
end


termination_flag = false;


end