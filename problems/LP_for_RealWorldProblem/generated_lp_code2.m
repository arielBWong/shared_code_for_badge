function [solutions, B, solution_one, solution_feasible] = generated_lp_code2(xu)

% Problem Definition
% xu = 0.75;              % Given constant
% Objective coefficients
c1 = [(0.5 + xu); 1];   % f1 = (0.5 + xu)*x1 + 1*x2
c2 = [1; 2];            % f2 = 1*x1 + 2*x2

% Decision variable bounds and equality constraint
lb = [0; 0];
ub = [1; 1];

Aeq = [1, 1];         % x1 + x2 = 1
beq = 1;

theta = atan((xu + 0.5 - 1)/1);  
d3 = (2-1) * sin(theta);
d4 = d3 * sin(theta);
d4 = d4+1;

% fprintf('inequal constraint is %0.2f \n', d4);


A = [-1, -2];
b = [-d4];

B = [];

num_solutions = 20;
options = optimoptions('linprog', 'Display', 'off');

% Step 1: Optimize secondary objective to find lower bound of epsilon
try
    [x_c2, fval_c2, exitflag_c2] = linprog(c2, A, b, Aeq, beq, lb, ub, options);
    if exitflag_c2 ~= 1
        disp('Secondary objective optimization did not succeed.');
        fval_c2 = NaN;
    end
catch ME
    disp('Error in optimizing secondary objective:');
    disp(ME.message);
    fval_c2 = NaN;
end

epsilon_lower = fval_c2; % Lower bound of epsilon

% Step 2: Optimize primary objective to find x*
try
    [x_c1, fval_c1, exitflag_c1] = linprog(c1, A, b, Aeq, beq, lb, ub, options);
    if exitflag_c1 ~= 1
        disp('Primary objective optimization did not succeed.');
        fval_c1 = NaN;
    end
catch ME
    disp('Error in optimizing primary objective:');
    disp(ME.message);
    fval_c1 = NaN;
end

% Step 3: Evaluate secondary objective at x_c1 to find upper bound of epsilon
try
    epsilon_upper = c2' * x_c1; % Upper bound of epsilon
catch ME
    disp('Error in evaluating secondary objective at x_c1:');
    disp(ME.message);
    epsilon_upper = NaN;
    solutions = [];
    return;
end

if isnan(epsilon_lower) || isnan(epsilon_upper)
    disp('Unable to determine epsilon range due to prior errors.');
else
    % Generate epsilon values between lower and upper bounds
    epsilon_values = linspace(epsilon_lower, epsilon_upper, num_solutions);

    % Storage for solutions and objective values
    solutions = [];
    objective1_values = [];
    objective2_values = [];

    % Step 4: Perform epsilon-constraint method
    for epsilon = epsilon_values
        % Add the epsilon constraint to the inequality constraints
        A_eps = [A; c2'];       % Append c2' to the constraint matrix
        b_eps = [b; epsilon];   % Append epsilon to RHS vector

        % Solve the primary objective under epsilon constraint
        try
            [x, fval, exitflag] = linprog(c1, A_eps, b_eps, Aeq, beq, lb, ub, options);
            if exitflag == 1
                solutions = [solutions, x];
                objective1_values = [objective1_values, fval];
                objective2_values = [objective2_values, c2' * x];
            else
                disp(['Epsilon constraint failed at epsilon = ', num2str(epsilon)]);
            end
        catch ME
            disp(['Error in solving epsilon-constraint optimization for epsilon = ', num2str(epsilon)]);
            disp(ME.message);
        end
    end


    for kk = 1:size(solutions, 2)
        solution = solutions(:, kk);
        C = Aeq * solution - beq;
        B = [B; C'];
        if any(C > 1e-3) || any(C < -1e-3)
            disp('Solution is not feasible')
        end
    end

    solutions = solutions';

    % for this problem, pick up one that solution that meet utility
    % function minimal value
    utility_vec = []; % to pick up the minimal utility function value
    for ii = 1:size(solutions, 1)
        f1 = sum(c1' .* solutions(ii, :));
        f2 = sum(c2' .* solutions(ii, :));

        v = (f1 - 1)^2 + (f2 - 1)^2;
        utility_vec = [utility_vec, v];
    end
    [~, idx] = sort(utility_vec);
    solution_one = solutions(idx(1), :);

    f2_separate = sum(c2' .* solution_one);
    f2_all = sum(repmat(c2', size(solutions, 1), 1) .* solutions, 2);
    feasible_idx = find(f2_all >= f2_separate);
    solution_feasible = solutions(feasible_idx, :);
    

%     % Plot Pareto front
%     figure;
%     plot(objective2_values, objective1_values, 'o-');
%     xlabel('Objective 2 (c2^T x)');
%     ylabel('Objective 1 (c1^T x)');
%     title('Pareto Front (Epsilon Constraint Method)');
%     grid on;
% 
%     % Display results
%     disp('Epsilon Constraint Solutions:');
%     disp(solutions);
%     disp('Objective Values:');
%     disp([objective1_values; objective2_values]);
%     close();

end