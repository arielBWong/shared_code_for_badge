function [solutions, B] = generated_lp_code(xu)

% Define objective coefficients for the two objectives
c1 = [-7; -4; -8];  % Coefficients for the primary objective
c2 = [-8; -7; -4];  % Coefficients for the secondary objective

% Define the constraint matrices
A = [-9, -4, 0; 10, -1, -2; 0, 1, 5];                                           % Coefficient matrix for inequalities
b = [61 - 3*xu(1) + 9*xu(2); 924 - 5*xu(1) - 9*xu(2); 420 - 3*xu(1) + 3*xu(2)]; % Right-hand side of the inequalities
B = [];

options = optimoptions('linprog', 'Display', 'none'); % Suppress display

% Step 1: Optimize secondary objective to find lower bound of epsilon
try
    [x_c2, fval_c2, exitflag_c2] = linprog(c2, A, b, [], [], zeros(size(c1)), [], options);
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
    [x_c1, fval_c1, exitflag_c1] = linprog(c1, A, b, [], [], zeros(size(c1)), [], options);
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
    epsilon_values = linspace(epsilon_lower, epsilon_upper, 20);

    % Storage for solutions and objective values
    solutions = [];
    objective1_values = [];
    objective2_values = [];

    % Step 4: Perform epsilon-constraint method
    for epsilon = epsilon_values
        % Add the epsilon constraint to the inequality constraints
        A_eps = [A; c2']; % Append c2' to the constraint matrix
        b_eps = [b; epsilon]; % Append epsilon to RHS vector

        % Solve the primary objective under epsilon constraint
        try
            [x, fval, exitflag] = linprog(c1, A_eps, b_eps, [], [], zeros(size(c1)), [], options);
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
        C = A * solution - b;
        B = [B; C'];
        if any(C > 1e-6)
            disp('Solution is not feasible')
        end

    end
    
    solutions = solutions';
    % % Plot Pareto front
    % figure;
    % plot(objective2_values, objective1_values, 'o-');
    % xlabel('Objective 2 (c2^T x)');
    % ylabel('Objective 1 (c1^T x)');
    % title('Pareto Front (Epsilon Constraint Method)');
    % grid on;
    %
    % % Display results
    % disp('Epsilon Constraint Solutions:');
    % disp(solutions);
    % disp('Objective Values:');
    % disp([objective1_values; objective2_values]);


end