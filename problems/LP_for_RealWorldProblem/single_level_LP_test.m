function solutions =  single_level_LP_test()
% this method is used to test a toy problem for LP solver
addpath(fullfile(pwd, 'methods', 'src','Utility'));
add_project_path()
prob = real_world1();
xu = [25, 65];


A = [-9, -4, 0; 10, -1, -2; 0, 1, 5];
b = [61 - 3*xu(1) + 9*xu(2); 924 - 5*xu(1) - 9*xu(2); 420 - 3*xu(1) + 3*xu(2)];

n = 20;

nd = [];
ndc = [];
solutions = [];



options = optimoptions('linprog','Display','none');

% find minimum value of object 1
f1_coefficient_vector = [-7, -4, -8];
[xf1, f1min] = linprog(f1_coefficient_vector, A, b, [], [], zeros(1,3), ones(1, 3)*1000, options);


if isempty(xf1)
    fprintf("[INFO] LL search find no solution \n")
    return;
end

% find minimum value of object 2
f2_coefficient_vector = [-8, -7, -4];
[xf2, f2min] = linprog(f2_coefficient_vector, A, b, [], [], zeros(1, 3), ones(1, 3)*1000, options);

extreme_f2 = sum(f2_coefficient_vector .* xf1'); % f2 equation, but f1 solution, gives max f2 on PF

if abs(extreme_f2 - f2min)< 1e-6
    fprintf("[INFO] This MO is not conflicting \n");
    solutions = xf2';
    return
end

% re-minimize f1 but with f2 turned into contraint.
for ii = 1:n
    additional_constraintb = f2min + (extreme_f2 - f2min) * ii/n;

    Anew = [A; f2_coefficient_vector];
    bnew = [b; additional_constraintb];

    [x, fval] = linprog(f1_coefficient_vector, Anew,  bnew, [], [], zeros(1, 3), ones(1, 3)*10, options);
   
    if isempty(x)
        continue;
    end

    solutions = [solutions; x'];
    

    [f,c] = prob.evaluate_l(xu, x');
    nd = [nd; f];
    ndc = [ndc; c]; 
end


plot(nd(:, 1), nd(:, 2), 'ro');
remove_project_path()
end

