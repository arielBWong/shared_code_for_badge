function [bestx, bestf, bestc, archive, external_return] = ea_solver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param, varargin)
% This ea solver deal with solutions with addons.
%--------------------------------------------------------------------------
% create argument parser
p = inputParser;
addRequired(p, 'funh_obj');
addRequired(p, 'num_xvar');
addRequired(p, 'lb');
addRequired(p, 'ub');
addRequired(p, 'initmatrix');
addRequired(p, 'funh_con');
addRequired(p, 'param');
addParameter(p, 'externalfunction', []);
addParameter(p, 'visualize', false);
addParameter(p, 'pf', []);
addParameter(p, 'termination_criterion', true);
addParameter(p, 'init_addon', []);
addParameter(p, 'xposition_indicator', []);
parse(p, funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param, varargin{:});
%------- interpret argument ---------------
funh_obj= p.Results.funh_obj;
num_xvar = p.Results.num_xvar;
lb = p.Results.lb;
ub = p.Results.ub;
initmatrix= p.Results.initmatrix;
funh_con = p.Results.funh_con;
param = p.Results.param;
external_funh= p.Results.externalfunction;
visualize = p.Results.visualize;
pf = p.Results.pf;
termination_criterion = p.Results.termination_criterion;
init_addon = p.Results.init_addon;
xposition_indicator = p.Results.xposition_indicator;
%-----------

record = false;
if visualize && record
    f1 = figure(2);
    obj = VideoWriter('moving.avi');
    obj.Quality= 100;
    obj.FrameRate = 25;
    open(obj);
end
if visualize 
    f1 = figure('Position',[100, 100, 600, 600]);
end

% Initialization of solutions
% deal with solutions passed in make sure initmatrix is unique
initmatrix = unique(initmatrix,'rows', 'stable');
n_init = size(initmatrix, 1);
n_rest = param.popsize - n_init;
X_pop =  unifrnd(repmat(lb,n_rest,1),repmat(ub,n_rest,1));
X_pop  = [X_pop; initmatrix];
% make sure initialization unique
X_pop = unique(X_pop, 'rows','stable'); 

archive_x = [];
archive_x = [archive_x; X_pop];

% active population
population = solutions();

% make sure addon has the same size as population size
if size(init_addon, 1) < param.popsize
    n = param.popsize - size(init_addon, 1);
    m = size(init_addon, 1);
    ids = randi(m, 1, n);
    init_addon = [init_addon; init_addon(ids, :)];
end
assert(size(X_pop, 1) == size(init_addon, 1), "solution size is not equal to matching addon size");

% constructing the first population
assemble_X = customized_solution_assemble(X_pop, init_addon, xposition_indicator);
F = funh_obj(assemble_X);
C = funh_con(assemble_X);
population.add(X_pop, [], F, [], C, [], init_addon);

% construct ND archive records
nd_archive = solutions();
nd_archive.copy(population)
nd_archive.nd_sort();

if visualize
%     F = nd_archive.FUs;
%     F = population.FUs;
%     plotMO(f1, F, pf, 0);
    plotMO_seperateFeasible(f1, population.FUs, population.FCs, pf, 1);
end

% start EA search 
for ii = 1:param.gen-1
    % fprintf("generation %d \n", ii+1);
    % generate child population
    child_X = generate_child_DE(lb, ub, population.XUs, param);
    add_on = population.addons;

    % make sure child_X has no repeated solutions, 2N->N rely on child_X
    % and old X has no repeated solutions
    [child_X, ia, ~] = unique(child_X, 'rows', 'stable');
    add_on = add_on(ia, :);
    [child_X, ia] = remove_repeated_solution(archive_x, child_X);
    add_on = add_on(ia, :);
    archive_x = [archive_x; child_X];

    % 2N population
    assemble_X = customized_solution_assemble(child_X, add_on, xposition_indicator);
    F = funh_obj(assemble_X);
    C = funh_con(assemble_X);
    population.add(child_X, [], F, [], C, [], add_on);

    nd_archive.add(child_X, [], F, [], C, [], add_on)
    nd_archive.nd_sort();

    % 2N->N population
    population.DSS_newpopulation(param.popsize, lb, ub);

    if visualize
        %         F = nd_archive.FUs;
        %         F = population.FUs;
        %         plotMO(f1, F, pf, ii+1);
        plotMO_seperateFeasible(f1, population.FUs, population.FCs, pf, ii+1);
    end
end

population.nd_sort();

bestx = population.XUs;
bestf = population.FUs;
bestc = population.FCs;
archive = archive_x;
external_return = [];

bestaddon = population.addons;
bestx =  customized_solution_assemble(bestx, bestaddon, xposition_indicator);

end

function plotMO_seperateFeasible(fighn, f, c, pf, gen)
clf(fighn);
plot(pf(:, 1), pf(:, 2), 'k.'); hold on;
feasible_id = c <= 0;
ff = f(feasible_id, :);
scatter(ff(:, 1), ff(:, 2), 40, 'red', 'filled');
infeasible_id = ~feasible_id;
fc = f(infeasible_id, :);
scatter(fc(:, 1), fc(:, 2), 50, 'blue', 'LineWidth',2);
title(num2str(gen));
pause(0.1);

end

function plotMO(fighn, pop, pf, gen)
clf(fighn);
plot(pf(:, 1), pf(:, 2), 'k.'); hold on;
scatter(pop(:, 1), pop(:, 2),  30, 'r', 'filled');
title(num2str(gen));
grid on;
pause(0.1);
end

function solutions =  customized_solution_assemble(x, addon, xposition_indicator)
assert(size(x, 1) == size(addon, 1), "solution size is not equal to matching addon size");
num_sols = size(x, 1);
num_vals = size(x, 2) + size(addon, 2);

solutions = zeros(num_sols, num_vals);
tmp_xposition_indicator = repmat(xposition_indicator, num_sols, 1);
tmp_addon_indicator = ~tmp_xposition_indicator;

solutions(tmp_xposition_indicator) = x;
solutions(tmp_addon_indicator) = addon;
end