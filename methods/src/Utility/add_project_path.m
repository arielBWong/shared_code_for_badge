function add_project_path()
%- --- algorithm path-----
problem_path = fullfile(pwd, 'problems');
addpath(problem_path);
problem_path = fullfile(pwd, 'problems', 'LP_for_RealWorldProblem');
addpath(problem_path);

utility_path = fullfile(pwd, 'methods');
addpath(utility_path);
utility_path = fullfile(pwd, 'methods','src');
addpath(utility_path);
utility_path = fullfile(pwd, 'methods','src', 'ULprocess');
addpath(utility_path);
utility_path = fullfile(pwd, 'methods','src','LLprocess');
addpath(utility_path);
utility_path = fullfile(pwd, 'methods', 'src', 'Utility');
addpath(utility_path);
utility_path = fullfile(pwd, 'methods', 'src', 'Compared_LLprediction');
addpath(utility_path);
utility_path = fullfile(pwd, 'methods', 'src', 'Classes');
addpath(utility_path);
utility_path = fullfile(pwd, 'methods', 'src', 'Decision_making');
addpath(utility_path);
utility_path = fullfile(pwd, 'methods', 'src', 'para_run_utility');
addpath(utility_path);
solver_path = fullfile(pwd, 'methods', 'globalsolver');
addpath(solver_path);
solver_path = fullfile(pwd, 'methods', 'globalsolver','ND_Sort');
addpath(solver_path);
solver_path = fullfile(pwd, 'post_process', 'rawdata_tocollection');
addpath(solver_path);
solver_path = fullfile(pwd, 'post_process', 'compactdata_totable');
addpath(solver_path);




%---- algorithm path ------
end