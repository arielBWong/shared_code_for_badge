function remove_project_path()
%- --- algorithm path-----
problem_path = fullfile(pwd, 'problems');
rmpath(problem_path);
utility_path = fullfile(pwd, 'methods');
rmpath(utility_path);
utility_path = fullfile(pwd, 'methods','src');
rmpath(utility_path);
utility_path = fullfile(pwd, 'methods','src', 'ULprocess');
rmpath(utility_path);
utility_path = fullfile(pwd, 'methods','src','LLprocess');
rmpath(utility_path);
utility_path = fullfile(pwd, 'methods', 'src', 'Utility');
rmpath(utility_path);
utility_path = fullfile(pwd, 'methods', 'src', 'Compared_LLprediction');
rmpath(utility_path);
utility_path = fullfile(pwd, 'methods', 'src', 'Classes');
rmpath(utility_path);
utility_path = fullfile(pwd, 'methods', 'src', 'Decision_making');
rmpath(utility_path);
utility_path = fullfile(pwd, 'methods', 'src', 'para_run_utility');
rmpath(utility_path);
solver_path = fullfile(pwd, 'methods', 'globalsolver');
rmpath(solver_path);
solver_path = fullfile(pwd, 'methods', 'globalsolver','ND_Sort');
rmpath(solver_path);
solver_path = fullfile(pwd, 'post_process', 'rawdata_tocollection');
rmpath(solver_path);
end