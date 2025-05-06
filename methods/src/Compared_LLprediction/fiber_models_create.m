function model_matrix = fiber_models_create(Active_pop)
% this function is to create fiber models of each XL variables
% determine the number of models
num_xlvar = size(Active_pop.XL, 2);
Active_pop = squeeze_pop(Active_pop);

[~, ia, ~] = unique(Active_pop.XU, 'rows', 'stable');
XU = Active_pop.XU(ia, :);

% ----- remove repeated solutions
% Determine size of training data extracted from each XL collection
nd_size_list = [];
for jj = 1: size(XU, 1)
    ia = ismember(Active_pop.XU, XU(jj, :), 'rows');
    nd_size_list = [nd_size_list, sum(ia)];
end
size_limit = min(nd_size_list);

%----------------------
% sortXL in accending order for each variable
sorted_XL = cell(1, size(XU, 1));
for ii = 1: size(XU, 1)
    XL = extract_XL(Active_pop, XU(ii, :), size_limit);
    [~, ix] = sort(XL(:, 1)) ;
    XL = XL(ix, :);                             % sort according to the first variable
    sorted_XL{ii} = XL;
end

% create models for  each variable
model_matrix = cell(size_limit, num_xlvar);
for xl_fiber = 1: size_limit  
    for xl_ii = 1:num_xlvar
        RHS = [];
        for xl_jj = 1:size(XU, 1)
            XL = sorted_XL{xl_jj};
            RHS = [RHS; XL(xl_fiber, xl_ii)];            
        end
        LHS = XU;
        model_matrix{xl_fiber, xl_ii} = quadApprox(RHS, LHS);
    end
end
end


function XL = extract_XL(Active_pop, xu, num)
ia = ismember(Active_pop.XU, xu, 'rows');
XL = Active_pop.XL(ia, :);
XL = XL(1:num, :); 
end


function Active_compactpop = squeeze_pop(Active_pop)
[~, ia, ~] = unique([Active_pop.XU, Active_pop.FU, Active_pop.FC], 'rows', 'stable');
Active_compactpop.XU = Active_pop.XU(ia, :);
Active_compactpop.FU = Active_pop.FU(ia, :);
Active_compactpop.FC = Active_pop.FC(ia, :);
Active_compactpop.XL = Active_pop.XL(ia, :);
Active_compactpop.FL = Active_pop.FL(ia, :);
Active_compactpop.FLC = Active_pop.FLC(ia, :);
end
