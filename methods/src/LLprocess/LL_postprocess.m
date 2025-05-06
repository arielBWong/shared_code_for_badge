function [XL, FL, FLC] = LL_postprocess(XL, FL, FLC)
% This function is used for indexing Uppper level sorting
% first unique solutions are kepted
% and then FLC is converted into one value (indicate cc or cv)
% 0 means cc (constraint comply), > 0 means constraint violation (cv)
% 
[XL,ind] = unique(XL, 'rows', 'stable');
FL = FL(ind, :);
if ~isempty(FLC)
    FLC = FLC(ind, :);
else
    FLC = zeros(size(XL, 1), 1); % convert all feasible
end

%This feasibiliy sum process, considered possible problem with negative value being
%too big so that after sum, FLC shows feasible, but infact infeasible
feasible_binary = FLC <= 0;
FLC(feasible_binary) = 0;
FLC = sum(FLC, 2);


% sort XL according to FL for network training
[~, idx] = sort(FL(:, 1));
XL = XL(idx, :);
FL = FL(idx, :);
FLC = FLC(idx, :);

end