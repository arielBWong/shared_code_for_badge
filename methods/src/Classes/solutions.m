classdef solutions < handle
    properties
        xu = []; 
        XL = [];
        XU = [];
        XUE = []; % expanded xu, for constructing training data for NN 
        FU = [];
        FL = [];
        FC = [];
        FLC = [];
        addon = [];
%         XUR = [];   % temporary variable for manuscript case study demonstration
%         XURD = [];   % temporary variable for manuscript case study demonstration
%         XUS = [];   % temporary variable for manuscript case study demonstartion 
    end

    methods
        function obj = solutions()
            % this class has dependent on pop_sort method
            obj.XU = [];
            obj.XL = [];
            obj.xu = [];
            obj.XUE = [];
            obj.FU = [];
            obj.FL = [];
            obj.FC = [];
            obj.FLC = [];
            obj.addon = [];
%             obj.XUR = [];  % Randomly generated r, <temporary variable>
%             obj.XURD = []; % Randomly generated r, then sort it,<temporary variable>
%             obj.XUS = [];  % with linspace r but shuffle it.,<temporary variable>
        end

        function add(obj, xu, XL, FU, FL, FC, FLC, varargin)
            % this forms list of cell array
            obj.xu = [obj.xu, {xu}];
            obj.XL = [obj.XL, {XL}];

            if isempty(XL)  % single level problem
                obj.XU = obj.xu;
                obj.XUE = [];
            else
                obj.XU = [obj.XU, {repmat(xu, size(XL, 1), 1)}];
                obj.XUE = [obj.XUE, {[repmat(xu, size(XL, 1), 1), linspace(0, 1, size(XL, 1))']}];
%                 r =  rand(size(XL, 1), 1); %<temporary variable>
%                 obj.XUR = [obj.XUR, {[repmat(xu, size(XL, 1), 1), r]}]; % Bring in random generated r %<temporary variable>
%                 r = sort(r);    %<temporary variable>            
%                 obj.XURD =  [obj.XURD, {[repmat(xu, size(XL, 1), 1), r]}]; % Bring in random generated but sorted r %<temporary variable>
%                 
%                 regular_spaced = linspace(0, 1, size(XL, 1))'; %<temporary variable>
%                 shuffle_id = randperm(size(XL, 1)); %<temporary variable>
%                 regular_butShuffled = regular_spaced(shuffle_id); %<temporary variable>
%                 obj.XUS = [obj.XUS, {[repmat(xu, size(XL, 1), 1), regular_butShuffled]}]; %<temporary variable>

            end
            obj.FU = [obj.FU, {FU}];
            obj.FC = [obj.FC, {FC}];
            obj.FL = [obj.FL, {FL}];
            obj.FLC = [obj.FLC, {FLC}];

            if ~isempty(varargin)
                obj.addon = [obj.addon, {varargin{1}}];
            end
        end

        function merge(obj,  new_solutions)
            obj.xu = [obj.xu, new_solutions.xu];
            obj.XL = [obj.XL, new_solutions.XL];
            obj.XU = [obj.XU, new_solutions.XU];
            obj.FU = [obj.FU, new_solutions.FU];
            obj.FL = [obj.FL, new_solutions.FL];
            obj.FC = [obj.FC, new_solutions.FC];
            obj.FLC = [obj.FLC, new_solutions.FLC];
            obj.XUE = [obj.XUE, new_solutions.XUE];
            obj.addon = [obj.addon, new_solutions.addon];
        end

        function copy(obj, another_solutions)
            if ~isempty(obj.xu)
                obj.XU = [];
                obj.XL = [];
                obj.xu = [];
                obj.XUE = [];
                obj.FU = [];
                obj.FL = [];
                obj.FC = [];
                obj.FLC = [];
            end
            obj.merge(another_solutions);
        end

        function nd_sort(obj, varargin) % only ND is kept in this sorting
            % convert cell array to matrix
            tmpFU = obj.FUs;
            tmpFL = obj.FLs;
            tmpFC = obj.FCs;
            tmpFLC = obj.FLCs;
            tmpXU = obj.XUs;
            tmpXL = obj.XLs;
            tmpAddon = obj.addons;

            
            % ND sort wrapper
            pop.X = tmpXU;
            pop.F = tmpFU;
            pop.C = tmpFC;
            pop.dummy1 = tmpXL;
            pop.dummy2 = tmpFLC;
            pop.dummy3 = tmpFL;
            pop.dummy4 = tmpAddon;

            [pop, front_idx] = pop_sort(pop);
            nd_front_id = front_idx == 1;

            % keep ND
            ndXU = pop.X(nd_front_id, :);
            if isempty(varargin)  % this condition only ND related part (truncate XL, FL, FCL, FU, FC)
                ndFU = pop.F(nd_front_id, :);
                if ~isempty(pop.C)
                    ndFC = pop.C(nd_front_id, :);
                else
                    ndFC = [];
                end

                if ~isempty(pop.dummy1)
                    ndXL = pop.dummy1(nd_front_id, :);
                else
                    ndXL = [];
                end

                if ~isempty(pop.dummy3)
                    ndFL = pop.dummy3(nd_front_id, :);
                else
                    ndFL = [];
                end

                if ~isempty(pop.dummy2)
                    ndFLC = pop.dummy2(nd_front_id, :);
                else
                    ndFLC = [];
                end

                if ~isempty(pop.dummy4)
                    ndAddon = pop.dummy4(nd_front_id, :);
                else
                    ndAddon = [];
                end

                % distribute back to object variables
                unique_ndXU = unique(ndXU, 'rows', 'stable');
                num_ndxu = size(unique_ndXU, 1);
                obj.clear_data;
                for ii = 1:num_ndxu
                    tmp_xu = unique_ndXU(ii, :);
                    ia = ismember(ndXU, tmp_xu, 'rows');
                    %xu, XL, FU, FL, FC, FLC
                    one_xu = tmp_xu;

                    if ~isempty(ndXL)
                        one_XL = ndXL(ia, :);
                    else
                        one_XL = [];
                    end

                    one_FU = ndFU(ia, :);

                    if ~isempty(ndFL)
                        one_FL = ndFL(ia,:);

                    else
                        one_FL = [];
                    end

                    if ~isempty(ndFC)
                        one_FC = ndFC(ia,:);
                    else
                        one_FC = [];
                    end

                    if ~isempty(ndFLC)
                        one_FLC = ndFLC(ia, :);
                    else
                        one_FLC = [];
                    end

                    if ~isempty(ndAddon)
                        one_addon = ndAddon(ia, :);
                    else
                        one_addon = [];
                    end

                    obj.add(one_xu, one_XL, one_FU, one_FL, one_FC, one_FLC, one_addon);
                end
            else % this condition keep all record related to ND xu
                unique_ndXU = unique(ndXU, 'rows', 'stable');
                num_ndxu = size(unique_ndXU, 1);
                obj.clear_data;
                for ii = 1:num_ndxu
                    tmp_xu = unique_ndXU(ii, :);
                    ia = ismember(tmpXU, tmp_xu, 'rows');

                    one_xu = tmp_xu;
                    if ~isempty(tmpXL)
                        one_XL = tmpXL(ia, :);
                    else
                        one_XL = [];
                    end

                    one_FU = tmpFU(ia, :);

                    if ~isempty(tmpFL)
                        one_FL = tmpFL(ia,:);
                    else
                        one_FL = [];
                    end

                    if ~isempty(tmpFC)
                        one_FC = tmpFC(ia,:);
                    else
                        one_FC = [];
                    end

                    if ~isempty(tmpFLC)
                        one_FLC = tmpFLC(ia, :);
                    else
                        one_FLC = [];
                    end
                    obj.add(one_xu, one_XL, one_FU, one_FL, one_FC, one_FLC);
                end
            end
        end

        function DSS_newpopulation(obj, popsize, lb, ub)
            % convert cell array to matrix
            tmpFU = obj.FUs;
            tmpFL = obj.FLs;
            tmpFC = obj.FCs;
            tmpFLC = obj.FLCs;
            tmpXU = obj.XUs;
            tmpXL = obj.XLs;
            tmpAddon = obj.addons;

            % ND sort wrapper
            pop.X = tmpXU;
            pop.F = tmpFU;
            pop.C = tmpFC;
            pop.dummy1 = tmpXL;
            pop.dummy2 = tmpFLC;
            pop.dummy3 = tmpFL;
            pop.dummy4 = tmpAddon;

            [pop, front_idx] = pop_sort(pop);
            nd_idx_binary = front_idx == 1;
            select_candidateID = find(front_idx == 1);
            ndXU = pop.X(nd_idx_binary, :);
            unique_ndXU = unique(ndXU, 'rows', 'stable');
            num_ndxu = size(unique_ndXU, 1);

            if num_ndxu > popsize  % DSS selection
                sparseXID = Sparse_Selection(select_candidateID, pop.X, lb, ub, popsize);
                select_XU_id = select_candidateID(sparseXID);
                % XU selected by select_XU_id should be unique (if not something is wrong)
            else
                tmp_uniqueXU = unique(pop.X, 'rows',  'stable');
                [~, ~, ib] = intersect(tmp_uniqueXU, pop.X, 'stable', 'rows');
                select_XU_id = ib(1:popsize);
            end


            newXU = pop.X(select_XU_id, :);
            tmp_newXU = unique(newXU, "rows", 'stable');
            assert(size(tmp_newXU, 1) == size(newXU, 1), 'DSS selection return repeated solutions after selection');

            obj.clear_data;
            for ii = 1:size(newXU, 1)
                tmp_xu = newXU(ii, :);
                ia = ismember(tmpXU, tmp_xu, 'rows'); % in 2N population, childX and old X have no repeated solution
                one_xu = tmp_xu;
                if ~isempty(tmpXL)
                    one_XL = tmpXL(ia, :);
                else
                    one_XL = [];
                end

                one_FU = tmpFU(ia, :);

                if ~isempty(tmpFL)
                    one_FL = tmpFL(ia,:);
                else
                    one_FL = [];
                end

                if ~isempty(tmpFC)
                    one_FC = tmpFC(ia,:);
                else
                    one_FC = [];
                end

                if ~isempty(tmpFLC)
                    one_FLC = tmpFLC(ia, :);
                else
                    one_FLC = [];
                end

                if ~isempty(tmpAddon)
                    one_addon = tmpAddon(ia, :);
                else
                    one_addon = [];
                end
                obj.add(one_xu, one_XL, one_FU, one_FL, one_FC, one_FLC, one_addon);

            end

        end

        
        function num_nd = nd_sort_no_reduction(obj)
            % This sorting is only for dual population search
            % convert cell array to matrix
            tmpFU = obj.FUs;
            tmpFL = obj.FLs;
            tmpFC = obj.FCs;
            tmpFLC = obj.FLCs;
            tmpXU = obj.XUs;
            tmpXL = obj.XLs;
            tmpAddon = obj.addons;

            % move repeated solutions
            [tmpXU, ia] = unique(tmpXU, 'rows', 'stable');
            tmpFU = tmpFU(ia, :);

            % tmpFL = tmpFL(ia, :);
            if ~isempty(tmpFC)
                tmpFC = tmpFC(ia, :);
            end
            % tmpFLC = tmpFLC(ia, :);
            % tmpXL = tmpXL(ia, :);

            % ND sort wrapper
            pop.X = tmpXU;
            pop.F = tmpFU;
            pop.C = tmpFC;
            pop.dummy1 = tmpXL;
            pop.dummy2 = tmpFLC;
            pop.dummy3 = tmpFL;
            pop.dummy4 = tmpAddon;

            [pop, front_idx] = pop_sort(pop);
            
            % assign pop back to solutions.
            obj.clear_data();            
            obj.add(pop.X, pop.dummy1, pop.F, pop.dummy3, pop.C, pop.dummy2);

            nd_id = front_idx == 1;
            num_nd = sum(nd_id);   
        end


        function FU =  FUs(obj)
            FU = cat(1, obj.FU{:});
        end

        function FL = FLs(obj)
            FL = cat(1, obj.FL{:});
        end

        function FC = FCs(obj)
            FC = cat(1, obj.FC{:});
        end

        function FLC = FLCs(obj)
            FLC = cat(1, obj.FLC{:});
        end

        function XU = XUs(obj)
            XU = cat(1, obj.XU{:});
        end

        function xu = xus(obj)
            xu = cat(1, obj.xu{:});
        end

        function XL = XLs(obj)
            XL = cat(1, obj.XL{:});
        end

        function addons = addons(obj)
            if isempty(obj.addon)
                addons = [];
            else
                addons = cat(1, obj.addon{:});
            end
        end


        function clear_data(obj)
            obj.XU = [];
            obj.XL = [];
            obj.xu = [];
            obj.XUE = [];
            obj.FU = [];
            obj.FL = [];
            obj.FC = [];
            obj.FLC = [];
            obj.addon = [];
        end

    end
end