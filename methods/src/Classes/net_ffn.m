classdef net_ffn < handle
    properties
        Trg_LHS; % LHS training data
        Trg_RHS; % RHS training data
        trg_control; % maximum training data
        net;
    end
    methods
        function obj = net_ffn(hiddenSizes, trg_control)
            obj.Trg_LHS = [];
            obj.Trg_RHS = [];
            obj.net = feedforwardnet(hiddenSizes);

%             obj.net.divideParam.trainRatio = 0.85;  % 70% training
%             obj.net.divideParam.valRatio = 0.15;   % 15% validation
%             obj.net.divideParam.testRatio = 0;  % 15% test
            
            obj.net.trainParam.showWindow = 0;
            obj.trg_control = trg_control;
        end

        function XL = predict(obj, varargin)
            % xu,num_output bu_LHS, bl_LHS, bu_RHS, bl_RHS           

            narginchk(5, 7);
            xu = varargin{1};
            num_output = varargin{2};
           
            if size(obj.Trg_LHS, 1) < 400 % for cases whether network has not been built
                xl_bu = varargin{5};
                xl_bl = varargin{6};
                XL = unifrnd(repmat(xl_bl, num_output, 1), repmat(xl_bu, num_output, 1));
                fprintf("[INFO] Network has not accumulated enough data, initialize with random values \n");

                return;
            end

            if nargin == 7
                bu_LHS = varargin{3};
                bl_LHS = varargin{4};
                bu_RHS = varargin{5};
                bl_RHS = varargin{6};
            else
                bu_LHS = varargin{3};
                bl_LHS = varargin{4};
            end
            
            norm_xu = (xu - bl_LHS ) ./ (bu_LHS - bl_LHS);

            input = [repmat(norm_xu, num_output, 1), linspace(0, 1, num_output)'];
            try
                Norm_XL = obj.net(input');
            catch e
                disp(e);
             
            end
            Norm_XL = Norm_XL';

            if nargin == 7
                XL = bl_RHS + (bu_RHS - bl_RHS) .* Norm_XL;
                maskLB = XL >= bl_RHS;
                maskUB = XL <= bu_RHS;
                XL =  XL .* maskLB .* maskUB + bl_RHS .* (~maskLB) + bu_RHS .* (~maskUB);
            else
                XL = norm_XL;
            end
        end


        function train(obj, bu_LHS, bl_LHS, bu_RHS, bl_RHS, verbose)
            % if size(obj.Trg_LHS, 1) < 400
            if size(obj.Trg_LHS, 1) < 200
                return
            end
            Norm_Trg_LHS = (obj.Trg_LHS - bl_LHS) ./ (bu_LHS - bl_LHS);
            Norm_Trg_RHS = (obj.Trg_RHS - bl_RHS) ./ (bu_RHS - bl_RHS);
            trg_num = size(Norm_Trg_LHS, 1);
            select_idx = trg_num : -1 : (trg_num -obj. trg_control+1);

            if trg_num > obj.trg_control
                Norm_Trg_LHS = Norm_Trg_LHS(select_idx, :);
                Norm_Trg_RHS = Norm_Trg_RHS(select_idx, :);
                obj.Trg_LHS = obj.Trg_LHS(select_idx, :);
                obj.Trg_RHS = obj.Trg_RHS(select_idx, :);
            end

            [obj.net, tr] = train(obj.net,Norm_Trg_LHS',Norm_Trg_RHS');

            if verbose
                fprintf('[INFO] generator NN train data size %d, \n', size(Norm_Trg_LHS, 1));
                fprintf('[INFO] generator NN train time: %0.2f \n', tr.time(end));
            end
        end

        function data_accumulation(obj, LHS, RHS, FLCs)
            % filter LHS and RHS that has inf in it
            [LHS, RHS] = data_filter(obj, LHS, RHS, FLCs);

            if size(LHS, 1) ~= size(RHS, 1)
                error('\n [INFO] training data size inconsistant \n');
            else
                obj.Trg_LHS = [obj.Trg_LHS; LHS];
                obj.Trg_RHS = [obj.Trg_RHS; RHS];
            end

        end

        function [LHS_out, RHS_out] = data_filter(obj, LHS, RHS, FLCs)
            % this function is to filter solutions with Inf objectives out

            % Check if the input matrices have the same number of rows
            if size(LHS, 1) ~= size(RHS, 1) % || size(LHS, 1) ~= size(FUs, 1)
                error('All input matrices must have the same number of rows.');
            end

            % Filter out LL solutions that are infeasible
            if ~isempty(FLCs)
                rowsWithInf = any(FLCs > 1e-6, 2);

                % Remove rows with Inf from LHS and RHS
                LHS_out = LHS(~rowsWithInf, :);
                RHS_out = RHS(~rowsWithInf, :);
            else
                LHS_out = LHS;
                RHS_out = RHS;
            end
        end

    end
end