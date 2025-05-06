classdef real_world3
    properties       
        name;
        ul_bu;
        ul_bl;
        ll_bu;
        ll_bl;
        
        ul_nobj = 2;
        ll_nobj = 2;
        ul_ncon = 1;
        ll_ncon = 1;

    end
    methods
        function obj = real_world3()
            obj.ul_bu = [5];
            obj.ul_bl = [0.5];
            obj.ll_bu = [1, 1];
            obj.ll_bl = [0, 0];
            obj.name = 'real_world3';
            % f = obj.PF_UL(1025);
            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            
            f = [];
            c = [];
        end
        
        function [f, c] = evaluate_l(obj, xu, xl)

            % LL variables   
            q1 = xl(:, 1);
            l1 = xl(:, 2);

            q2 = xl(:, 3);
            l2 = xl(:, 4);

            % UL variables
            w = xu(:, 1);
            t = xu(:, 2);
            r1 = xu(:, 3);
            r2 = xu(:, 4);
            g1 = xu(:, 5);
            g2 = xu(:, 6);


            f1 =  0.7 * q1 - q1 .* t - 0.3 * (45 - q1).^2 + (r1 - q1) .* (0.9 - 0.01 .* (r1 + r2 - q1 - q2)) ...
                - 0.2 * (0.2 * q1 - l1).^2 + (g1 - l1) .* (0.8 - 0.01 * (g1 + g2 - l1 - l2));

            f2 =  0.8 .* q2 - q2 .* t - 0.2 .* (47 - q2).^2 + (r2 - q2) .* (0.9 - 0.01 .* (r1 + r2 - q1 - q2)) ...
                - 0.1 .* (0.3 .* q2 - l2).^2 + (g2 - l2) .* (0.8 - 0.01 .* (g1 + g2 -l1 -l2));

            f = [-f1, -f2];
            
            c1 = c1 + cl2 - 20;
           
            

            c = [];


        end
        
        function pf = PF_LL(obj, n, xu)
           pf = [];
        end

        function ps = PS_LL(obj, n, xu)
           ps = [];
        end
        
        function pf = PF_UL(obj, n)
            PF = [];
  
        end
        
       
        
        function literature_hvmean = return_meanHV(obj)
            literature_hvmean = NaN;
        end
        
        function literature_hvstd = return_stdHV(obj)
            literature_hvstd = NaN;
        end

        function thv = target_hv(obj)
            thv = NaN;
        end

        function tigd = target_igd(obj)
             tigd = NaN;
        end

        function solution = generated_lp_code2(obj, xu)
            % Problem Definition
            % xu = 0.75;              % Given constant
            % Objective coefficients
            c1 = [(0.5 + xu); 1];   % f1 = (0.5 + xu)*x1 + 1*x2
            c2 = [1; 2];            % f2 = 1*x1 + 2*x2

            % Decision variable bounds and equality constraint
            lb = [0; 0];
            ub = [1; 1];

            A = [];
            b = [];

            x1 = linspace(0, 1, 1000);
            x2 = 1 - x1;
            solutions = [x1', x2'];

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
            solution= solutions(idx(1), :);
            end
        end
        
  
end