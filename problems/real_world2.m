classdef real_world2
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
        function obj = real_world2()
            obj.ul_bu = [5];
            obj.ul_bl = [0.5];
            obj.ll_bu = [1, 1];
            obj.ll_bl = [0, 0];
            obj.name = 'realworld2';
            % f = obj.PF_UL(1025);
            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            f1 = -xu .* xl(:, 1);
            f2 = xl(:, 1) + 2 * xl(:, 2);
            f = [f1, f2];
            c = [];
        end
        
        function [f, c] = evaluate_l(obj, xu, xl)
            f1 = (0.5 + xu) .* xl(:, 1) + 1 .* xl(:, 2);
            f2 = xl(:, 1) + 2 .* xl(:, 2);
            f = [f1, f2];
            c1 = xl(:, 1) + xl(:, 2) - 1 - 1e-3;
            c2 = 1 - 1e-3 - xl(:, 1) - xl(:, 2);

            theta = atan((xu(:, 1) + 0.5 - 1)/1);
            d3 = (2-1) .* sin(theta);
            d4 = d3 .* sin(theta);
            d4 = d4+1;

            c3 = d4 - xl(:, 1) - 2 .* xl(:, 2);

            c = [c1, c2, c3];


        end
        
        function pf = PF_LL(obj, n, xu)
           pf = [];
        end

        function ps = PS_LL(obj, n, xu)
           ps = [];
        end
        
        function pf = PF_UL(obj, n)
            xu = linspace(0.5, sqrt(5)/2, 2*n);
            pf = [];
 
            for ii = 1:2*n
                xl = obj.generated_lp_code2(xu(ii));
                fu = obj.evaluate_u(xu(ii), xl);
                pf = [pf; fu];
            end

            pf = sparse_pf_selection(pf, n);
            % plot(pf(:, 1), pf(:, 2), 'o-'); 
            
            % savename = 'real_world2_ULPF1025.mat';
            % save(savename, "pf");
  

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