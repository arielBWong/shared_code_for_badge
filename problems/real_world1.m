classdef real_world1
    properties       
        name;
        ul_bu;
        ul_bl;
        ll_bu;
        ll_bl;
        
        ul_nobj = 2;
        ll_nobj = 2;
        ul_ncon = 2;
        ll_ncon = 3;

    end
    methods
        function obj = real_world1()
            obj.ul_bu = [1e3, 1e3];
            obj.ul_bl = [0, 0];
            obj.ll_bu = [1e3, 1e3, 1e3];
            obj.ll_bl = [0, 0, 0];
            obj.name = 'real_world1';
            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            f1 = xu(:, 1) + 9 * xu(:, 2) + 10 * xl(:, 1) + xl(:, 2) + 3* xl(:, 3);
            f2 = 9 * xu(:, 1) + 2 * xu(:, 2) + 2 * xl(:, 1) + 7 * xl(:, 2) + 4*xl(:, 3);

            f = [f1, f2];
            f = f .* -1;

            c1 = 3 * xu(:, 1) + 9 * xu(:, 2) + 9 * xl(:, 1) + 5 * xl(:, 2) + 3* xl(:, 3) - 1039;
            c2 = -4 * xu(:, 1) - xu(:, 2) + 3 * xl(:, 1)- 3* xl(:, 2) + 2*xl(:, 3) - 94;
            c = [c1, c2];
        end
        
        function [f, c] = evaluate_l(obj, xu, xl)
            f1 = 4 * xu(:, 1) + 6 * xu(:, 2) + 7 * xl(:, 1) + 4 *xl(:, 2) + 8* xl(:, 3);
            f2 = 6 * xu(:, 1) + 4 * xu(:, 2) + 8 * xl(:, 1) + 7 * xl(:, 2) + 4 * xl(:, 3);

            f = [f1, f2];
            f = f .* -1;

            c1 = 3 * xu(:, 1) - 9 * xu(:, 2) - 9 * xl(:, 1) - 4 *xl(:, 2)  - 61;
            c2 = 5 * xu(:, 1) + 9 * xu(:, 2) + 10 * xl(:, 1) - xl(:, 2) - 2 * xl(:, 3) - 924;
            c3 = 3 * xu(:, 1) - 3 * xu(:, 2) + 1 * xl(:, 2) + 5 *xl(:, 3) - 420;
            c  = [c1, c2, c3];
          
        end
        
        function pf = PF_LL(obj, n, xu)
           pf = [];
        end

        function ps = PS_LL(obj, n, xu)
           ps = [];
        end
        
        function pf = PF_UL(obj, n)
           pf = [-480, -1850];
%            savename = ('real_world1_ULPF1025.mat');
%            save(savename, 'pf');
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
        
        
    end
end