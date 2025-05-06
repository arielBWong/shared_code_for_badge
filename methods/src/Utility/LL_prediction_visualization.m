function LL_prediction_visualization(gen, prob, net,Params)
% Visualization for reviewer 1 comment 3

if gen == 0
    nu = 40;
    filename = sprintf('%s_test_data.mat', prob.name);
    filename = fullfile(pwd, 'post_process', 'prediction_test', filename);
    load(filename);        % load xu
    my_igd = ones(1, 21) .* -1;
    c_igd = ones(1, 21) .* -1;

    for ii = 1:21
        filename = sprintf('%s_%d_prediction_result_no_%d.mat', prob.name, nu, ii);
        filename = fullfile(pwd, 'post_process', 'prediction_test/cG-BLEMO', filename);
        load(filename);    % get xu, fl, flc, xl from cGAN record

        % load net
        filename = sprintf('%s_psp_net_%d.mat', prob.name, nu);
        filename = fullfile(pwd, 'post_process','prediction_test', filename);
        load(filename);

        % load xu
        fl_cgan = fl;
        xu = TESTXU(ii, :);
        [my_igd(ii), c_igd(ii)] = visualization(xu, fl_cgan, flc, net, prob, 0, Params, ii+nu);
    end

    filename = sprintf('%s_test_igd_%d.mat', prob.name, nu);
    filename = fullfile(pwd, 'post_process','prediction_test', filename);
    save(filename, "my_igd", "c_igd");

    [sorted, id] = sort(my_igd);
    fprintf("[INFO] median IGD PSP %0.4f, id %d \n", sorted(11), id(11));

    [sorted, id] = sort(c_igd);
    fprintf("[INFO] median IGD cG %0.4f, id %d \n", sorted(11), id(11));

    return 
end


if gen < 2
    filename = sprintf('%s_gen1_prediction_record.mat', prob.name);
    filename = fullfile(pwd, 'post_process', 'LL_predication_visualization', filename);
    load(filename); % get xu, fl, flc, xl from cGAN record
    % load cgan fl predication result.
    fl_cgan = fl;
    visualization(xu, fl_cgan, flc, net, prob, 1, Params);

else
    filename = sprintf('%s_gen8_prediction_record.mat', prob.name);
    filename = fullfile(pwd, 'post_process', 'LL_predication_visualization', filename);
    load(filename); % get xu, fl, flc, xl from cGAN record
    % load cgan fl predication result.
    fl_cgan = fl;
    visualization(xu, fl_cgan, flc, net, prob, 2, Params);

end

end




function [my_igd, c_igd] = visualization(xu, fl_cgan, flc_cgan, net, prob, gen, Params, id)

pf = prob.PF_LL(1025, xu);
XL = net.predict(xu, Params.LL_popsize, prob.ul_bu, prob.ul_bl, prob.ll_bu, prob.ll_bl);
[FL, FLC] = prob.evaluate_l(repmat(xu, 20, 1), XL);

%----------calculate igd------------------------------
ll_ideal = min(pf, [], 1);
ll_nadir = max(pf, [], 1);
norm_FL = (FL - ll_ideal) ./ (ll_nadir - ll_ideal);
norm_pf = (pf - ll_ideal) ./ (ll_nadir - ll_ideal);
my_igd = mean(min(pdist2(norm_pf, norm_FL), [], 2));

norm_cFL = (fl_cgan - ll_ideal) ./ (ll_nadir - ll_ideal);
c_igd = mean(min(pdist2(norm_pf, norm_cFL), [], 2));
%--------------------------------------------------------


f1 = figure(1);

scatter(pf(:, 1), pf(:, 2), 20, [0 0.4470 0.7410],'filled', 'DisplayName', 'LL PF');
hold on;
scatter(FL(:, 1), FL(:, 2), 80,[0.8500 0.3250 0.0980],  'filled', 'DisplayName', 'PSP-BLEMO');
scatter(fl_cgan(:, 1), fl_cgan(:, 2), 80,[0.9290 0.6940 0.1250], 'filled', 'DisplayName','cG-BLEMO');
% ax = gca;
% ax.XAxis.FontSize = 12;
% ax.YAxis.FontSize = 12;
legend('FontSize', 20);
xlabel('f^{1}','FontSize',24);
ylabel('f^{2}', 'FontSize',24, 'Rotation',0);
box on;

%-----inset 
inset = false;
if inset
    ax = axes('Position', [0.2 0.6 0.25 0.25]); % [left bottom width height]

    % inset_fl_cgan = fl_cgan(fl_cgan(:, 1) < 15 & fl_cgan(:, 2)< 15, :); % inset_fl_cgan = fl_cgan(fl_cgan(:, 1) < 5 & fl_cgan(:, 2)< 10, :); % 5/10 for DS3
    inset_fl_cgan = fl_cgan(fl_cgan(:, 1) < 5 & fl_cgan(:, 2)< 10, :); % 5/10 for DS3
    scatter(pf(:, 1), pf(:, 2), 20, [0 0.4470 0.7410], 'filled', 'DisplayName', 'LL PS'); hold on;
    scatter(inset_fl_cgan(:, 1), inset_fl_cgan(:, 2), 20,[0.9290 0.6940 0.1250], 'filled', 'DisplayName','cG-BLEMO');

    % inset_FL = FL(FL(:, 1)< 15 & FL(:, 2)< 15, :);  %inset_FL = FL(FL(:, 1)< 5 & FL(:, 2)< 10, :); % 5/10 for DS3
    inset_FL = FL(FL(:, 1)< 5 & FL(:, 2)< 10, :); % 5/10 for DS3
    scatter(inset_FL(:, 1), inset_FL(:, 2), 20,[0.8500 0.3250 0.0980], 'filled', 'DisplayName', 'PSP-BLEMO');
    xlabel('f^{1}','FontSize',14);
    ylabel('f^{2}', 'FontSize',14, 'Rotation',0);
    box on; % Draw box around the inset
end
% savename = sprintf('%s_test_gen%d.fig', prob.name, gen);
% savename = fullfile(pwd, 'post_process', 'LL_predication_visualization', savename);
% saveas(f1, savename);
% 
% savename = sprintf('%s_test_gen%d.png', prob.name, gen);
% savename = fullfile(pwd,  'post_process', 'LL_predication_visualization', savename);
% saveas(f1, savename);


savename = sprintf('%s_test_gen%d_%d.fig', prob.name, gen, id);
savename = fullfile(pwd, 'post_process', 'prediction_test/visualization/800', savename);
saveas(f1, savename);

savename = sprintf('%s_test_gen%d_%d.png', prob.name, gen, id);
savename = fullfile(pwd,  'post_process', 'prediction_test/visualization/800', savename);
saveas(f1, savename);

close(f1);
end