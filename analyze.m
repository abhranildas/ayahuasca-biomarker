%% 1. Parse data
fname = 'biomarkers.xlsx';
full_sheet = readtable(fname);

% Drop markers we are not analysing
full_sheet(:, {'ast','rbc','Hct','leu','seg_abs','seg_rel','mon_rel', ...
    'eos_rel','neu_rel','lym_rel','totalchol'}) = [];

% Drop implausible values (likely measurement error)
full_sheet.corti_sal(full_sheet.corti_sal>1500) = nan;
full_sheet.alt(full_sheet.alt>100) = nan;

% Log-transform skewed markers toward normality. The fuller candidate list is
% kept (commented) for reference; presently only corti_sal is transformed.
% log_markers = {'alt','ast_alt','neu_abs','eos_abs','lym_abs','mon_abs','plate','SII',...
% 'glucose','triglycerides','corti_plasm','il6','bdnf','corti_sal'};
log_markers = {'corti_sal'};

% Replace non-positive values with each column's smallest positive value (so
% the log is defined), then log-transform in place.
temp_mat = full_sheet{:, log_markers};
temp_mat_nan = temp_mat;
temp_mat_nan(temp_mat_nan <= 0) = nan;
min_vals = min(temp_mat_nan, [], 1, 'omitnan');
for c = 1:size(temp_mat, 2)
    temp_mat(temp_mat(:, c) <= 0, c) = min_vals(c);
end
full_sheet{:, log_markers} = log(temp_mat);

% Marker names (exclude MADRS, which is the clinical score, not a biomarker)
all_num_cols = full_sheet(:, vartype('numeric')).Properties.VariableNames;
markers = all_num_cols(~strcmpi(all_num_cols, 'MADRS'));
n_markers = numel(markers);

% Split by timepoint, then by group/treatment. The *_full tables retain MADRS.
before = full_sheet(strcmpi(full_sheet.timepoint, 'before'), :);
after = full_sheet(strcmpi(full_sheet.timepoint, 'after'), :);

% Baseline: H0 (healthy), D0 (depressed)
H0_full = before(strcmpi(before.group, 'H'), :);
D0_full = before(strcmpi(before.group, 'D'), :);

% Post-treatment: Ha/Da (Ayahuasca), Hp/Dp (Placebo)
Ha_full = after(strcmpi(after.group, 'H') & strcmpi(after.treatment, 'Ayahuasca'), :);
Hp_full = after(strcmpi(after.group, 'H') & strcmpi(after.treatment, 'Placebo'), :);
Da_full = after(strcmpi(after.group, 'D') & strcmpi(after.treatment, 'Ayahuasca'), :);
Dp_full = after(strcmpi(after.group, 'D') & strcmpi(after.treatment, 'Placebo'), :);

% Restrict to the marker columns only
H0 = H0_full(:, markers);
D0 = D0_full(:, markers);
Ha = Ha_full(:, markers);
Hp = Hp_full(:, markers);
Da = Da_full(:, markers);
Dp = Dp_full(:, markers);

% Z-score every group against the healthy-baseline mean and SD
[sd_H0,mu_H0] = std(H0,0,'omitnan');
H0_z=(H0-mu_H0)./sd_H0;
D0_z=(D0-mu_H0)./sd_H0;
Ha_z=(Ha-mu_H0)./sd_H0;
Hp_z=(Hp-mu_H0)./sd_H0;
Da_z=(Da-mu_H0)./sd_H0;
Dp_z=(Dp-mu_H0)./sd_H0;

%% 2. Greedily accumulate markers for H0/D0 separation
% Run first so the individual-marker figure below can use the greedy top-2.
num_runs = 100;
K_folds = 5;
% all_ranks(run, marker): the sequence position each marker took in that run
all_ranks = zeros(num_runs, n_markers);
for i_run = 1:num_runs
    fprintf('Running greedy accumulation %d / %d...\n', i_run, num_runs);

    % Start each run with every marker available
    rest_markers = markers;
    rest_H0_z = H0_z{:,:};
    rest_D0_z = D0_z{:,:};

    greedy_markers_run = cell(1, n_markers);
    greedy_H0_z = [];
    greedy_D0_z = [];

    % Add markers one at a time, each step picking the one that minimizes
    % single-shot 5-fold CV error when appended to the current panel.
    for i_greedy = 1:n_markers
        num_markers_rest = size(rest_H0_z, 2);
        err_check = nan(1, num_markers_rest);

        for i_check = 1:num_markers_rest
            test_H = [greedy_H0_z, rest_H0_z(:, i_check)];
            test_D = [greedy_D0_z, rest_D0_z(:, i_check)];

            [samp_err, ~] = cv_classify_error(test_H, test_D, K_folds, 1);
            err_check(i_check) = samp_err;
        end

        [~, idx_best] = min(err_check);

        % Append the winner, then remove it from the pool
        greedy_markers_run{i_greedy} = rest_markers{idx_best};
        greedy_H0_z = [greedy_H0_z, rest_H0_z(:, idx_best)];
        greedy_D0_z = [greedy_D0_z, rest_D0_z(:, idx_best)];
        rest_markers(idx_best) = [];
        rest_H0_z(:, idx_best) = [];
        rest_D0_z(:, idx_best) = [];
    end

    % Record where each marker landed in this run's sequence
    [~, ranks_this_run] = ismember(markers, greedy_markers_run);
    all_ranks(i_run, :) = ranks_this_run;
end

% Stable ordering: markers sorted by mean greedy rank (best first)
mean_ranks = mean(all_ranks, 1);
rank_ci = prctile(all_ranks, [12.5 87.5], 1);  % 75% interval, empirical (ranks are bounded, not normal)
[sorted_mean_ranks, sort_idx] = sort(mean_ranks);
sorted_rank_ci_lo = sorted_mean_ranks - rank_ci(1, sort_idx);
sorted_rank_ci_hi = rank_ci(2, sort_idx) - sorted_mean_ranks;
stable_markers = markers(sort_idx);
top2 = stable_markers(1:2);

% Cumulative performance of the stable panel (markers 1..i at each step)
[~, stable_idx] = ismember(stable_markers, markers);
stable_H0_z = H0_z{:, stable_idx};
stable_D0_z = D0_z{:, stable_idx};
cum_err_mean = nan(1, n_markers);
cum_acc_ci_lo = nan(1, n_markers);
cum_acc_ci_hi = nan(1, n_markers);
cum_linear_err_mean = nan(1, n_markers);
cum_linear_acc_ci_lo = nan(1, n_markers);
cum_linear_acc_ci_hi = nan(1, n_markers);
n_perms = 1000;
cum_null_floor = nan(1, n_markers);
cum_null_upper_999 = nan(1, n_markers);
cum_null_upper_99 = nan(1, n_markers);
cum_null_upper_95 = nan(1, n_markers);
for i = 1:n_markers
    % Accumulate markers 1 through i
    eval_H = stable_H0_z(:, 1:i);
    eval_D = stable_D0_z(:, 1:i);

    % Evaluate this fixed panel with multi-rep 5-fold CV (QDA, then linear)
    [cum_err_mean(i), ~, rep_errs] = cv_classify_error(eval_H, eval_D, K_folds, num_runs);
    acc_ci = prctile(1 - rep_errs, [12.5 87.5]);
    cum_acc_ci_lo(i) = acc_ci(1); cum_acc_ci_hi(i) = acc_ci(2);
    [cum_linear_err_mean(i), ~, rep_errs_linear] = cv_classify_error(eval_H, eval_D, K_folds, num_runs, 'linear', true);
    linear_acc_ci = prctile(1 - rep_errs_linear, [12.5 87.5]);
    cum_linear_acc_ci_lo(i) = linear_acc_ci(1); cum_linear_acc_ci_hi(i) = linear_acc_ci(2);

    % Permutation test for this cumulative step -> 1-tailed null CIs
    fprintf('Evaluating cumulative step %d/%d and running %d permutations...\n', i, n_markers, n_perms);
    null_accs = perm_null(eval_H, eval_D, K_folds, n_perms);
    [cum_null_floor(i), cum_null_upper_999(i), cum_null_upper_99(i), cum_null_upper_95(i)] = null_ci(null_accs);
end
% Calculate Cumulative Accuracies
cum_acc = 1 - cum_err_mean;
cum_linear_acc = 1 - cum_linear_err_mean;

% Figure: mean greedy rank (top) and cumulative panel accuracy (bottom)
figure('Color', '#DAF2FB');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
x_axis = 1:n_markers;

% Marker labels (significance is shown by the bottom-tile marker colour)
custom_labels = strrep(stable_markers, '_', '\_');

% Top tile: mean greedy rank per marker (dots + black error bars)
ax1 = nexttile; hold on;
errorbar(x_axis, sorted_mean_ranks, sorted_rank_ci_lo, sorted_rank_ci_hi, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 0.75, 'CapSize', 0);
plot(x_axis, sorted_mean_ranks, '-ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
ylabel('mean rank');
set(gca, 'XTick', x_axis, 'XTickLabel', custom_labels, 'TickLabelInterpreter', 'tex', ...
    'TickDir', 'out', 'FontSize', 13, 'xlim', [0 n_markers+1], 'ylim', [0 n_markers+1], 'ytick', [1 n_markers], 'Color', 'w');
xtickangle(90);
box off;

% Bottom tile: cumulative accuracy vs the permutation null bands
ax2 = nexttile; hold on;
plot([0, n_markers+1], [0.5, 0.5], '-', 'LineWidth', 1, 'Color', [1 .85 .6]);  % chance level
rect_w = 0.6;
for k = 1:n_markers
    draw_ci_bands(k, rect_w/2, cum_null_floor(k), ...
        cum_null_upper_999(k), cum_null_upper_99(k), cum_null_upper_95(k));
end
% Linear accuracy (grey line); error bars are the 75% interval across CV reps
errorbar(x_axis, cum_linear_acc, cum_linear_acc - cum_linear_acc_ci_lo, cum_linear_acc_ci_hi - cum_linear_acc, ...
    '-o', 'Color', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerSize', 4, 'LineWidth', 1, 'CapSize', 0);
% Quadratic accuracy: black line and error bars throughout; markers filled
% black where significant (beats the 95% null), hollow (white face) otherwise
sig = cum_acc > cum_null_upper_95;
plot(x_axis, cum_acc, '-k', 'LineWidth', 1);
errorbar(x_axis(sig), cum_acc(sig), cum_acc(sig) - cum_acc_ci_lo(sig), cum_acc_ci_hi(sig) - cum_acc(sig), ...
    'ok', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', 1, 'CapSize', 0);
errorbar(x_axis(~sig), cum_acc(~sig), cum_acc(~sig) - cum_acc_ci_lo(~sig), cum_acc_ci_hi(~sig) - cum_acc(~sig), ...
    'ok', 'LineStyle', 'none', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 4, 'LineWidth', 1, 'CapSize', 0);
box off;
ylabel('combined separation');
set(gca, 'XTick', [], 'xlim', [0 n_markers+1], 'ylim',[.47 .9],'ytick', [0.5 .8], 'yticklabel', {'50%', '80%'}, ...
    'FontSize', 13, 'Color', 'w');
linkaxes([ax1, ax2], 'x');
sgtitle('Greedy marker ranking: H vs D separation');

%% 3. Visualize the greedy top-2 combination in 2D
% Uses the best 2-marker combination identified by the greedy search above,
% so this plot is dynamic rather than hard-coded to a specific pair.
results_top2 = classify_normals(H0{:, top2}, D0{:, top2}, 'input_type','samp','samp_balance',true,'prior_1',0.5,'samp_opt',0,'plotmode',0);

% Dynamic axis limits from the raw data range (5% padding on each side)
xdat = [H0{:, top2{1}}; D0{:, top2{1}}];
ydat = [H0{:, top2{2}}; D0{:, top2{2}}];
padx = 0.05 * (max(xdat,[],'omitnan') - min(xdat,[],'omitnan'));
pady = 0.05 * (max(ydat,[],'omitnan') - min(ydat,[],'omitnan'));
ax = [min(xdat,[],'omitnan')-padx, max(xdat,[],'omitnan')+padx, ...
      min(ydat,[],'omitnan')-pady, max(ydat,[],'omitnan')+pady];

figure('Color', '#DAF2FB');
hold on
fcontour(@(x,y) arrayfun(@(x0,y0) results_top2.post_1([x0; y0]), x, y), ax, ...
    'Fill', 'on', 'MeshDensity', 200, 'LevelList', linspace(0, 1, 100));

cb = colorbarpzn(0, 1, 'full', 0.5, 'colorP', [0.8 0.8 1], 'colorN', [1 0.8 0.8]);
cb.Ticks = [0, 0.5, 1];
cbTitle = title(cb, '$P(H | \mathbf{m})$');
cbTitle.Interpreter = 'latex';

plot(H0{:, top2{1}}, H0{:, top2{2}}, 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 4)
plot(D0{:, top2{1}}, D0{:, top2{2}}, 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 4)
plot_boundary(results_top2.norm_bd, 2, 'plot_type', 'line');
% plot the best linear boundary as well
plot_boundary(results_top2.norm_linear_bd, 2, 'plot_type', 'line', 'line_color', .5*[1 1 1]);
axis(ax)
xlabel(top2{1}, 'Interpreter', 'none')
ylabel(top2{2}, 'Interpreter', 'none')
set(gca, 'fontsize', 13, 'xtick', [], 'ytick', [])
title('Top-2 marker H vs D boundary');

% 5-fold CV test accuracy for the top-2 (quadratic & linear)
K_folds = 5; n_reps = 50;
cv_err = cv_classify_error(H0{:, top2}, D0{:, top2}, K_folds, n_reps);
linear_cv_err = cv_classify_error(H0{:, top2}, D0{:, top2}, K_folds, n_reps, 'linear', true);
fprintf('Top-2 (%s, %s) CV accuracy: QDA %.1f%%, linear %.1f%%\n', ...
    top2{1}, top2{2}, 100*(1-cv_err), 100*(1-linear_cv_err));

%% 4. Separate H0 vs D0 by individual markers
% Class-balanced H-vs-D accuracy for each biomarker, plus a joint top-2 column.
ind_samp_err=nan(1,n_markers);       % quadratic (QDA) CV error
ind_samp_acc_ci_lo=nan(1,n_markers);
ind_samp_acc_ci_hi=nan(1,n_markers);
ind_linear_err=nan(1,n_markers);     % linear-boundary CV error
ind_linear_acc_ci_lo=nan(1,n_markers);
ind_linear_acc_ci_hi=nan(1,n_markers);
bds=nan(n_markers,2);                % the two quadratic-boundary roots per marker
bds_alpha=nan(n_markers,2);          % per-root opacity for the boundary overlay
post_funcs=cell(1,n_markers);        % posterior P(H|x) handle per marker
K_folds = 5;
n_reps=50;

% Permutation null bounds per marker
n_perms = 1000;
null_floor = nan(1, n_markers);
null_upper_999 = nan(1, n_markers);
null_upper_99 = nan(1, n_markers);
null_upper_95 = nan(1, n_markers);

for i=1:n_markers
    marker=markers{i};
    fprintf('Evaluating %s (%d/%d) and running %d permutations...\n', marker, i, n_markers, n_perms);

    % Extract this marker's column by name (robust to column reordering)
    col_H = H0_z{:, marker};
    col_D = D0_z{:, marker};

    % Full-data fit; used for the boundary overlay in the top tile
    results=classify_normals(col_H, col_D, 'input_type','samp','samp_balance',true,'prior_1',0.5,'plotmode',0,'samp_opt',0);
    post_funcs{i} = results.post_1;

    % Quadratic boundary coefficients and their two roots
    q2 = results.norm_bd.q2;
    q1 = results.norm_bd.q1;
    q0 = results.norm_bd.q0;
    current_bds = sort(roots([q2 q1 q0]));
    if numel(current_bds) < 2, current_bds(end+1) = nan; end
    bds(i,:) = current_bds;

    % Tangent line to the quadratic boundary at root 1, scored as its own classifier
    r1 = current_bds(1);
    lin_bd1.q2 = 0;
    lin_bd1.q1 = 2 * q2 * r1 + q1;
    lin_bd1.q0 = -lin_bd1.q1 * r1;
    res1 = classify_normals(col_H, col_D, 'dom', lin_bd1, 'input_type', 'samp', 'samp_balance', true, 'prior_1', 0.5, 'plotmode', 0, 'samp_opt', 0);

    % Same for root 2
    r2 = current_bds(2);
    lin_bd2.q2 = 0;
    lin_bd2.q1 = 2 * q2 * r2 + q1;
    lin_bd2.q0 = -lin_bd2.q1 * r2;
    res2 = classify_normals(col_H, col_D, 'dom', lin_bd2, 'input_type', 'samp', 'samp_balance', true, 'prior_1', 0.5, 'plotmode', 0, 'samp_opt', 0);

    % Compute excess of better linear boundary & incremental quadratic accuracy
    acc1 = 1 - res1.samp_err;
    acc2 = 1 - res2.samp_err;
    acc_quad = 1 - results.samp_err;

    if acc1 >= acc2
        best_idx = 1;
        other_idx = 2;
        best_lin_acc = acc1;
    else
        best_idx = 2;
        other_idx = 1;
        best_lin_acc = acc2;
    end

    ex_best = max(0, best_lin_acc - 0.5);
    ex_incremental = max(0, acc_quad - best_lin_acc);

    if ex_best > 0
        % Scale so the best boundary is totally opaque (1.0)
        bds_alpha(i, best_idx) = 1.0;
        % The other boundary is proportional to the incremental accuracy
        bds_alpha(i, other_idx) = ex_incremental / ex_best;
    else
        % Fallback if neither performs better than random chance
        bds_alpha(i, 1) = 0;
        bds_alpha(i, 2) = 0;
    end

    % Cross-validated test error (quadratic, then linear boundary)
    [ind_samp_err(i), ~, rep_errs] = cv_classify_error(col_H, col_D, K_folds, n_reps);
    acc_ci = prctile(1 - rep_errs, [12.5 87.5]);
    ind_samp_acc_ci_lo(i) = acc_ci(1); ind_samp_acc_ci_hi(i) = acc_ci(2);
    [ind_linear_err(i), ~, rep_errs_linear] = cv_classify_error(col_H, col_D, K_folds, n_reps, 'linear', true);
    linear_acc_ci = prctile(1 - rep_errs_linear, [12.5 87.5]);
    ind_linear_acc_ci_lo(i) = linear_acc_ci(1); ind_linear_acc_ci_hi(i) = linear_acc_ci(2);

    % Permutation null (quadratic) -> 1-tailed CIs
    null_accs = perm_null(col_H, col_D, K_folds, n_perms);
    [null_floor(i), null_upper_999(i), null_upper_99(i), null_upper_95(i)] = null_ci(null_accs);
end

% Sort markers by CV test accuracy (best first) and reorder per-marker arrays
[~, idx] = sort(ind_samp_err);
ind_samp_err = ind_samp_err(idx);
ind_samp_acc_ci_lo = ind_samp_acc_ci_lo(idx);
ind_samp_acc_ci_hi = ind_samp_acc_ci_hi(idx);
ind_linear_err = ind_linear_err(idx);
ind_linear_acc_ci_lo = ind_linear_acc_ci_lo(idx);
ind_linear_acc_ci_hi = ind_linear_acc_ci_hi(idx);
ind_markers = markers(idx);
bds = bds(idx, :);
bds_alpha = bds_alpha(idx, :);
post_funcs = post_funcs(idx);

% ...and the per-marker null bounds
null_floor = null_floor(idx);
null_upper_999 = null_upper_999(idx);
null_upper_99 = null_upper_99(idx);
null_upper_95 = null_upper_95(idx);

% Reorder the z-scored data tables to the new marker order
H0_z = H0_z(:, ind_markers);
D0_z = D0_z(:, ind_markers);
Ha_z = Ha_z(:, ind_markers);
Hp_z = Hp_z(:, ind_markers);
Da_z = Da_z(:, ind_markers);
Dp_z = Dp_z(:, ind_markers);

% Joint "top 2" column: the best 2-marker combination from the greedy search
joint_markers = stable_markers(1:2);
fprintf('Classification top-2 (greedy): %s + %s\n', joint_markers{1}, joint_markers{2});
joint_H_data = H0{:, joint_markers};
joint_D_data = D0{:, joint_markers};
[joint_quad_err, ~, rep_errs_joint] = cv_classify_error(joint_H_data, joint_D_data, K_folds, n_reps);
joint_quad_acc_ci = prctile(1 - rep_errs_joint, [12.5 87.5]);
joint_quad_acc_ci_lo = joint_quad_acc_ci(1); joint_quad_acc_ci_hi = joint_quad_acc_ci(2);
[joint_lin_err, ~, rep_errs_joint_linear] = cv_classify_error(joint_H_data, joint_D_data, K_folds, n_reps, 'linear', true);
joint_lin_acc_ci = prctile(1 - rep_errs_joint_linear, [12.5 87.5]);
joint_lin_acc_ci_lo = joint_lin_acc_ci(1); joint_lin_acc_ci_hi = joint_lin_acc_ci(2);

% Permutation null -> 1-tailed CIs for the joint model
joint_null_accs = perm_null(joint_H_data, joint_D_data, K_folds, n_perms);
[joint_null_floor, joint_null_upper_999, joint_null_upper_99, joint_null_upper_95] = null_ci(joint_null_accs);

% Joint decision values mapped to posterior P(H|x), rescaled to [-3, 3] so the
% joint column shares the y-axis of the individual-marker violins.
results_joint = classify_normals(joint_H_data, joint_D_data, 'input_type','samp','samp_balance',true,'prior_1',0.5,'samp_opt',0,'plotmode',0);
joint_dv_H0 = results_joint.samp_dv{1};
joint_dv_D0 = results_joint.samp_dv{2};
post_joint_dv_H0 = 6 * sigmoid(joint_dv_H0) - 3;
post_joint_dv_D0 = 6 * sigmoid(joint_dv_D0) - 3;
post_bd_0          = 0;

% =========================================================================
% --- CREATE TEMPORARY PLOTTING VARIABLES (DO NOT OVERWRITE MASTER ARRAYS) ---
% =========================================================================
plot_samp_err        = [joint_quad_err, ind_samp_err(:)'];
plot_samp_acc_ci_lo  = [joint_quad_acc_ci_lo, ind_samp_acc_ci_lo(:)'];
plot_samp_acc_ci_hi  = [joint_quad_acc_ci_hi, ind_samp_acc_ci_hi(:)'];
plot_linear_err      = [joint_lin_err, ind_linear_err(:)'];
plot_linear_acc_ci_lo = [joint_lin_acc_ci_lo, ind_linear_acc_ci_lo(:)'];
plot_linear_acc_ci_hi = [joint_lin_acc_ci_hi, ind_linear_acc_ci_hi(:)'];
plot_markers       = [{'top 2'}, ind_markers(:)'];
plot_post_funcs    = [{[]}, post_funcs(:)']; % Placeholder for joint marker

% Plotting arrays for all CIs
plot_null_floor     = [joint_null_floor, null_floor(:)'];
plot_null_upper_999 = [joint_null_upper_999, null_upper_999(:)'];
plot_null_upper_99  = [joint_null_upper_99, null_upper_99(:)'];
plot_null_upper_95  = [joint_null_upper_95, null_upper_95(:)'];
plot_bds           = [[post_bd_0, post_bd_0]; bds];
plot_bds_alpha     = [[1, 0]; bds_alpha];

% Per-column data for the top tile (joint column first, then markers)
plot_H_data = [{post_joint_dv_H0}, num2cell(H0_z{:,:}, 1)];
plot_D_data = [{post_joint_dv_D0}, num2cell(D0_z{:,:}, 1)];
n_total = length(plot_markers);
x = 1:n_total;
x_H = x - .14;   % healthy points sit slightly left of centre, depressed slightly right
x_D = x + .14;

figure('Color', '#DAF2FB');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Top tile: means ± SD (H vs D) ---
nexttile; hold on
% Define pure colors
pure_blue = [0 0 1];
pure_red  = [1 0 0];
jit_width = 0.12; % Controls how wide the jitter spreads

for i_marker = 1:n_total
    % Background band: a vertical P(H) colour gradient behind each column
    xSpan_col = x(i_marker) + [-0.28, 0.28];
    y_grid = linspace(-4, 5, 200)';

    if i_marker == 1
        % For joint marker, Y is already scaled P(H): y = 6*P - 3  =>  P = (y+3)/6
        P_H = max(0, min(1, (y_grid + 3) / 6));
    else
        % For individual markers, use the saved post_1 function
        func = plot_post_funcs{i_marker};
        P_H = arrayfun(@(v) func(v), y_grid);
    end

    % Map Probability to RGB (Light Blue -> White -> Light Red)
    RGB = zeros(length(y_grid), 3);
    c_P = [0.8 0.8 1]; % Light Blue (Healthy)
    c_Z = [1 1 1];     % White
    c_N = [1 0.8 0.8]; % Light Red (Depressed)
    for j = 1:length(y_grid)
        if P_H(j) >= 0.5
            f = (P_H(j) - 0.5) * 2;
            RGB(j,:) = f * c_P + (1 - f) * c_Z;
        else
            f = P_H(j) * 2;
            RGB(j,:) = f * c_Z + (1 - f) * c_N;
        end
    end

    % Draw the vertical color gradient using surf
    [X_surf, Y_surf] = meshgrid(xSpan_col, y_grid);
    Z_surf = zeros(size(X_surf)); % Placed at bottom
    C_surf = repmat(reshape(RGB, [length(y_grid), 1, 3]), [1, 2, 1]);
    surf(X_surf, Y_surf, Z_surf, C_surf, 'EdgeColor', 'none');

    % This column's values, NaNs removed
    yH_val = plot_H_data{i_marker};
    yD_val = plot_D_data{i_marker};
    yH_val = yH_val(~isnan(yH_val));
    yD_val = yD_val(~isnan(yD_val));

    muH = mean(yH_val); sdH = std(yH_val);
    muD = mean(yD_val); sdD = std(yD_val);

    % Mean +/- SD bars
    plot([x_H(i_marker), x_H(i_marker)], [muH-sdH, muH+sdH], '-', 'Color', [pure_blue .3], 'LineWidth', 3);
    plot([x_D(i_marker), x_D(i_marker)], [muD-sdD, muD+sdD], '-', 'Color', [pure_red .3], 'LineWidth', 3);

    % Jittered scatter: healthy (blue), depressed (red)
    xH_jit = x_H(i_marker) + (rand(size(yH_val)) - 0.5) * jit_width;
    xD_jit = x_D(i_marker) + (rand(size(yD_val)) - 0.5) * jit_width;
    scatter(xH_jit, yH_val, 2, pure_blue, 'filled', 'MarkerEdgeColor', 'none');
    scatter(xD_jit, yD_val, 2, pure_red, 'filled', 'MarkerEdgeColor', 'none');

    % Boundary lines, opacity = how much that boundary contributes (clamped to [0,1])
    xSpan = x(i_marker)+.3*[-1 1];
    a1 = max(0, min(1, plot_bds_alpha(i_marker, 1)));
    a2 = max(0, min(1, plot_bds_alpha(i_marker, 2)));

    plot(xSpan, plot_bds(i_marker,1)*[1 1], '-', 'color', [0 0 0 a1], 'LineWidth', 1);
    plot(xSpan, plot_bds(i_marker,2)*[1 1], '-', 'color', [0 0 0 a2], 'LineWidth', 1);
end

% Tick labels: black if significant (beats 95% null), grey otherwise;
% joint column is bold only when significant
custom_labels = cell(1, n_total);
for k = 1:n_total
    safe_name = strrep(plot_markers{k}, '_', '\_');
    if (1 - plot_samp_err(k)) > plot_null_upper_95(k)
        if k == 1
            custom_labels{k} = sprintf('\\bf{\\color{black}%s}', safe_name);
        else
            custom_labels{k} = sprintf('\\color{black}%s', safe_name);
        end
    else
        custom_labels{k} = sprintf('\\color[rgb]{0.7,0.7,0.7}%s', safe_name);
    end
end

set(gca, 'XTick', x, 'XTickLabel', custom_labels, 'TickLabelInterpreter', 'tex', ...
    'xlim', [0 n_total+1], 'ylim', [-4 5], 'YTick', [-4 0 4], 'TickDir', 'out', ...
    'fontsize', 13, 'Color', 'w', 'Layer', 'top'); % 'top' keeps the axis frame above the surf band
xtickangle(90);
ylabel('marker values');
cb = colorbarpzn(0, 1, 'full', 0.5, 'colorP', [0.8 0.8 1], 'colorN', [1 0.8 0.8]);
cb.Ticks = [0, 0.5, 1];
cbTitle = title(cb, '$P(H | \mathbf{m})$');
cbTitle.Interpreter = 'latex';

% --- Bottom tile: sorted individual test accuracy ---
nexttile; hold on;
plot([0, n_total+1], [0.5, 0.5], '-', 'LineWidth', .5, 'Color', [1 .85 .6]);  % chance level

% Permutation null bands
rect_w_bot = 0.6;
for k = 1:n_total
    draw_ci_bands(k, rect_w_bot/2, plot_null_floor(k), ...
        plot_null_upper_999(k), plot_null_upper_99(k), plot_null_upper_95(k));
end

% Linear accuracy (grey line); error bars are the 75% interval across CV reps
lacc = 1 - plot_linear_err;
errorbar(x, lacc, lacc - plot_linear_acc_ci_lo, plot_linear_acc_ci_hi - lacc, '-o', 'Color', .5*[1 1 1], ...
    'MarkerFaceColor', .5*[1 1 1], 'MarkerSize', 4, 'CapSize', 0, 'LineWidth', 1);
% Quadratic accuracy: black line and error bars; markers filled black where
% significant (beats the 95% null), hollow (white face) otherwise
qacc = 1 - plot_samp_err;
sig = qacc > plot_null_upper_95;
plot(x, qacc, '-k', 'LineWidth', 1);
errorbar(x(sig), qacc(sig), qacc(sig) - plot_samp_acc_ci_lo(sig), plot_samp_acc_ci_hi(sig) - qacc(sig), ...
    'ok', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'CapSize', 0, 'LineWidth', 1);
errorbar(x(~sig), qacc(~sig), qacc(~sig) - plot_samp_acc_ci_lo(~sig), plot_samp_acc_ci_hi(~sig) - qacc(~sig), ...
    'ok', 'LineStyle', 'none', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 4, 'CapSize', 0, 'LineWidth', 1);

set(gca, 'XAxisLocation', 'top', ...
    'TickDir', 'out', ...
    'XTick', x, ...
    'XTickLabel', [], ...
    'xlim', [0 n_total+1], ...
    'ylim', [.35 .9], ...
    'ytick', [.5 .8], ...
    'yticklabel', {'50%','80%'}, ...
    'fontsize', 13, 'Color', 'w');
box off
ylabel('separation');
sgtitle('Per-marker H vs D separation');

%% 5. Greedily accumulate restoration markers (Depressed)
num_runs = 100;
K_folds = 5;
n_perms = 1000;

% Align each depressed baseline subject with their post-treatment row
[~, actual_bA, aA] = intersect(D0_full.id, Da_full.id, 'stable');
[~, actual_bP, aP] = intersect(D0_full.id, Dp_full.id, 'stable');
dep_preA = D0_z{actual_bA, :}; dep_postA = Da_z{aA, :};
dep_preP = D0_z{actual_bP, :}; dep_postP = Dp_z{aP, :};

% Greedy search; reference is the healthy baseline (recovery moves toward it)
stable_markers_rec = greedy_restoration(dep_preA, dep_postA, dep_preP, dep_postP, ...
    H0_z{:, :}, true, ind_markers, num_runs, K_folds, n_perms);
[~, dep_top2_idx] = ismember(stable_markers_rec(1:2), ind_markers);

%% 6. Individual marker restorations of Depressed
fprintf('Depressed restoration top-2 (greedy): %s + %s\n', stable_markers_rec{1}, stable_markers_rec{2});

% Joint column uses the greedy top-2; colours = dark/light red (Aya/Placebo)
% Reuses the baseline/post-treatment alignment computed for the greedy search above.
restoration_figure(dep_preA, dep_postA, dep_preP, dep_postP, ...
    H0_z{:, :}, D0_z{:, :}, true, dep_top2_idx, [.8 0 0], [1 .6 .6], ind_markers);

%% 7. Treatment trajectories in top-2 baseline-classification space.
% Visualizes, per group and arm, how subjects move relative to the fixed
% baseline (H0 vs D0) boundary in the top-2 marker space (results_top2 from
% Section 3, i.e. top2 = {corti_sal, crp}) between baseline and post-treatment.

% Raw (unscaled) baseline vs post-treatment values in the top-2 marker space,
% aligned per subject id (same alignment logic as Section 5)
[~, bA_idx, aA_idx] = intersect(D0_full.id, Da_full.id, 'stable');
[~, bP_idx, aP_idx] = intersect(D0_full.id, Dp_full.id, 'stable');
pre_A = D0{bA_idx, top2};  post_A = Da{aA_idx, top2};
pre_P = D0{bP_idx, top2};  post_P = Dp{aP_idx, top2};

[~, bA_idx_H, aA_idx_H] = intersect(H0_full.id, Ha_full.id, 'stable');
[~, bP_idx_H, aP_idx_H] = intersect(H0_full.id, Hp_full.id, 'stable');
H_pre_A = H0{bA_idx_H, top2};  H_post_A = Ha{aA_idx_H, top2};
H_pre_P = H0{bP_idx_H, top2};  H_post_P = Hp{aP_idx_H, top2};

% Posterior P(healthy | m) at baseline and post-treatment, from the fixed
% baseline boundary, for depressed (pD_*) and healthy (pH_*) subjects
pD_pre_A  = arrayfun(@(x,y) results_top2.post_1([x;y]), pre_A(:,1),  pre_A(:,2));
pD_post_A = arrayfun(@(x,y) results_top2.post_1([x;y]), post_A(:,1), post_A(:,2));
pD_pre_P  = arrayfun(@(x,y) results_top2.post_1([x;y]), pre_P(:,1),  pre_P(:,2));
pD_post_P = arrayfun(@(x,y) results_top2.post_1([x;y]), post_P(:,1), post_P(:,2));
pH_pre_A  = arrayfun(@(x,y) results_top2.post_1([x;y]), H_pre_A(:,1),  H_pre_A(:,2));
pH_post_A = arrayfun(@(x,y) results_top2.post_1([x;y]), H_post_A(:,1), H_post_A(:,2));
pH_pre_P  = arrayfun(@(x,y) results_top2.post_1([x;y]), H_pre_P(:,1),  H_pre_P(:,2));
pH_post_P = arrayfun(@(x,y) results_top2.post_1([x;y]), H_post_P(:,1), H_post_P(:,2));

% Treatment-trajectory plots: baseline -> post-treatment, over the baseline
% H/D boundary, dot at the tip (the 'after' point). The other group is shown
% static at its baseline spot (filled dot) for context.
pure_blue = [0 0 1]; pure_red = [1 0 0];
all_x = [pre_A(:,1); post_A(:,1); pre_P(:,1); post_P(:,1); H_pre_A(:,1); H_post_A(:,1); H_pre_P(:,1); H_post_P(:,1); H0{:,top2{1}}; D0{:,top2{1}}];
all_y = [pre_A(:,2); post_A(:,2); pre_P(:,2); post_P(:,2); H_pre_A(:,2); H_post_A(:,2); H_pre_P(:,2); H_post_P(:,2); H0{:,top2{2}}; D0{:,top2{2}}];
padx = 0.05 * (max(all_x,[],'omitnan') - min(all_x,[],'omitnan'));
pady = 0.05 * (max(all_y,[],'omitnan') - min(all_y,[],'omitnan'));
ax65 = [min(all_x,[],'omitnan')-padx, max(all_x,[],'omitnan')+padx, ...
        min(all_y,[],'omitnan')-pady, max(all_y,[],'omitnan')+pady];

% 2x2 grid (Depressed/Healthy x Ayahuasca/Placebo), as in Section 10: each
% quadrant shows one group+arm's trajectories, with the other group shown
% static at baseline for context, and its own colorbar-overlaid P(H) trace.
figure('Color', '#DAF2FB');
sgtitle('Treatment trajectories');

light_red = [1 .5 .5]; light_blue = [.5 .5 1];

subplot(2,2,1);
cb1 = plot_trajectory_subplot(ax65, results_top2, H0{:,top2{1}}, H0{:,top2{2}}, pure_blue, ...
    pre_A, post_A, light_red, pure_red, top2, true);
title('Ayahuasca');

subplot(2,2,2);
cb2 = plot_trajectory_subplot(ax65, results_top2, H0{:,top2{1}}, H0{:,top2{2}}, pure_blue, ...
    pre_P, post_P, light_red, pure_red, top2, false);
title('Placebo');

subplot(2,2,3);
cb3 = plot_trajectory_subplot(ax65, results_top2, D0{:,top2{1}}, D0{:,top2{2}}, pure_red, ...
    H_pre_A, H_post_A, light_blue, pure_blue, top2, false);

subplot(2,2,4);
cb4 = plot_trajectory_subplot(ax65, results_top2, D0{:,top2{1}}, D0{:,top2{2}}, pure_red, ...
    H_pre_P, H_post_P, light_blue, pure_blue, top2, false);

% Force layout to settle before reading colorbar positions, then overlay each
% quadrant's own P(H) trace on its colorbar
drawnow;
plot_colorbar_trace(cb1, pD_pre_A, pD_post_A, [1 .7 .7], pure_red, pure_red, 1.5, pure_red);
plot_colorbar_trace(cb2, pD_pre_P, pD_post_P, [1 .7 .7], pure_red, pure_red, 1.5, pure_red);
plot_colorbar_trace(cb3, pH_pre_A, pH_post_A, [.75 .75 1], 'none', [.3 .3 1], 2, pure_blue);
plot_colorbar_trace(cb4, pH_pre_P, pH_post_P, [.75 .75 1], 'none', [.3 .3 1], 2, pure_blue);

%% 8. Greedily accumulate restoration markers (Healthy)
num_runs = 100;
K_folds = 5;
n_perms = 1000;

% Align each healthy baseline subject with their post-treatment row
[~, actual_bA_H, aA_H] = intersect(H0_full.id, Ha_full.id, 'stable');
[~, actual_bP_H, aP_H] = intersect(H0_full.id, Hp_full.id, 'stable');
heal_preA = H0_z{actual_bA_H, :}; heal_postA = Ha_z{aA_H, :};
heal_preP = H0_z{actual_bP_H, :}; heal_postP = Hp_z{aP_H, :};

% Greedy search; reference is the depressed baseline
stable_markers_rec = greedy_restoration(heal_preA, heal_postA, heal_preP, heal_postP, ...
    D0_z{:, :}, false, ind_markers, num_runs, K_folds, n_perms);
[~, heal_top2_idx] = ismember(stable_markers_rec(1:2), ind_markers);

%% 9. Individual marker restorations of Healthy
fprintf('Healthy restoration top-2 (greedy): %s + %s\n', stable_markers_rec{1}, stable_markers_rec{2});

% Joint column uses the greedy top-2; colours = dark/light blue (Aya/Placebo)
% Reuses the baseline/post-treatment alignment computed for the greedy search above.
restoration_figure(heal_preA, heal_postA, heal_preP, heal_postP, ...
    H0_z{:, :}, D0_z{:, :}, false, heal_top2_idx, [0 .2 .6], [.6 .8 1], ind_markers);

%% 10. All Treatment Trajectories
H0 = before(strcmpi(before.group, 'H'), :);
D0 = before(strcmpi(before.group, 'D'), :);
Ha = after(strcmpi(after.group, 'H') & strcmpi(after.treatment, 'Ayahuasca'), :);
Hp = after(strcmpi(after.group, 'H') & strcmpi(after.treatment, 'Placebo'), :);
Da = after(strcmpi(after.group, 'D') & strcmpi(after.treatment, 'Ayahuasca'), :);
Dp = after(strcmpi(after.group, 'D') & strcmpi(after.treatment, 'Placebo'), :);
marker_1 = 'crp';
marker_2 = 'creatinine';
ix = find(strcmpi(ind_markers, marker_1));
iy = find(strcmpi(ind_markers, marker_2));

% Build per-subject baseline->post point pairs (depressed & healthy, both arms)
D_Ax0=[]; D_Ay0=[]; D_Ax1=[]; D_Ay1=[]; D_Px0=[]; D_Py0=[]; D_Px1=[]; D_Py1=[];
for i = 1:size(D0,1)
    x0 = D0_z{i,ix}; y0 = D0_z{i,iy};
    if ~isnan(x0) && ~isnan(y0)
        if strcmpi(D0.treatment{i}, 'Ayahuasca')
            ia = find(strcmpi(Da.id, D0.id{i}));
            if ~isempty(ia) && ~isnan(Da_z{ia,ix}) && ~isnan(Da_z{ia,iy})
                D_Ax0(end+1)=x0; D_Ay0(end+1)=y0; D_Ax1(end+1)=Da_z{ia,ix}; D_Ay1(end+1)=Da_z{ia,iy};
            end
        else
            ip = find(strcmpi(Dp.id, D0.id{i}));
            if ~isempty(ip) && ~isnan(Dp_z{ip,ix}) && ~isnan(Dp_z{ip,iy})
                D_Px0(end+1)=x0; D_Py0(end+1)=y0; D_Px1(end+1)=Dp_z{ip,ix}; D_Py1(end+1)=Dp_z{ip,iy};
            end
        end
    end
end

H_Ax0=[]; H_Ay0=[]; H_Ax1=[]; H_Ay1=[]; H_Px0=[]; H_Py0=[]; H_Px1=[]; H_Py1=[];
for i = 1:size(H0,1)
    x0 = H0_z{i,ix}; y0 = H0_z{i,iy};
    if ~isnan(x0) && ~isnan(y0)
        if strcmpi(H0.treatment{i}, 'Ayahuasca')
            ia = find(strcmpi(Ha.id, H0.id{i}));
            if ~isempty(ia) && ~isnan(Ha_z{ia,ix}) && ~isnan(Ha_z{ia,iy})
                H_Ax0(end+1)=x0; H_Ay0(end+1)=y0; H_Ax1(end+1)=Ha_z{ia,ix}; H_Ay1(end+1)=Ha_z{ia,iy};
            end
        else
            ip = find(strcmpi(Hp.id, H0.id{i}));
            if ~isempty(ip) && ~isnan(Hp_z{ip,ix}) && ~isnan(Hp_z{ip,iy})
                H_Px0(end+1)=x0; H_Py0(end+1)=y0; H_Px1(end+1)=Hp_z{ip,ix}; H_Py1(end+1)=Hp_z{ip,iy};
            end
        end
    end
end

% Baseline H-vs-D boundary and global axis limits
dist_1 = [[H_Ax0, H_Px0]', [H_Ay0, H_Py0]']; % Healthy baseline
dist_2 = [[D_Ax0, D_Px0]', [D_Ay0, D_Py0]']; % Depressed baseline
results = classify_normals(dist_1, dist_2, 'input_type', 'samp', 'samp_balance', true, 'prior_1', 0.5, 'samp_opt', 0, 'plotmode', 0);
baseline_bd = results.norm_bd;

% Calculate dynamic axes limits over ALL points
X_all = [D_Ax0 D_Ax1 D_Px0 D_Px1 H_Ax0 H_Ax1 H_Px0 H_Px1];
Y_all = [D_Ay0 D_Ay1 D_Py0 D_Py1 H_Ay0 H_Ay1 H_Py0 H_Py1];
ax = [min(X_all)-0.5, max(X_all)+0.5, min(Y_all)-0.5, max(Y_all)+0.5];

% Figure setup (2x2: depressed/healthy x Ayahuasca/placebo)
figure('Color', '#DAF2FB');
sgtitle('Treatment trajectories');
main_ax = zeros(1, 4);
cb_arr = gobjects(1, 4); % Pre-allocate for colorbar handles

% Generate unified background contours for all 4 subplots
for idx = 1:4
    main_ax(idx) = subplot(2, 2, idx); hold on;
    fcontour(@(x,y) arrayfun(@(x0,y0) results.post_1([x0; y0]), x, y), ax, ...
        'Fill', 'on', 'MeshDensity', 200, 'LevelList', linspace(0, 1, 100));
    cb_arr(idx) = colorbarpzn(0, 1, 'full', 0.5, 'colorP', [0.8 0.8 1], 'colorN', [1 0.8 0.8]);
    plot_boundary(baseline_bd, 2, 'plot_type', 'line', 'line_color', [0 0 0]);
    cb_arr(idx).Ticks = 0:.25:1;
    cb_arr(idx).TickLength = 0.03;
    if idx==1
        xlabel(marker_1, 'Interpreter', 'none'); ylabel(marker_2, 'Interpreter', 'none');
        cbTitle = title(cb_arr(idx), '$P(H | \mathbf{m})$');
        cbTitle.Interpreter = 'latex';
        cb_arr(idx).TickLabels = {'0%','25%','50%','75%','100%'};
    else
        cb_arr(idx).TickLabels = {};
    end
    if idx==1, title('Ayahuasca'); elseif idx==2, title('Placebo'); end
    box on; axis(ax); axis square;
    set(gca, 'XTick', [], 'YTick', [], 'fontsize', 13, 'Color', 'w');
end

% CRITICAL: Force layout to settle before reading colorbar positions
drawnow;
cb_ax = zeros(1, 4);
for idx = 1:4
    cb_ax(idx) = axes('Position', cb_arr(idx).Position, 'Color', 'none', ...
        'XLim', [0 1], 'YLim', [0 1], 'XTick', [], 'YTick', [], ...
        'XColor', 'none', 'YColor', 'none', 'Box', 'off');
    uistack(cb_ax(idx), 'top'); % Force overlay axes above colorbar
end

% Setup Subplot 1 (Depressed Ayahuasca)
axes(main_ax(1));
plot([H_Ax0, H_Px0], [H_Ay0, H_Py0], 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 3, 'MarkerEdgeColor', 'none'); % static healthy baseline (shown for context)
hDA_traces = plot(nan, nan, '-', 'Color', [1 .5 .5], 'linewidth', 1.25);
hDA = plot(D_Ax0, D_Ay0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3);
axes(cb_ax(1)); hold on;
N_DA = length(D_Ax0);
jitter = linspace(.1,.9,N_DA);
Prob_DA0 = arrayfun(@(x,y) results.post_1([x; y]), D_Ax0, D_Ay0);
hDA_cb_traces = plot(reshape([jitter; jitter; nan(1, N_DA)], [], 1), nan(3*N_DA, 1), '-', 'color', [1 .7 .7], 'LineWidth', .5);
hDA_cb_pts = plot(jitter, Prob_DA0, 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 1.5);
hDA_cb_mean_trace = plot([0.5; 0.5; nan], nan(3, 1), 'r-', 'LineWidth', 1.5);
hDA_cb_mean_pt = plot(0.5, mean(Prob_DA0), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);

% Setup Subplot 2 (Depressed Placebo)
axes(main_ax(2));
plot([H_Ax0, H_Px0], [H_Ay0, H_Py0], 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 3, 'MarkerEdgeColor', 'none'); % static healthy baseline (shown for context)
hDP_traces = plot(nan, nan, '-', 'Color', [1 .5 .5], 'linewidth', 1.25);
hDP = plot(D_Px0, D_Py0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3);
axes(cb_ax(2)); hold on;
N_DP = length(D_Px0);
jitter = linspace(.1,.9,N_DP);
Prob_DP0 = arrayfun(@(x,y) results.post_1([x; y]), D_Px0, D_Py0);
hDP_cb_traces = plot(reshape([jitter; jitter; nan(1, N_DP)], [], 1), nan(3*N_DP, 1), '-', 'color', [1 .7 .7], 'LineWidth', .5);
hDP_cb_pts = plot(jitter, Prob_DP0, 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 1.5);
hDP_cb_mean_trace = plot([0.5; 0.5; nan], nan(3, 1), 'r-', 'LineWidth', 1.5);
hDP_cb_mean_pt = plot(0.5, mean(Prob_DP0), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);

% Setup Subplot 3 (Healthy Ayahuasca)
axes(main_ax(3));
plot([D_Ax0, D_Px0], [D_Ay0, D_Py0], 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3, 'MarkerEdgeColor', 'none'); % static depressed baseline (shown for context)
hHA_traces = plot(nan, nan, '-', 'Color', [.5 .5 1], 'linewidth', 1.25);
hHA = plot(H_Ax0, H_Ay0, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 3);
axes(cb_ax(3)); hold on;
N_HA = length(H_Ax0);
jitter = linspace(.1,.9,N_HA);
Prob_HA0 = arrayfun(@(x,y) results.post_1([x; y]), H_Ax0, H_Ay0);
hHA_cb_traces = plot(reshape([jitter; jitter; nan(1, N_HA)], [], 1), nan(3*N_HA, 1), '-', 'color', [.75 .75 1], 'LineWidth', .5);
hHA_cb_pts = plot(jitter, Prob_HA0, 'o', 'Color', 'none', 'MarkerFaceColor', [.3 .3 1], 'MarkerSize', 2);
hHA_cb_mean_trace = plot([0.5; 0.5; nan], nan(3, 1), 'b-', 'LineWidth', 1.5);
hHA_cb_mean_pt = plot(0.5, mean(Prob_HA0), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

% Setup Subplot 4 (Healthy Placebo)
axes(main_ax(4));
plot([D_Ax0, D_Px0], [D_Ay0, D_Py0], 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3, 'MarkerEdgeColor', 'none'); % static depressed baseline (shown for context)
hHP_traces = plot(nan, nan, '-', 'Color', [.5 .5 1], 'linewidth', 1.25);
hHP = plot(H_Px0, H_Py0, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 3);
axes(cb_ax(4)); hold on;
N_HP = length(H_Px0);
jitter = linspace(.1,.9,N_HP);
Prob_HP0 = arrayfun(@(x,y) results.post_1([x; y]), H_Px0, H_Py0);
hHP_cb_traces = plot(reshape([jitter; jitter; nan(1, N_HP)], [], 1), nan(3*N_HP, 1), '-', 'color', [.75 .75 1], 'LineWidth', .5);
hHP_cb_pts = plot(jitter, Prob_HP0, 'o', 'Color', 'none', 'MarkerFaceColor', [.3 .3 1], 'MarkerSize', 2);
hHP_cb_mean_trace = plot([0.5; 0.5; nan], nan(3, 1), 'b-', 'LineWidth', 1.5);
hHP_cb_mean_pt = plot(0.5, mean(Prob_HP0), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

% Animation loop: interpolate each point from baseline to post
gif_filename = 'treatment_trajectories.gif';
t_vals = linspace(0, 1, 50);

for f = 1:length(t_vals)
    t = t_vals(f);
    if ~ishandle(hDA), break; end

    % Update Depressed Ayahuasca
    c_DAx = D_Ax0 + t*(D_Ax1-D_Ax0); c_DAy = D_Ay0 + t*(D_Ay1-D_Ay0);
    set(hDA, 'XData', c_DAx, 'YData', c_DAy);
    set(hDA_traces, 'XData', reshape([D_Ax0; c_DAx; nan(1, N_DA)], [], 1), 'YData', reshape([D_Ay0; c_DAy; nan(1, N_DA)], [], 1));
    Prob_DA = arrayfun(@(x,y) results.post_1([x; y]), c_DAx, c_DAy);
    set(hDA_cb_pts, 'YData', Prob_DA);
    set(hDA_cb_traces, 'YData', reshape([Prob_DA0; Prob_DA; nan(1, N_DA)], [], 1));
    set(hDA_cb_mean_pt, 'YData', mean(Prob_DA));
    set(hDA_cb_mean_trace, 'YData', [mean(Prob_DA0); mean(Prob_DA); nan]);

    % Update Depressed Placebo
    c_DPx = D_Px0 + t*(D_Px1-D_Px0); c_DPy = D_Py0 + t*(D_Py1-D_Py0);
    set(hDP, 'XData', c_DPx, 'YData', c_DPy);
    set(hDP_traces, 'XData', reshape([D_Px0; c_DPx; nan(1, N_DP)], [], 1), 'YData', reshape([D_Py0; c_DPy; nan(1, N_DP)], [], 1));
    Prob_DP = arrayfun(@(x,y) results.post_1([x; y]), c_DPx, c_DPy);
    set(hDP_cb_pts, 'YData', Prob_DP);
    set(hDP_cb_traces, 'YData', reshape([Prob_DP0; Prob_DP; nan(1, N_DP)], [], 1));
    set(hDP_cb_mean_pt, 'YData', mean(Prob_DP));
    set(hDP_cb_mean_trace, 'YData', [mean(Prob_DP0); mean(Prob_DP); nan]);

    % Update Healthy Ayahuasca
    c_HAx = H_Ax0 + t*(H_Ax1-H_Ax0); c_HAy = H_Ay0 + t*(H_Ay1-H_Ay0);
    set(hHA, 'XData', c_HAx, 'YData', c_HAy);
    set(hHA_traces, 'XData', reshape([H_Ax0; c_HAx; nan(1, N_HA)], [], 1), 'YData', reshape([H_Ay0; c_HAy; nan(1, N_HA)], [], 1));
    Prob_HA = arrayfun(@(x,y) results.post_1([x; y]), c_HAx, c_HAy);
    set(hHA_cb_pts, 'YData', Prob_HA);
    set(hHA_cb_traces, 'YData', reshape([Prob_HA0; Prob_HA; nan(1, N_HA)], [], 1));
    set(hHA_cb_mean_pt, 'YData', mean(Prob_HA));
    set(hHA_cb_mean_trace, 'YData', [mean(Prob_HA0); mean(Prob_HA); nan]);

    % Update Healthy Placebo
    c_HPx = H_Px0 + t*(H_Px1-H_Px0); c_HPy = H_Py0 + t*(H_Py1-H_Py0);
    set(hHP, 'XData', c_HPx, 'YData', c_HPy);
    set(hHP_traces, 'XData', reshape([H_Px0; c_HPx; nan(1, N_HP)], [], 1), 'YData', reshape([H_Py0; c_HPy; nan(1, N_HP)], [], 1));
    Prob_HP = arrayfun(@(x,y) results.post_1([x; y]), c_HPx, c_HPy);
    set(hHP_cb_pts, 'YData', Prob_HP);
    set(hHP_cb_traces, 'YData', reshape([Prob_HP0; Prob_HP; nan(1, N_HP)], [], 1));
    set(hHP_cb_mean_pt, 'YData', mean(Prob_HP));
    set(hHP_cb_mean_trace, 'YData', [mean(Prob_HP0); mean(Prob_HP); nan]);

    drawnow;

    % GIF Export Logic
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    if f == 1
        delay = 1.0;
    elseif f == length(t_vals)
        delay = 1.0;
    else
        delay = 0.02;
    end

    if f == 1
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', delay);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay);
    end
end

%% 11. Baseline MADRS vs baseline marker values
% Only the depressed group filled out the MADRS at baseline (healthy subjects
% have none), so this is depressed-only -- both arms combined, since treatment
% assignment hasn't happened yet at baseline. Markers are z-scored against the
% healthy-baseline mean/SD (D0_z, in ind_markers order); correlation is
% invariant to that affine rescaling.
%
% Per marker: cross-validated R^2 of baseline-marker -> baseline-MADRS, and a
% permutation null distribution of R^2 (same LOOCV + permutation scheme as
% Section 13).
madrs0 = D0_full.MADRS;

cv_r2_vals        = nan(1, n_markers);
cv_p_vals         = nan(1, n_markers);
null_r2_lower_999 = nan(1, n_markers);  % common floor for the CI bands
null_r2_upper_999 = nan(1, n_markers);
null_r2_upper_99  = nan(1, n_markers);
null_r2_upper_95  = nan(1, n_markers);
n_perms           = 1000;

for i = 1:n_markers
    v = ~isnan(D0_z{:, i}) & ~isnan(madrs0);
    all_x = D0_z{v, i};
    all_y = madrs0(v);

    % Require enough points to do LOOCV safely
    if length(all_x) >= 5
        % LOOCV R^2 from out-of-fold predictions
        cvp = cvpartition(length(all_y), 'LeaveOut');
        y_pred_true = nan(size(all_y));
        for fold = 1:cvp.NumTestSets
            train_idx = training(cvp, fold);
            test_idx = test(cvp, fold);
            p_cv = polyfit(all_x(train_idx), all_y(train_idx), 1);
            y_pred_true(test_idx) = polyval(p_cv, all_x(test_idx));
        end
        SS_res_true = sum((all_y - y_pred_true).^2);
        SS_tot_true = sum((all_y - mean(all_y)).^2);
        true_r2 = 1 - (SS_res_true / SS_tot_true);
        cv_r2_vals(i) = true_r2;

        % Permutation null: shuffle baseline MADRS and recompute LOOCV R^2
        null_r2_dist = nan(1, n_perms);
        for p_idx = 1:n_perms
            shuffled_y = all_y(randperm(length(all_y)));
            y_pred_null = nan(size(all_y));

            for fold = 1:cvp.NumTestSets
                train_idx = training(cvp, fold);
                test_idx = test(cvp, fold);
                p_cv_null = polyfit(all_x(train_idx), shuffled_y(train_idx), 1);
                y_pred_null(test_idx) = polyval(p_cv_null, all_x(test_idx));
            end

            SS_res_null = sum((shuffled_y - y_pred_null).^2);
            SS_tot_null = sum((shuffled_y - mean(shuffled_y)).^2);
            null_r2_dist(p_idx) = 1 - (SS_res_null / SS_tot_null);
        end

        % p-value and null CI bounds
        cv_p_vals(i) = sum(null_r2_dist >= true_r2) / n_perms;
        null_r2_lower_999(i) = prctile(null_r2_dist, 0.1);
        null_r2_upper_999(i) = prctile(null_r2_dist, 99.9);
        null_r2_upper_99(i)  = prctile(null_r2_dist, 99);
        null_r2_upper_95(i)  = prctile(null_r2_dist, 95);
    end
end

% Sort markers by decreasing Cross-Validated R^2 (used by both figures below)
[~, sort_idx] = sort(cv_r2_vals, 'descend', 'MissingPlacement', 'last');

% FIGURE 1: Scatter Plots of Baseline MADRS vs Baseline Marker Value
figure('Position', [100 100 1400 900], 'Color', '#DAF2FB');
t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:n_markers
    i = sort_idx(k);
    v = ~isnan(D0_z{:, i}) & ~isnan(madrs0);
    x = D0_z{v, i}; y = madrs0(v);

    nexttile; hold on;
    xline(0, 'k');
    scatter(x, y, 8, 'r', 'filled');
    if numel(x) > 2
        p_line = polyfit(x, y, 1);
        xl = xlim();
        y_fit = polyval(p_line, xl);
        plot(xl, y_fit, '-r', 'LineWidth', 1);
        xlim(xl);
    end

    % Significant if the CV R^2 beats chance (no directional hypothesis here,
    % unlike Section 13, so there's no slope-sign condition)
    if cv_r2_vals(i) > 0 && cv_p_vals(i) < 0.05
        t_color = 'k';
    else
        t_color = [0.6 0.6 0.6]; % Gray
    end
    title(ind_markers{i}, 'Interpreter', 'none', 'FontSize', 10, 'Color', t_color);
    box off;

    if k == 1
        set(gca, 'XTick', 0, 'YTick', [0 40], 'Color', 'w');
    else
        set(gca, 'XTick', 0, 'YTick', [], 'Color', 'w');
    end
end
xlabel(t, 'baseline marker (z)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(t, 'baseline MADRS', 'FontSize', 12, 'FontWeight', 'bold');
title(t, 'Supp fig: Baseline MADRS vs raw baseline marker values');

% FIGURE 2: Significance Summary (True R^2 vs Null 95/99/99.9% CIs)
figure('Position', [150 150 1200 450], 'Color', '#DAF2FB');
hold on;
title('Supp fig: Baseline raw marker \rightarrow baseline MADRS: CV R^2 significance');
x_ax = 1:n_markers;
y_true_perc = cv_r2_vals(sort_idx) * 100; % Convert to percentage

% Create separate semi-transparent rectangles for each marker's Null CIs
rect_w = 0.6; % Width of the rectangles
for k = 1:n_markers
    orig_idx = sort_idx(k);
    draw_ci_bands(k, rect_w/2, null_r2_lower_999(orig_idx) * 100, ...
        null_r2_upper_999(orig_idx) * 100, null_r2_upper_99(orig_idx) * 100, null_r2_upper_95(orig_idx) * 100);
end

% True cross-validated R^2 (black, connected)
plot(x_ax, y_true_perc, '-ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

% Formatting
yline(0, 'k-'); % Baseline of zero variance explained

% Generate dynamic tick labels using TeX to color code significance
custom_labels = cell(1, n_markers);
for k = 1:n_markers
    orig_idx = sort_idx(k);
    % Escape underscores for TeX interpreter
    safe_name = strrep(ind_markers{orig_idx}, '_', '\_');

    if cv_r2_vals(orig_idx) > 0 && cv_p_vals(orig_idx) < 0.05
        custom_labels{k} = sprintf('\\color{black}%s', safe_name);
    else
        custom_labels{k} = sprintf('\\color{gray}%s', safe_name);
    end
end

set(gca, 'XTick', x_ax, 'XTickLabel', custom_labels, 'TickLabelInterpreter', 'tex', ...
    'Color', 'w', 'FontSize', 11, 'TickDir', 'out', 'xlim', [0 n_markers+1],'ylim',[-50 40],'ytick',[-50 0 10]);
xtickangle(90);
ytickformat('%g%%'); % Formats the Y-axis ticks with a percentage sign
ylabel('CV R^2');
box off;

% Per marker: same as above, but correlating baseline MADRS with the
% univariate posterior P(H|marker) instead of the raw marker value. The
% posterior comes from a single-marker Gaussian LLR comparing the
% healthy-baseline distribution N(0,1) (true by construction of the
% z-score) against the depressed-baseline distribution N(md,sd) fit to
% D0_z -- both taken alone, with no higher-dimensional / joint-marker
% information.
% The direction of the hypothesis (higher P(H) -> lower MADRS) is built into
% the model itself via fit_neg_slope: a fold whose fit comes out positive
% falls back to the training mean. The CV R^2 then measures direction and
% predictive strength at once, so no separate slope-sign condition is needed.
cv_r2_vals_post        = nan(1, n_markers);
cv_p_vals_post         = nan(1, n_markers);
null_r2_lower_999_post = nan(1, n_markers);
null_r2_upper_999_post = nan(1, n_markers);
null_r2_upper_99_post  = nan(1, n_markers);
null_r2_upper_95_post  = nan(1, n_markers);

for i = 1:n_markers
    v = ~isnan(D0_z{:, i}) & ~isnan(madrs0);
    all_z = D0_z{v, i};
    all_y = madrs0(v);

    % Require enough points to do LOOCV safely
    if length(all_z) >= 5
        % The LLR's depressed-baseline mean/SD are refit inside each fold, on
        % the training subjects only, so a held-out subject's own marker value
        % never contributes to the P(H) transform applied to it. This depends
        % only on the marker values, not on MADRS, so the same per-fold P(H)
        % (all_x_train/all_x_test) is reused for both the true and the
        % permuted-MADRS null below.
        cvp = cvpartition(length(all_y), 'LeaveOut');
        all_x_train = cell(1, cvp.NumTestSets);
        all_x_test  = cell(1, cvp.NumTestSets);
        for fold = 1:cvp.NumTestSets
            train_idx = training(cvp, fold);
            test_idx = test(cvp, fold);
            md_fold = mean(all_z(train_idx), 'omitnan');
            sd_fold = std(all_z(train_idx), 'omitnan');
            LLR_fold = @(z) log(normpdf(z, 0, 1)) - log(normpdf(z, md_fold, sd_fold));
            all_x_train{fold} = sigmoid(LLR_fold(all_z(train_idx)));
            all_x_test{fold}  = sigmoid(LLR_fold(all_z(test_idx)));
        end

        % LOOCV R^2 from out-of-fold predictions
        y_pred_true = nan(size(all_y));
        for fold = 1:cvp.NumTestSets
            train_idx = training(cvp, fold);
            test_idx = test(cvp, fold);
            p_cv = fit_neg_slope(all_x_train{fold}, all_y(train_idx));
            y_pred_true(test_idx) = polyval(p_cv, all_x_test{fold});
        end
        SS_res_true = sum((all_y - y_pred_true).^2);
        SS_tot_true = sum((all_y - mean(all_y)).^2);
        true_r2 = 1 - (SS_res_true / SS_tot_true);
        cv_r2_vals_post(i) = true_r2;

        % Permutation null: shuffle baseline MADRS and recompute LOOCV R^2
        null_r2_dist = nan(1, n_perms);
        for p_idx = 1:n_perms
            shuffled_y = all_y(randperm(length(all_y)));
            y_pred_null = nan(size(all_y));

            for fold = 1:cvp.NumTestSets
                train_idx = training(cvp, fold);
                test_idx = test(cvp, fold);
                p_cv_null = fit_neg_slope(all_x_train{fold}, shuffled_y(train_idx));
                y_pred_null(test_idx) = polyval(p_cv_null, all_x_test{fold});
            end

            SS_res_null = sum((shuffled_y - y_pred_null).^2);
            SS_tot_null = sum((shuffled_y - mean(shuffled_y)).^2);
            null_r2_dist(p_idx) = 1 - (SS_res_null / SS_tot_null);
        end

        % p-value and null CI bounds
        cv_p_vals_post(i) = sum(null_r2_dist >= true_r2) / n_perms;
        null_r2_lower_999_post(i) = prctile(null_r2_dist, 0.1);
        null_r2_upper_999_post(i) = prctile(null_r2_dist, 99.9);
        null_r2_upper_99_post(i)  = prctile(null_r2_dist, 99);
        null_r2_upper_95_post(i)  = prctile(null_r2_dist, 95);
    end
end

% Sort markers by decreasing direction-constrained CV R^2 (used by both figures below)
[~, sort_idx_post] = sort(cv_r2_vals_post, 'descend', 'MissingPlacement', 'last');

% FIGURE 3: Scatter Plots of Baseline MADRS vs Baseline Posterior P(H|marker)
figure('Position', [100 100 1400 900], 'Color', '#DAF2FB');
t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:n_markers
    i = sort_idx_post(k);

    md = mean(D0_z{:, i}, 'omitnan');
    sd = std(D0_z{:, i}, 'omitnan');
    LLR = @(z) log(normpdf(z, 0, 1)) - log(normpdf(z, md, sd));

    v = ~isnan(D0_z{:, i}) & ~isnan(madrs0);
    x = sigmoid(LLR(D0_z{v, i})); y = madrs0(v);

    nexttile; hold on;
    xline(0.5, 'k');
    scatter(x, y, 8, 'r', 'filled');
    if numel(x) > 2
        p_line = polyfit(x, y, 1);
        xl = xlim();
        y_fit = polyval(p_line, xl);
        plot(xl, y_fit, '-r', 'LineWidth', 1);
        xlim(xl);
    end

    % The negative-slope requirement is already inside the CV R^2, so beating
    % the permutation null is the whole significance criterion
    if cv_p_vals_post(i) < 0.05
        t_color = 'k';
    else
        t_color = [0.6 0.6 0.6]; % Gray
    end
    title(ind_markers{i}, 'Interpreter', 'none', 'FontSize', 10, 'Color', t_color);
    box off;

    if k == 1
        set(gca, 'XTick', [0 1], 'YTick', [0 40], 'Color', 'w');
    else
        set(gca, 'XTick', [0 1], 'YTick', [], 'Color', 'w');
    end
end
xlabel(t, 'P(H|marker)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(t, 'baseline MADRS', 'FontSize', 12, 'FontWeight', 'bold');
title(t, 'Baseline MADRS vs baseline posterior P(H|marker)');

% FIGURE 4: Significance Summary (True R^2 vs Null 95/99/99.9% CIs)
figure('Position', [150 150 1200 450], 'Color', '#DAF2FB');
hold on;
title('Baseline P(H|marker) \rightarrow baseline MADRS: CV R^2 significance');
x_ax = 1:n_markers;
y_true_perc = cv_r2_vals_post(sort_idx_post) * 100; % Convert to percentage

% Create separate semi-transparent rectangles for each marker's Null CIs
rect_w = 0.6; % Width of the rectangles
for k = 1:n_markers
    orig_idx = sort_idx_post(k);
    draw_ci_bands(k, rect_w/2, null_r2_lower_999_post(orig_idx) * 100, ...
        null_r2_upper_999_post(orig_idx) * 100, null_r2_upper_99_post(orig_idx) * 100, null_r2_upper_95_post(orig_idx) * 100);
end

% True cross-validated R^2 (black, connected)
plot(x_ax, y_true_perc, '-ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

% Formatting
yline(0, 'k-'); % Baseline of zero variance explained

% Generate dynamic tick labels using TeX to color code significance
custom_labels = cell(1, n_markers);
for k = 1:n_markers
    orig_idx = sort_idx_post(k);
    % Escape underscores for TeX interpreter
    safe_name = strrep(ind_markers{orig_idx}, '_', '\_');

    if cv_p_vals_post(orig_idx) < 0.05
        custom_labels{k} = sprintf('\\color{black}%s', safe_name);
    else
        custom_labels{k} = sprintf('\\color{gray}%s', safe_name);
    end
end

set(gca, 'XTick', x_ax, 'XTickLabel', custom_labels, 'TickLabelInterpreter', 'tex', ...
    'Color', 'w', 'FontSize', 11, 'TickDir', 'out', 'xlim', [0 n_markers+1],'ylim',[-50 40],'ytick',[-50 0 10]);
xtickangle(90);
ytickformat('%g%%'); % Formats the Y-axis ticks with a percentage sign
ylabel('CV R^2');
box off;

%% 12. Greedy marker panel -> baseline MADRS (joint posterior P(H|panel))
% Section 11 treats each marker on its own; here markers are combined into a
% panel and MADRS is predicted from a single joint posterior P(H|panel),
% built the same way as Section 11 Figs 3/4 but from a multivariate Gaussian
% fit to the panel's JOINT distribution, rather than one univariate Gaussian
% per marker averaged together. A forward greedy search grows the panel one
% marker at a time, adding at each step whichever remaining marker most
% improves the panel's leave-one-out CV R^2 -- so the search is driven by
% out-of-sample performance, never by the in-sample fit, which would simply
% add every marker.
%
% Healthy is fit as a multivariate Gaussian from all H0 subjects -- they never
% enter the MADRS regression's CV folds, so this carries no leakage. Depressed
% is fit as a multivariate Gaussian refit inside every CV fold from that
% fold's training subjects only, exactly as in the univariate case above, so a
% held-out subject's own marker values never contribute to the joint Gaussian
% used to transform it.
%
% Each panel keeps whichever subjects have all of *that panel's* markers
% present, exactly as Section 11 keeps a marker's own available subjects --
% not the subjects with every one of the 25 markers present. So a singleton
% panel is scored on the same subjects (and the same LOOCV folds) as its
% Section 11 counterpart, and only larger panels lose subjects to missingness.
v0 = ~isnan(madrs0);
Z_greedy = D0_z{v0, :};
H0Z_all  = H0_z{:, :};
y_greedy = madrs0(v0);
n_allmarkers = sum(~any(isnan(Z_greedy), 2));

% Panels much larger than n/3 predictors give an unstable depressed-covariance
% estimate within each training fold (and an exactly singular one once the
% panel reaches n-1); n_allmarkers is the worst case (all markers required at
% once), so capping against it keeps every panel size safe.
max_panel = min(n_markers, floor(n_allmarkers / 3));
fprintf('Greedy MADRS regression: %d MADRS cases (%d once all %d markers are required), panels up to %d markers.\n', ...
    numel(y_greedy), n_allmarkers, n_markers, max_panel);

panel_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Real data: the order markers entered, and the panel CV R^2 after each addition
[greedy_order, cum_r2_greedy] = greedy_forward_r2_joint(Z_greedy, H0Z_all, y_greedy, max_panel, panel_cache);
greedy_names = ind_markers(greedy_order);

% Permutation null: the panel order is fixed from the real data above (as in
% Section 2/4's stable panel), so the null for each cumulative panel size k
% shuffles MADRS through that *same* fixed k-marker panel -- exactly Section
% 11's per-marker null, just applied to a growing panel instead of one marker.
% It does not re-run marker selection under permutation, so it does not
% correct for the optimism of having picked that panel by greedy search.
% n_perms is larger than Section 11's (1000): each permutation here is just a
% cached-transform linear-regression refit, not a Gaussian refit, so it's
% cheap -- and the 99.9th-percentile band needs many samples in the tail to
% stop jumping around from run to run.
n_perms = 20000;
null_curves_greedy = nan(n_perms, max_panel);
for p_idx = 1:n_perms
    if mod(p_idx, 2000) == 0
        fprintf('  greedy permutation %d / %d...\n', p_idx, n_perms);
    end
    y_perm = y_greedy(randperm(numel(y_greedy)));
    for k = 1:max_panel
        xt = panel_cache(mat2str(sort(greedy_order(1:k))));
        null_curves_greedy(p_idx, k) = score_from_transform(xt, y_perm);
    end
end

% Per panel size: 1-tailed p-value and the null CI bounds (percentile floor, as
% in Section 11, rather than the raw minimum -- CV R^2 has no lower limit)
cum_p_greedy          = nan(1, max_panel);
null_r2_lower_999_gr  = nan(1, max_panel);
null_r2_upper_999_gr  = nan(1, max_panel);
null_r2_upper_99_gr   = nan(1, max_panel);
null_r2_upper_95_gr   = nan(1, max_panel);
for k = 1:max_panel
    cum_p_greedy(k)         = sum(null_curves_greedy(:, k) >= cum_r2_greedy(k)) / n_perms;
    null_r2_lower_999_gr(k) = prctile(null_curves_greedy(:, k), 0.1);
    null_r2_upper_999_gr(k) = prctile(null_curves_greedy(:, k), 99.9);
    null_r2_upper_99_gr(k)  = prctile(null_curves_greedy(:, k), 99);
    null_r2_upper_95_gr(k)  = prctile(null_curves_greedy(:, k), 95);
end

[best_r2, best_k] = max(cum_r2_greedy);
fprintf('Best panel: %d marker(s), CV R^2 = %.3f, p = %.3f (%s)\n', ...
    best_k, best_r2, cum_p_greedy(best_k), strjoin(greedy_names(1:best_k), ', '));

% FIGURE: cumulative panel CV R^2 against the permutation null bands. Markers
% are filled where the panel beats the 95% null, hollow otherwise.
figure('Position', [150 150 1200 450], 'Color', '#DAF2FB');
hold on;
x_ax = 1:max_panel;
rect_w = 0.6;
for k = 1:max_panel
    draw_ci_bands(k, rect_w/2, null_r2_lower_999_gr(k) * 100, ...
        null_r2_upper_999_gr(k) * 100, null_r2_upper_99_gr(k) * 100, null_r2_upper_95_gr(k) * 100);
end
yline(0, 'k-');
sig = cum_r2_greedy > null_r2_upper_95_gr;
plot(x_ax, cum_r2_greedy * 100, '-k', 'LineWidth', 1);
plot(x_ax(sig),  cum_r2_greedy(sig) * 100,  'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
plot(x_ax(~sig), cum_r2_greedy(~sig) * 100, 'ok', 'MarkerFaceColor', 'w', 'MarkerSize', 5);
ylabel('panel CV R^2');
ytickformat('%g%%');
set(gca, 'XTick', x_ax, 'XTickLabel', strrep(greedy_names, '_', '\_'), 'TickLabelInterpreter', 'tex', ...
    'TickDir', 'out', 'FontSize', 12, 'xlim', [0 max_panel+1], 'ylim', [-10 15], 'Color', 'w');
xtickangle(90);
box off;
title('Greedy marker panel \rightarrow baseline MADRS');

%% 13. Change in MADRS vs change in P(H)
% Align baseline with post-treatment per arm, then take the MADRS change
[~, bA_D, aA_D] = intersect(D0_full.id, Da_full.id, 'stable');
[~, bP_D, aP_D] = intersect(D0_full.id, Dp_full.id, 'stable');
dMADRS_A = Da_full.MADRS(aA_D) - D0_full.MADRS(bA_D);
dMADRS_P = Dp_full.MADRS(aP_D) - D0_full.MADRS(bP_D);

% Per marker: cross-validated R^2 of delta-P(H) -> delta-MADRS and a
% permutation null distribution of R^2. The direction of the hypothesis (more
% restoration -> more MADRS improvement) is built into the model itself via
% fit_neg_slope: a fold whose fit comes out positive falls back to the
% training mean. The CV R^2 then measures direction and predictive strength at
% once, so no separate slope-sign condition is needed.
cv_r2_vals        = nan(1, n_markers);
cv_p_vals         = nan(1, n_markers);
null_r2_lower_999 = nan(1, n_markers);  % common floor for the CI bands
null_r2_upper_999 = nan(1, n_markers);
null_r2_upper_99  = nan(1, n_markers);
null_r2_upper_95  = nan(1, n_markers);
n_perms           = 1000;

for i = 1:n_markers
    % Per-marker LLR boundary from the baseline depressed stats (vs N(0,1) healthy)
    md = mean(D0_z{:, i}, 'omitnan');
    sd = std(D0_z{:, i}, 'omitnan');
    LLR = @(z) log(normpdf(z, 0, 1)) - log(normpdf(z, md, sd));

    % Change in P(H) from before to after (Ayahuasca arm only)
    dA = sigmoid(LLR(Da_z{aA_D, i})) - sigmoid(LLR(D0_z{bA_D, i}));

    vA = ~isnan(dA) & ~isnan(dMADRS_A);

    all_x = dA(vA);
    all_y = dMADRS_A(vA);

    % Require enough points to do LOOCV safely
    if length(all_x) >= 5
        % LOOCV R^2 from out-of-fold predictions
        cvp = cvpartition(length(all_y), 'LeaveOut');
        y_pred_true = nan(size(all_y));
        for fold = 1:cvp.NumTestSets
            train_idx = training(cvp, fold);
            test_idx = test(cvp, fold);
            p_cv = fit_neg_slope(all_x(train_idx), all_y(train_idx));
            y_pred_true(test_idx) = polyval(p_cv, all_x(test_idx));
        end
        SS_res_true = sum((all_y - y_pred_true).^2);
        SS_tot_true = sum((all_y - mean(all_y)).^2);
        true_r2 = 1 - (SS_res_true / SS_tot_true);
        cv_r2_vals(i) = true_r2;

        % Permutation null: shuffle delta-MADRS and recompute LOOCV R^2
        null_r2_dist = nan(1, n_perms);
        for p_idx = 1:n_perms
            shuffled_y = all_y(randperm(length(all_y)));
            y_pred_null = nan(size(all_y));

            for fold = 1:cvp.NumTestSets
                train_idx = training(cvp, fold);
                test_idx = test(cvp, fold);
                p_cv_null = fit_neg_slope(all_x(train_idx), shuffled_y(train_idx));
                y_pred_null(test_idx) = polyval(p_cv_null, all_x(test_idx));
            end

            SS_res_null = sum((shuffled_y - y_pred_null).^2);
            SS_tot_null = sum((shuffled_y - mean(shuffled_y)).^2);
            null_r2_dist(p_idx) = 1 - (SS_res_null / SS_tot_null);
        end

        % p-value and null CI bounds
        cv_p_vals(i) = sum(null_r2_dist >= true_r2) / n_perms;
        null_r2_lower_999(i) = prctile(null_r2_dist, 0.1);
        null_r2_upper_999(i) = prctile(null_r2_dist, 99.9);
        null_r2_upper_99(i)  = prctile(null_r2_dist, 99);
        null_r2_upper_95(i)  = prctile(null_r2_dist, 95);
    end
end

% Sort markers by decreasing direction-constrained CV R^2
[~, sort_idx] = sort(cv_r2_vals, 'descend', 'MissingPlacement', 'last');

% FIGURE 1: Scatter Plots of Delta MADRS vs Delta P(H)
figure('Position', [100 100 1400 900], 'Color', '#DAF2FB');
t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:n_markers
    i = sort_idx(k);

    md = mean(D0_z{:, i}, 'omitnan');
    sd = std(D0_z{:, i}, 'omitnan');
    LLR = @(z) log(normpdf(z, 0, 1)) - log(normpdf(z, md, sd));

    % Calculate delta P(H) (Ayahuasca arm only)
    dA = sigmoid(LLR(Da_z{aA_D, i})) - sigmoid(LLR(D0_z{bA_D, i}));

    vA = ~isnan(dA) & ~isnan(dMADRS_A);

    plot_dA = dA(vA); plot_dMA = dMADRS_A(vA);

    nexttile; hold on;

    xline(0, 'k');
    yline(0, 'k');

    scatter(plot_dA, plot_dMA, 8, 'r', 'filled');

    all_x = plot_dA;
    all_y = plot_dMA;
    if length(all_x) > 2
        p_line = polyfit(all_x, all_y, 1);
        xl = xlim();
        y_fit = polyval(p_line, xl);
        plot(xl, y_fit, '-r', 'LineWidth', 1);
        xlim(xl);
    end

    % The negative-slope requirement is already inside the CV R^2, so beating
    % the permutation null is the whole significance criterion
    if cv_p_vals(i) < 0.05
        t_color = 'k';
    else
        t_color = [0.6 0.6 0.6]; % Gray
    end
    title(ind_markers{i}, 'Interpreter', 'none', 'FontSize', 10, 'Color', t_color);

    box off;

    if k == 1
        set(gca, 'XTick', 0, 'YTick', [-40 0], 'Color', 'w');
    else
        set(gca, 'XTick', 0, 'YTick', [], 'Color', 'w');
    end
end
xlabel(t, '\Delta P(H)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(t, '\Delta MADRS', 'FontSize', 12, 'FontWeight', 'bold');
title(t, '\Delta MADRS vs \Delta P(H) per marker (Ayahuasca)');

% FIGURE 2: Significance Summary (True R^2 vs Null 95/99/99.9% CIs)
figure('Position', [150 150 1200 450], 'Color', '#DAF2FB');
hold on;
title('\Delta P(H) \rightarrow \Delta MADRS: CV R^2 significance');
x_ax = 1:n_markers;
y_true_perc = cv_r2_vals(sort_idx) * 100; % Convert to percentage

% Create separate semi-transparent rectangles for each marker's Null CIs
rect_w = 0.6; % Width of the rectangles
for k = 1:n_markers
    orig_idx = sort_idx(k);
    draw_ci_bands(k, rect_w/2, null_r2_lower_999(orig_idx) * 100, ...
        null_r2_upper_999(orig_idx) * 100, null_r2_upper_99(orig_idx) * 100, null_r2_upper_95(orig_idx) * 100);
end

% True cross-validated R^2 (black, connected)
plot(x_ax, y_true_perc, '-ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

% Formatting
yline(0, 'k-'); % Baseline of zero variance explained

% Generate dynamic tick labels using TeX to color code significance
custom_labels = cell(1, n_markers);
for k = 1:n_markers
    orig_idx = sort_idx(k);
    % Escape underscores for TeX interpreter
    safe_name = strrep(ind_markers{orig_idx}, '_', '\_');

    % Apply the compound significance condition for the axis label color
    if cv_p_vals(orig_idx) < 0.05
        custom_labels{k} = sprintf('\\color{black}%s', safe_name);
    else
        custom_labels{k} = sprintf('\\color{gray}%s', safe_name);
    end
end

set(gca, 'XTick', x_ax, 'XTickLabel', custom_labels, 'TickLabelInterpreter', 'tex', ...
    'Color', 'w', 'FontSize', 11, 'TickDir', 'out', 'xlim', [0 n_markers+1],'ylim',[-50 40],'ytick',[-50 0 10]);
xtickangle(90);
ytickformat('%g%%'); % Formats the Y-axis ticks with a percentage sign
ylabel('CV R^2');
box off;

%% 13.5 Nested LOOCV elastic-net: predict baseline MADRS from all marker baseline P(H)
% Same question as Section 12's greedy panel -- does baseline MADRS relate to
% where subjects sit on the healthy-depressed spectrum -- but letting a
% regression freely weight all n_markers P(H) values at once, rather than
% growing a panel by greedy forward selection. Predictors are each marker's
% P(H|marker) (its H-D posterior), not the raw marker value.
%
% With n_markers predictors and only a few dozen depressed subjects, ordinary
% least squares would have as many (or more) free parameters as data points:
% it would fit the training data perfectly and explain nothing out of sample.
% Elastic net adds a combined L1+L2 penalty on the coefficients, shrinking
% them all and zeroing out the least useful ones, so the model stays
% identifiable and testable at this sample size. (Alpha = 0.5 below splits
% the penalty evenly between L1 (lasso, sparsity) and L2 (ridge, stability
% under correlated markers); lambda -- how hard to shrink -- is tuned by the
% inner CV rather than fixed by hand.)
%
% Healthy is fixed from all H0 subjects (they never enter these folds, so no
% leakage risk). Depressed is refit inside every outer LOOCV fold from that
% fold's training subjects only, exactly as in Sections 11/12, so a held-out
% subject's own marker value never contributes to the P(H) transform used to
% score them. The inner 5-fold CV only tunes the elastic-net lambda, so it
% reuses that already-leakage-free training-fold feature matrix.
v_b = ~any(isnan(D0_z{:,:}), 2) & ~isnan(madrs0);
Z_b = D0_z{v_b, :};
y_b = madrs0(v_b);
n_b = numel(y_b);

mu_H_b = mean(H0_z{:,:}, 1, 'omitnan');
sd_H_b = std(H0_z{:,:}, 0, 1, 'omitnan');

fprintf('Running Nested Elastic Net CV (baseline MADRS) on %d complete cases...\n', n_b);

y_pred_enet_b = nan(n_b, 1);
coeff_matrix_b = zeros(n_markers, n_b);

% Outer LOOCV loop; the inner 5-fold CV (below) tunes lambda
for fold = 1:n_b
    train_idx = setdiff(1:n_b, fold);
    test_idx = fold;

    % Per-fold, per-marker P(H): depressed refit from training subjects only
    mu_D_fold = mean(Z_b(train_idx, :), 1);
    sd_D_fold = std(Z_b(train_idx, :), 0, 1);
    LLR_fold = @(Z) log(normpdf(Z, mu_H_b, sd_H_b)) - log(normpdf(Z, mu_D_fold, sd_D_fold));

    X_train = sigmoid(LLR_fold(Z_b(train_idx, :)));
    X_test  = sigmoid(LLR_fold(Z_b(test_idx, :)));
    Y_train = y_b(train_idx);

    % Freeze seed ONLY for the inner CV folds so the lambda tuning is deterministic
    rng(fold);

    % Inner CV: Tune lambda using 5-fold CV. Alpha = 0.5 specifies Elastic Net.
    [B, FitInfo] = lasso(X_train, Y_train, 'CV', 5, 'Alpha', 0.5, 'Standardize', true);

    % Select the lambda that minimizes Mean Squared Error on the inner folds
    idx_opt = FitInfo.IndexMinMSE;
    opt_B = B(:, idx_opt);
    opt_Intercept = FitInfo.Intercept(idx_opt);

    coeff_matrix_b(:, fold) = opt_B;
    y_pred_enet_b(test_idx) = X_test * opt_B + opt_Intercept;
end

% Final multivariate test R^2
SS_res_b = sum((y_b - y_pred_enet_b).^2);
SS_tot_b = sum((y_b - mean(y_b)).^2);
r2_cv_enet_b = 1 - SS_res_b / SS_tot_b;

fprintf('Nested LOOCV Elastic Net (baseline) Test R^2: %.3f\n', r2_cv_enet_b);

% --- Plotting the Selected Features ---

% Calculate how often each marker was selected (coefficient != 0)
selection_freq_b = sum(coeff_matrix_b ~= 0, 2) / n_b * 100;

% Sort markers by selection frequency (most selected first)
[sorted_freq_b, sort_idx_b] = sort(selection_freq_b, 'descend');
sorted_markers_b = ind_markers(sort_idx_b);
sorted_coeffs_b = coeff_matrix_b(sort_idx_b, :);

figure('Position', [150 150 1200 800], 'Color', '#DAF2FB');
t_enet_b = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_enet_b, 'Elastic-net: baseline MADRS from all marker baseline P(H)');

% Subplot 1: Selection Frequency
ax1 = nexttile; hold on;
bar(1:n_markers, sorted_freq_b, 'FaceColor', [0.4 0.6 0.8], 'EdgeColor', 'k');
yline(50, 'k--', 'LineWidth', 1); % 50% reference line
set(gca, 'XTick', 1:n_markers, 'XTickLabel', [], 'Color', 'w', ...
    'TickDir', 'out', 'FontSize', 12, 'xlim', [0 n_markers+1], 'ylim', [0 105]);
ylabel('Selection Frequency (%)');
title(sprintf('Elastic Net Feature Selection (Nested CV R^2 = %.2f)', r2_cv_enet_b));
box off;

% Subplot 2: Coefficient Weights Distribution
ax2 = nexttile; hold on;
yline(0, 'k-', 'LineWidth', 1);
for k = 1:n_markers
    weights = sorted_coeffs_b(k, sorted_coeffs_b(k,:) ~= 0);
    if ~isempty(weights)
        x_jitter = k + (rand(1, length(weights)) - 0.5) * 0.2;
        scatter(x_jitter, weights, 15, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
        plot([k-0.3, k+0.3], [mean(weights), mean(weights)], 'r-', 'LineWidth', 2);
    end
end
set(gca, 'XTick', 1:n_markers, 'XTickLabel', sorted_markers_b, 'TickLabelInterpreter', 'none', ...
    'Color', 'w', 'TickDir', 'out', 'FontSize', 11, 'xlim', [0 n_markers+1]);
xtickangle(90);
ylabel('Coefficient Weight');
box off;

%% 14. Change in MADRS vs raw marker shift
% Same subject alignment as above; predictor is the raw Z-scored change
% (post minus pre) in each marker rather than the derived delta-P(H).

cv_r2_raw             = nan(1, n_markers);
cv_p_raw              = nan(1, n_markers);
null_r2_lower_999_raw = nan(1, n_markers);
null_r2_upper_999_raw = nan(1, n_markers);
null_r2_upper_99_raw  = nan(1, n_markers);
null_r2_upper_95_raw  = nan(1, n_markers);

for i = 1:n_markers
    dA_raw = Da_z{aA_D, i} - D0_z{bA_D, i};

    vA = ~isnan(dA_raw) & ~isnan(dMADRS_A);

    all_x = dA_raw(vA);
    all_y = dMADRS_A(vA);

    if length(all_x) >= 5
        cvp = cvpartition(length(all_y), 'LeaveOut');
        y_pred_true = nan(size(all_y));
        for fold = 1:cvp.NumTestSets
            train_idx = training(cvp, fold);
            test_idx  = test(cvp, fold);
            p_cv = polyfit(all_x(train_idx), all_y(train_idx), 1);
            y_pred_true(test_idx) = polyval(p_cv, all_x(test_idx));
        end
        SS_res_true = sum((all_y - y_pred_true).^2);
        SS_tot_true = sum((all_y - mean(all_y)).^2);
        true_r2 = 1 - SS_res_true / SS_tot_true;
        cv_r2_raw(i) = true_r2;

        null_r2_dist = nan(1, n_perms);
        for p_idx = 1:n_perms
            shuffled_y = all_y(randperm(length(all_y)));
            y_pred_null = nan(size(all_y));
            for fold = 1:cvp.NumTestSets
                train_idx = training(cvp, fold);
                test_idx  = test(cvp, fold);
                p_cv_null = polyfit(all_x(train_idx), shuffled_y(train_idx), 1);
                y_pred_null(test_idx) = polyval(p_cv_null, all_x(test_idx));
            end
            SS_res_null = sum((shuffled_y - y_pred_null).^2);
            SS_tot_null = sum((shuffled_y - mean(shuffled_y)).^2);
            null_r2_dist(p_idx) = 1 - SS_res_null / SS_tot_null;
        end

        cv_p_raw(i)              = sum(null_r2_dist >= true_r2) / n_perms;
        null_r2_lower_999_raw(i) = prctile(null_r2_dist, 0.1);
        null_r2_upper_999_raw(i) = prctile(null_r2_dist, 99.9);
        null_r2_upper_99_raw(i)  = prctile(null_r2_dist, 99);
        null_r2_upper_95_raw(i)  = prctile(null_r2_dist, 95);
    end
end

[~, sort_idx_raw] = sort(cv_r2_raw, 'descend', 'MissingPlacement', 'last');

% FIGURE 1: Scatter Plots of Delta MADRS vs Raw Marker Delta
figure('Position', [100 100 1400 900], 'Color', '#DAF2FB');
t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:n_markers
    i = sort_idx_raw(k);

    dA_raw = Da_z{aA_D, i} - D0_z{bA_D, i};

    vA = ~isnan(dA_raw) & ~isnan(dMADRS_A);

    plot_dA = dA_raw(vA); plot_dMA = dMADRS_A(vA);

    nexttile; hold on;
    xline(0, 'k');
    yline(0, 'k');
    scatter(plot_dA, plot_dMA, 8, 'r', 'filled');

    all_x = plot_dA;
    all_y = plot_dMA;
    if length(all_x) > 2
        p_line = polyfit(all_x, all_y, 1);
        xl = xlim();
        y_fit = polyval(p_line, xl);
        plot(xl, y_fit, '-r', 'LineWidth', 1);
        xlim(xl);
    end

    if cv_r2_raw(i) > 0 && cv_p_raw(i) < 0.05
        t_color = 'k';
    else
        t_color = [0.6 0.6 0.6];
    end
    title(ind_markers{i}, 'Interpreter', 'none', 'FontSize', 10, 'Color', t_color);
    box off;
    if k == 1
        set(gca, 'XTick', 0, 'YTick', [-40 0], 'Color', 'w');
    else
        set(gca, 'XTick', 0, 'YTick', [], 'Color', 'w');
    end
end
xlabel(t, '\Delta Z-score', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(t, '\Delta MADRS', 'FontSize', 12, 'FontWeight', 'bold');
title(t, '\Delta MADRS vs raw marker shift per marker (Ayahuasca)');

% FIGURE 2: Significance Summary (True R^2 vs Null 95/99/99.9% CIs)
figure('Position', [150 150 1200 450], 'Color', '#DAF2FB');
hold on;
title('Raw shift \rightarrow \Delta MADRS: CV R^2 significance');
x_ax = 1:n_markers;
y_true_perc_raw = cv_r2_raw(sort_idx_raw) * 100;

rect_w = 0.6;
for k = 1:n_markers
    orig_idx = sort_idx_raw(k);
    draw_ci_bands(k, rect_w/2, null_r2_lower_999_raw(orig_idx) * 100, ...
        null_r2_upper_999_raw(orig_idx) * 100, null_r2_upper_99_raw(orig_idx) * 100, ...
        null_r2_upper_95_raw(orig_idx) * 100);
end

plot(x_ax, y_true_perc_raw, '-ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
yline(0, 'k-');

custom_labels = cell(1, n_markers);
for k = 1:n_markers
    orig_idx = sort_idx_raw(k);
    safe_name = strrep(ind_markers{orig_idx}, '_', '\_');
    if cv_r2_raw(orig_idx) > 0 && cv_p_raw(orig_idx) < 0.05
        custom_labels{k} = sprintf('\\color{black}%s', safe_name);
    else
        custom_labels{k} = sprintf('\\color{gray}%s', safe_name);
    end
end

set(gca, 'XTick', x_ax, 'XTickLabel', custom_labels, 'TickLabelInterpreter', 'tex', ...
    'Color', 'w', 'FontSize', 11, 'TickDir', 'out', 'xlim', [0 n_markers+1], 'ylim', [-50 40], 'ytick', [-50 0 10]);
xtickangle(90);
ytickformat('%g%%');
ylabel('CV R^2');
box off;

%% 15. Nested LOOCV elastic-net: predict delta-MADRS from all marker delta-P(H)
% Build the feature matrix X (per-marker delta-P(H)) and target Y (delta-MADRS)
all_dA = nan(length(aA_D), n_markers);
all_dP = nan(length(aP_D), n_markers);

for i = 1:n_markers
    md = mean(D0_z{:, i}, 'omitnan');
    sd = std(D0_z{:, i}, 'omitnan');
    LLR = @(z) log(normpdf(z, 0, 1)) - log(normpdf(z, md, sd));

    all_dA(:, i) = LLR(Da_z{aA_D, i}) - LLR(D0_z{bA_D, i});
    all_dP(:, i) = LLR(Dp_z{aP_D, i}) - LLR(D0_z{bP_D, i});
end

X_all = [all_dA; all_dP];
Y_all = [dMADRS_A; dMADRS_P];

% Elastic Net cannot handle NaNs. Isolate complete cases.
valid_rows = ~any(isnan(X_all), 2) & ~isnan(Y_all);
X_clean = X_all(valid_rows, :);
Y_clean = Y_all(valid_rows);
N_clean = length(Y_clean);

fprintf('Running Nested Elastic Net CV on %d complete cases...\n', N_clean);

y_pred_enet = nan(N_clean, 1);
coeff_matrix = zeros(n_markers, N_clean);

% Outer LOOCV loop; the inner 5-fold CV (below) tunes lambda
for fold = 1:N_clean
    train_idx = setdiff(1:N_clean, fold);
    test_idx = fold;

    X_train = X_clean(train_idx, :);
    Y_train = Y_clean(train_idx);
    X_test  = X_clean(test_idx, :);

    % Freeze seed ONLY for the inner CV folds so the lambda tuning is deterministic
    rng(fold);

    % Inner CV: Tune lambda using 5-fold CV. Alpha = 0.5 specifies Elastic Net.
    [B, FitInfo] = lasso(X_train, Y_train, 'CV', 5, 'Alpha', 0.5, 'Standardize', true);

    % Select the lambda that minimizes Mean Squared Error on the inner folds
    idx_opt = FitInfo.IndexMinMSE;
    opt_B = B(:, idx_opt);
    opt_Intercept = FitInfo.Intercept(idx_opt);

    % Record the chosen coefficients for this specific LOOCV fold
    coeff_matrix(:, fold) = opt_B;

    % Predict the hidden test patient
    y_pred_enet(test_idx) = X_test * opt_B + opt_Intercept;
end

% Final multivariate test R^2
SS_res_enet = sum((Y_clean - y_pred_enet).^2);
SS_tot_enet = sum((Y_clean - mean(Y_clean)).^2);
r2_cv_enet = 1 - (SS_res_enet / SS_tot_enet);

fprintf('Nested LOOCV Elastic Net Test R^2: %.3f\n', r2_cv_enet);

% --- Plotting the Selected Features ---

% Calculate how often each marker was selected (coefficient != 0)
selection_freq = sum(coeff_matrix ~= 0, 2) / N_clean * 100;

% Sort markers by selection frequency (most selected first)
[sorted_freq, sort_idx_enet] = sort(selection_freq, 'descend');
sorted_markers = ind_markers(sort_idx_enet);
sorted_coeffs = coeff_matrix(sort_idx_enet, :);

figure('Position', [150 150 1200 800], 'Color', '#DAF2FB');
t_enet = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_enet, 'Elastic-net: \Delta MADRS from all marker \Delta LLR');

% Subplot 1: Selection Frequency
ax1 = nexttile; hold on;
bar(1:n_markers, sorted_freq, 'FaceColor', [0.4 0.6 0.8], 'EdgeColor', 'k');
yline(50, 'k--', 'LineWidth', 1); % 50% reference line
set(gca, 'XTick', 1:n_markers, 'XTickLabel', [], 'Color', 'w', ...
    'TickDir', 'out', 'FontSize', 12, 'xlim', [0 n_markers+1], 'ylim', [0 105]);
ylabel('Selection Frequency (%)');
title(sprintf('Elastic Net Feature Selection (Nested CV R^2 = %.2f)', r2_cv_enet));
box off;

% Subplot 2: Coefficient Weights Distribution
ax2 = nexttile; hold on;
% Draw a zero line
yline(0, 'k-', 'LineWidth', 1);

% Plot boxplots for the coefficients (ignoring exact zeros for clarity on the spread)
for k = 1:n_markers
    % Extract non-zero weights for this marker across all folds
    weights = sorted_coeffs(k, sorted_coeffs(k,:) ~= 0);

    if ~isempty(weights)
        % Add slight jitter for visual clarity if there are few points
        x_jitter = k + (rand(1, length(weights)) - 0.5) * 0.2;
        scatter(x_jitter, weights, 15, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

        % Plot the mean weight as a red dash
        plot([k-0.3, k+0.3], [mean(weights), mean(weights)], 'r-', 'LineWidth', 2);
    end
end

set(gca, 'XTick', 1:n_markers, 'XTickLabel', sorted_markers, 'TickLabelInterpreter', 'none', ...
    'Color', 'w', 'TickDir', 'out', 'FontSize', 11, 'xlim', [0 n_markers+1]);
xtickangle(90);
ylabel('Coefficient Weight');
box off;

%% 16. Pooled distribution per marker (normality / skewness check)
% Wide figure so the one-row subplots aren't squished
figure('Position', [100, 100, 1800, 400]);

% Extract all numeric data directly from the 'all' table
all_numeric = full_sheet(:, vartype('numeric'));

% Loop through each marker to plot in its own subplot
for i = 1:n_markers
    % Create a subplot in a 1-row grid
    subplot(1, n_markers, i);
    hold on;

    % Extract current marker data and remove NaNs
    y_vals = all_numeric{:, i};
    y_vals(isnan(y_vals)) = [];

    % Generate random horizontal jitter around a fixed center (x = 1)
    % Jitter width is set to 0.4 (-0.2 to +0.2)
    x_jitter = 1 + 0.4 * (rand(length(y_vals), 1) - 0.5);

    % Plot as small, semi-transparent dots
    plot(x_jitter, y_vals, 'o', 'Color', [0.3 0.3 0.3 0.4], ...
        'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerSize', 1);

    % Formatting and Aesthetics for this individual subplot
    % Put the marker name as the title instead of an X-tick label to save space
    title(markers{i}, 'Interpreter', 'none', 'FontSize', 10);
    set(gca, 'xlim', [0.5, 1.5], 'XTick', [], 'TickDir', 'out', 'FontSize', 9);
    box off;

    % Add the Y-axis label only to the first plot
    if i == 1
        ylabel('Pooled Values');
    end
end

% Add an overall super-title for the figure
sgtitle('Pooled Distribution per Marker (Independent Y-Axes)');

%% 17. Correlation matrix of all biomarkers
% All subjects, raw (log-transformed) values; Pearson r is scale-invariant
X = full_sheet{:, ind_markers};
[R, P] = corr(X, 'Type','Pearson', 'Rows','pairwise');
R(logical(eye(n_markers)))=nan;     % blank the diagonal

figure
h = heatmap(ind_markers, ind_markers, R, ...
    'Colormap', parula, ...
    'ColorLimits', [-1 1]);                         % fix color scale
title('Pairwise Pearson correlation (z-scored biomarkers)');
xlabel('Biomarker'); ylabel('Biomarker');

%% 18. Add a third marker (creatinine): 3D boundary + rotation GIF
% Fixed corti_sal/crp/creatinine triplet (not the greedy panel); classify_normals
% draws the 3D boundary here since plotmode is left at its default.
results_3=classify_normals(H0{:,{'corti_sal','crp','creatinine'}},D0{:,{'corti_sal','crp','creatinine'}},'input_type','samp','samp_balance',true,'prior_1',0.5,'samp_opt',0);

% 5-fold CV test accuracy for the 3-marker panel
K_folds = 5;
n_reps = 50;
[cv_err_3, cv_err_sd_3] = cv_classify_error(H0{:,{'corti_sal','crp','creatinine'}}, D0{:,{'corti_sal','crp','creatinine'}}, K_folds, n_reps);

axis normal
xlabel('corti_sal','Interpreter','none', 'Rotation', 0, 'HorizontalAlignment', 'right')
ylabel('CRP', 'Rotation', 0, 'HorizontalAlignment', 'right')
zlabel('creatinine','Interpreter','none', 'Rotation', 0, 'HorizontalAlignment', 'right');

set(gca,'fontsize',13,'xticklabel',{},'yticklabel',{},'zticklabel',{},'color','w')
set(gcf, 'Color', [218 242 251]/255);

title(sprintf('accuracy = %.0f\\%%', 100*(1-cv_err_3)))

% Spin the 3D view and append each frame to an animated GIF
axis vis3d;      % freeze aspect ratio for smooth rotation
filename = 'rotating_plot.gif';
for az = 0:1:360
    view(az, 12);
    drawnow;
    frame = getframe(gcf);
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256, 'nodither');
    if az == 0
        imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.03);   % first frame loops forever
    else
        imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.03);
    end
end

%% 19. Local functions
% These are available to every script section regardless of run order.

% Transform a log-likelihood ratio (LLR) into a posterior probability P(H|x).
function p = sigmoid(x)
    p = 1 ./ (1 + exp(-x));
end

% Degree-1 fit constrained to the hypothesised negative slope: if the fitted
% slope comes out positive (wrong direction) the fit is discarded in favour of
% the intercept-only model, i.e. the training mean. Encoding the direction in
% the model this way lets a single cross-validated R^2 measure both direction
% and predictive strength, with no separate slope-sign test afterwards.
function p = fit_neg_slope(x, y)
    p = polyfit(x, y, 1);
    if p(1) > 0
        p = [0, mean(y)];
    end
end

% Forward greedy selection for the joint-P(H) MADRS regression: repeatedly
% append whichever remaining marker most improves the panel's leave-one-out CV
% R^2. `cache` maps a panel (marker-index set, order-independent) to its
% per-fold P(H) transform, so a panel is never fit twice -- panels of size 1
% are always the full candidate list and get hit on every subsequent call, and
% the true-data run and every permutation share the same cache since the
% transform depends only on marker values, not on y. Returns the chosen
% markers in the order they entered, and the panel's CV R^2 after each
% addition.
function [order, r2_curve] = greedy_forward_r2_joint(Z_all, H0Z_all, y_all, max_panel, cache)
    rest = 1:size(Z_all, 2);
    order = nan(1, max_panel);
    r2_curve = nan(1, max_panel);
    chosen = [];
    for step = 1:max_panel
        r2_check = nan(1, numel(rest));
        for c = 1:numel(rest)
            S = sort([chosen, rest(c)]);
            key = mat2str(S);
            if isKey(cache, key)
                xt = cache(key);
            else
                xt = panel_fold_transform(Z_all, H0Z_all, S);
                cache(key) = xt;
            end
            r2_check(c) = score_from_transform(xt, y_all);
        end
        [r2_curve(step), idx_best] = max(r2_check);
        chosen = [chosen, rest(idx_best)];
        order(step) = rest(idx_best);
        rest(idx_best) = [];
    end
end

% Per-fold P(H|panel) for one marker panel (column set S): healthy is fit once
% from all H0 subjects with every marker in S present (they never enter the
% MADRS CV folds, so this is not leakage); depressed is refit inside each fold
% from that fold's training subjects only, so a held-out subject's own marker
% values never contribute to the joint Gaussian used to transform it.
%
% Subjects are kept if this panel's markers are all present for them, not if
% every one of the 25 markers is -- so a singleton panel uses the same
% subjects (and LOOCV folds) as its Section 11 counterpart, and only larger
% panels lose subjects to missingness. xt.v records which rows of the
% (MADRS-complete) input this panel kept, so its own LOOCV partition, built
% from just those rows, lines up with them.
function xt = panel_fold_transform(Z_all, H0Z_all, S)
    Z = Z_all(:, S);
    H0Z = H0Z_all(:, S);
    k = numel(S);
    reg = 1e-5 * eye(k);

    v = ~any(isnan(Z), 2);
    Z = Z(v, :);

    hv = ~any(isnan(H0Z), 2);
    H0Z = H0Z(hv, :);
    mu_H = mean(H0Z, 1);
    Sigma_H = cov(H0Z) + reg;

    cvp = cvpartition(sum(v), 'LeaveOut');
    xt.v   = v;
    xt.cvp = cvp;
    xt.train = cell(1, cvp.NumTestSets);
    xt.test  = cell(1, cvp.NumTestSets);
    for fold = 1:cvp.NumTestSets
        train_idx = training(cvp, fold);
        test_idx  = test(cvp, fold);
        mu_D = mean(Z(train_idx, :), 1);
        Sigma_D = cov(Z(train_idx, :)) + reg;
        LLR_train = log(mvnpdf(Z(train_idx, :), mu_H, Sigma_H)) - log(mvnpdf(Z(train_idx, :), mu_D, Sigma_D));
        LLR_test  = log(mvnpdf(Z(test_idx, :),  mu_H, Sigma_H)) - log(mvnpdf(Z(test_idx, :),  mu_D, Sigma_D));
        xt.train{fold} = sigmoid(LLR_train);
        xt.test{fold}  = sigmoid(LLR_test);
    end
end

% LOOCV R^2 of MADRS ~ P(H|panel), given a panel's precomputed per-fold
% transform. xt.v selects this panel's subjects out of the full (MADRS-
% complete) y vector, and xt.cvp is the LOOCV partition built from exactly
% those subjects. The direction (higher P(H) -> lower MADRS) is enforced
% fold-by-fold via fit_neg_slope, as in Section 11 Figs 3/4.
function r2 = score_from_transform(xt, y_all)
    y = y_all(xt.v);
    cvp = xt.cvp;
    y_pred = nan(size(y));
    for fold = 1:cvp.NumTestSets
        train_idx = training(cvp, fold);
        test_idx  = test(cvp, fold);
        p_cv = fit_neg_slope(xt.train{fold}, y(train_idx));
        y_pred(test_idx) = polyval(p_cv, xt.test{fold});
    end
    SS_res = sum((y - y_pred).^2);
    SS_tot = sum((y - mean(y)).^2);
    r2 = 1 - SS_res / SS_tot;
end

% 1-tailed null-distribution summary: floor (minimum) plus the 99.9/99/95
% percentile upper bounds. (The lower bound of every band is the same floor.)
function [floor_val, u999, u99, u95] = null_ci(d)
    floor_val = min(d);
    u999 = prctile(d, 99.9);
    u99  = prctile(d, 99);
    u95  = prctile(d, 95);
end

% Draw the three nested null-distribution CI rectangles (99.9/99/95%) centred
% at xc with half-width hw, all anchored to a common floor.
function draw_ci_bands(xc, hw, floor_val, u999, u99, u95)
    xr = [xc-hw, xc+hw, xc+hw, xc-hw];
    c  = [1 .85 .6];
    patch(xr, [floor_val, floor_val, u999, u999], c, 'FaceAlpha', .1, 'EdgeColor', c);
    patch(xr, [floor_val, floor_val, u99,  u99 ], c, 'FaceAlpha', .2, 'EdgeColor', c);
    patch(xr, [floor_val, floor_val, u95,  u95 ], c, 'FaceAlpha', .3, 'EdgeColor', c);
end

% Draw each row's baseline->post-treatment segment as a plain solid line (a
% light tint of the group color) with a solid tip marker at the post
% position -- the static end-state of Section 10's animated trajectories.
function plot_static_trajectory(pre, post, light_color, solid_color, marker_size)
    for i = 1:size(pre, 1)
        if any(isnan(pre(i,:))) || any(isnan(post(i,:))), continue; end
        plot([pre(i,1), post(i,1)], [pre(i,2), post(i,2)], '-', 'Color', light_color, 'LineWidth', 1.25);
    end
    plot(post(:,1), post(:,2), 'o', 'MarkerFaceColor', solid_color, 'MarkerEdgeColor', solid_color, 'MarkerSize', marker_size);
end

% Overlay one group's baseline->post P(H) values directly on its colorbar
% (jittered along x, same [0,1] scale as the colorbar itself) -- the static
% end-state of Section 10's animated colorbar traces. Individual subjects are
% small dots; a fatter marker at x=0.5 shows the group average.
function plot_colorbar_trace(cb, pre_prob, post_prob, line_color, point_edge, point_face, point_size, mean_color)
    n = numel(pre_prob);
    jit = linspace(.1, .9, n)';

    drawnow;  % settle the layout so cb.Position is accurate
    cb_ax = axes('Position', cb.Position, 'Color', 'none', 'XLim', [0 1], 'YLim', [0 1], ...
        'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
    uistack(cb_ax, 'top');
    hold(cb_ax, 'on');

    for i = 1:n
        if isnan(pre_prob(i)) || isnan(post_prob(i)), continue; end
        plot([jit(i), jit(i)], [pre_prob(i), post_prob(i)], '-', 'Color', line_color, 'LineWidth', .5);
    end
    plot(jit, post_prob(:), 'o', 'Color', point_edge, 'MarkerFaceColor', point_face, 'MarkerSize', point_size);

    plot([0.5, 0.5], [mean(pre_prob, 'omitnan'), mean(post_prob, 'omitnan')], '-', 'Color', mean_color, 'LineWidth', 1.5);
    plot(0.5, mean(post_prob, 'omitnan'), 'o', 'MarkerFaceColor', mean_color, 'MarkerEdgeColor', mean_color, 'MarkerSize', 5);
end

% One quadrant of the 2x2 treatment-trajectory grid (Section 7): background
% P(H) contour and baseline boundary, the other group shown static at baseline
% for context, and this group's trajectories. Axis/colorbar labels are only
% drawn when show_labels is true (matches Section 10, which labels just the
% first quadrant).
function cb = plot_trajectory_subplot(ax65, results_top2, ctx_x, ctx_y, ctx_color, ...
        pre, post, light_color, solid_color, top2, show_labels)
    hold on
    fcontour(@(x,y) arrayfun(@(x0,y0) results_top2.post_1([x0; y0]), x, y), ax65, ...
        'Fill', 'on', 'MeshDensity', 200, 'LevelList', linspace(0, 1, 100));
    cb = colorbarpzn(0, 1, 'full', 0.5, 'colorP', [0.8 0.8 1], 'colorN', [1 0.8 0.8]);
    plot_boundary(results_top2.norm_bd, 2, 'plot_type', 'line', 'line_color', [0 0 0]);

    plot(ctx_x, ctx_y, 'o', 'MarkerFaceColor', ctx_color, 'MarkerEdgeColor', 'none', 'MarkerSize', 3);
    plot_static_trajectory(pre, post, light_color, solid_color, 3);

    axis(ax65); axis square; box on;
    set(gca, 'XTick', [], 'YTick', [], 'fontsize', 13, 'Color', 'w');
    cb.Ticks = 0:.25:1;
    cb.TickLength = 0.03;
    if show_labels
        xlabel(top2{1}, 'Interpreter', 'none'); ylabel(top2{2}, 'Interpreter', 'none');
        cbTitle = title(cb, '$P(H | \mathbf{m})$');
        cbTitle.Interpreter = 'latex';
        cb.TickLabels = {'0%', '25%', '50%', '75%', '100%'};
    else
        cb.TickLabels = {};
    end
end


% Permutation null distribution of class-balanced CV accuracy: shuffle the
% pooled rows, re-split into the original group sizes, and CV-classify.
function null_accs = perm_null(H, D, K_folds, n_perms)
    null_accs = nan(1, n_perms);
    all_data = [H; D];
    n_H = size(H, 1);
    for p = 1:n_perms
        idx_perm = randperm(size(all_data, 1));
        fake_H = all_data(idx_perm(1:n_H), :);
        fake_D = all_data(idx_perm(n_H+1:end), :);
        null_accs(p) = 1 - cv_classify_error(fake_H, fake_D, K_folds, 1);
    end
end

% Greedy accumulation of restoration markers, run num_runs times to get a
% stable marker ordering by mean greedy rank, then a cumulative-excess figure.
% preA/postA/preP/postP/ref are numeric matrices (rows x n_markers, in
% ind_markers column order). ref is the opposite-class baseline z-data.
% Returns the markers sorted by mean greedy rank (best first).
function stable_markers_rec = greedy_restoration(preA_arr, postA_arr, preP_arr, postP_arr, ref_arr, is_eval_D, ind_markers, num_runs, K_folds, n_perms)
    n_markers = numel(ind_markers);

    % --- The greedy accumulation loop ---
    all_ranks_rec = zeros(num_runs, n_markers);
    for i_run = 1:num_runs
        fprintf('Running greedy restoration accumulation %d / %d...\n', i_run, num_runs);

        rest_markers = ind_markers;
        rest_ref  = ref_arr;
        rest_preA = preA_arr; rest_postA = postA_arr;
        rest_preP = preP_arr; rest_postP = postP_arr;

        greedy_markers_run = cell(1, n_markers);
        greedy_ref = [];
        greedy_preA = []; greedy_postA = [];
        greedy_preP = []; greedy_postP = [];

        for i_greedy = 1:n_markers
            num_markers_rest = size(rest_ref, 2);
            pi_check = nan(1, num_markers_rest);

            for i_check = 1:num_markers_rest
                test_ref   = [greedy_ref, rest_ref(:, i_check)];
                test_preA  = [greedy_preA, rest_preA(:, i_check)];
                test_postA = [greedy_postA, rest_postA(:, i_check)];
                test_preP  = [greedy_preP, rest_preP(:, i_check)];
                test_postP = [greedy_postP, rest_postP(:, i_check)];

                [pi_mu, ~, ~] = excess_restore_cv(test_preA, test_postA, test_preP, test_postP, test_ref, is_eval_D, K_folds, 1, 0);
                pi_check(i_check) = pi_mu;
            end

            % For recovery, we want to MAXIMIZE the excess P(H) shift
            [~, idx_best] = max(pi_check);

            greedy_markers_run{i_greedy} = rest_markers{idx_best};
            greedy_ref   = [greedy_ref, rest_ref(:, idx_best)];
            greedy_preA  = [greedy_preA, rest_preA(:, idx_best)];
            greedy_postA = [greedy_postA, rest_postA(:, idx_best)];
            greedy_preP  = [greedy_preP, rest_preP(:, idx_best)];
            greedy_postP = [greedy_postP, rest_postP(:, idx_best)];

            rest_markers(idx_best) = [];
            rest_ref(:, idx_best)   = [];
            rest_preA(:, idx_best) = []; rest_postA(:, idx_best) = [];
            rest_preP(:, idx_best) = []; rest_postP(:, idx_best) = [];
        end

        [~, ranks_this_run] = ismember(ind_markers, greedy_markers_run);
        all_ranks_rec(i_run, :) = ranks_this_run;
    end

    % --- Analysis & Reordering ---
    mean_ranks_rec = mean(all_ranks_rec, 1);
    rank_ci_rec = prctile(all_ranks_rec, [12.5 87.5], 1);  % 75% interval, empirical (ranks are bounded, not normal)
    [sorted_mean_ranks_rec, sort_idx_rec] = sort(mean_ranks_rec);
    sorted_rank_ci_lo_rec = sorted_mean_ranks_rec - rank_ci_rec(1, sort_idx_rec);
    sorted_rank_ci_hi_rec = rank_ci_rec(2, sort_idx_rec) - sorted_mean_ranks_rec;
    stable_markers_rec = ind_markers(sort_idx_rec);

    % --- Evaluate Cumulative Performance of Stable Sequence ---
    [~, stable_idx_rec] = ismember(stable_markers_rec, ind_markers);
    stable_ref  = ref_arr(:, stable_idx_rec);
    stable_preA = preA_arr(:, stable_idx_rec);
    stable_postA = postA_arr(:, stable_idx_rec);
    stable_preP = preP_arr(:, stable_idx_rec);
    stable_postP = postP_arr(:, stable_idx_rec);

    cum_pi_mean = nan(1, n_markers);
    cum_pi_ci_lo = nan(1, n_markers);
    cum_pi_ci_hi = nan(1, n_markers);
    cum_null_floor = nan(1, n_markers);
    cum_null_upper_999 = nan(1, n_markers);
    cum_null_upper_99 = nan(1, n_markers);
    cum_null_upper_95 = nan(1, n_markers);

    for i = 1:n_markers
        eval_ref   = stable_ref(:, 1:i);
        eval_preA  = stable_preA(:, 1:i);
        eval_postA = stable_postA(:, 1:i);
        eval_preP  = stable_preP(:, 1:i);
        eval_postP = stable_postP(:, 1:i);

        fprintf('Evaluating cumulative step %d/%d with %d permutations...\n', i, n_markers, n_perms);
        [mean_p, ~, null_pi, rep_p] = excess_restore_cv(eval_preA, eval_postA, eval_preP, eval_postP, eval_ref, is_eval_D, K_folds, 50, n_perms);
        cum_pi_mean(i) = mean_p;
        pi_ci = prctile(rep_p, [12.5 87.5]);
        cum_pi_ci_lo(i) = pi_ci(1); cum_pi_ci_hi(i) = pi_ci(2);

        [cum_null_floor(i), cum_null_upper_999(i), cum_null_upper_99(i), cum_null_upper_95(i)] = null_ci(null_pi);
    end

    % Marker labels (significance is shown by the bottom-tile marker colour)
    custom_labels = strrep(stable_markers_rec, '_', '\_');

    % --- Plotting ---
    figure('Color', '#DAF2FB');
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    if is_eval_D, group_name = 'Depressed'; else, group_name = 'Healthy'; end
    sgtitle(sprintf('%s restoration: greedy marker accumulation', group_name));
    x_axis = 1:n_markers;

    % Top tile: mean greedy rank (dots + black error bars)
    ax1 = nexttile; hold on;
    errorbar(x_axis, sorted_mean_ranks_rec, sorted_rank_ci_lo_rec, sorted_rank_ci_hi_rec, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 0.75, 'CapSize', 0);
    plot(x_axis, sorted_mean_ranks_rec, '-ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    ylabel('mean rank');
    set(gca, 'XTick', x_axis, 'XTickLabel', custom_labels, 'TickLabelInterpreter', 'tex', ...
        'TickDir', 'out', 'FontSize', 13, 'xlim', [0 n_markers+1], 'ylim', [0 n_markers+1], 'ytick', [1 n_markers], 'Color', 'w');
    xtickangle(90);
    box off;

    % Bottom tile: cumulative excess vs the permutation null bands
    ax2 = nexttile; hold on;
    plot([0, n_markers+1], [0, 0], '-', 'LineWidth', .5, 'Color', [1 .85 .6]);
    rect_w = 0.8;
    for k = 1:n_markers
        draw_ci_bands(x_axis(k), rect_w/2, cum_null_floor(k), ...
            cum_null_upper_999(k), cum_null_upper_99(k), cum_null_upper_95(k));
    end
    % Black line and error bars throughout; markers filled black where
    % significant (beats the 95% null), hollow (white face) otherwise
    sig = cum_pi_mean > cum_null_upper_95;
    plot(x_axis, cum_pi_mean, '-k', 'LineWidth', 1);
    errorbar(x_axis(sig), cum_pi_mean(sig), cum_pi_mean(sig) - cum_pi_ci_lo(sig), cum_pi_ci_hi(sig) - cum_pi_mean(sig), ...
        'ok', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', 1, 'CapSize', 0);
    errorbar(x_axis(~sig), cum_pi_mean(~sig), cum_pi_mean(~sig) - cum_pi_ci_lo(~sig), cum_pi_ci_hi(~sig) - cum_pi_mean(~sig), ...
        'ok', 'LineStyle', 'none', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 4, 'LineWidth', 1, 'CapSize', 0);
    box off;
    ylabel('combined excess restoration');
    % ylim sized asymmetrically: top must clear the null-band ceiling
    % (u999)/CI/point estimates without clipping; bottom may crop the
    % noisy, single-min null floor rather than waste space on it.
    cum_hi = max([cum_null_upper_999, cum_pi_ci_hi, cum_pi_mean], [], 'omitnan');
    cum_lo = min([cum_pi_ci_lo, cum_pi_mean], [], 'omitnan');
    cum_pad = 0.1 * (cum_hi - cum_lo);
    set(gca, 'XTick', [], 'xlim', [0 n_markers+1], 'ylim', [cum_lo - cum_pad, cum_hi + cum_pad], ...
        'YTick', [0 round(cum_hi, 1)], 'FontSize', 13, 'Color', 'w');
    linkaxes([ax1, ax2], 'x');
end

% Individual-marker restoration figure (Delta P(H) per marker + CV excess).
% evalPre/Post A/P are numeric matrices (rows x n_markers, ind_markers order)
% for the cohort being evaluated. baselineH/baselineD are the full baseline
% z-matrices used for the per-marker LLR and the joint covariance. top2_idx
% are the two marker-column indices (from the greedy ranking) for the joint.
% cA/cP are the Ayahuasca/Placebo colors.
function restoration_figure(evalPreA, evalPostA, evalPreP, evalPostP, baselineH, baselineD, is_eval_D, top2_idx, cA, cP, ind_markers)
    n_markers = numel(ind_markers);
    n_perms = 1000;

    % Reference data for the excess CV is the opposite-class baseline.
    if is_eval_D, ref_all = baselineH; else, ref_all = baselineD; end

    scatter_A = cell(1, n_markers); scatter_P = cell(1, n_markers);
    cross_A = cell(1, n_markers);   cross_P = cell(1, n_markers);
    mu_A = zeros(1, n_markers);     mu_P = zeros(1, n_markers);
    pi_scores = nan(1, n_markers);
    pi_scores_ci_lo = nan(1, n_markers);
    pi_scores_ci_hi = nan(1, n_markers);
    pi_null_floor = nan(1, n_markers);
    pi_null_upper_999 = nan(1, n_markers);
    pi_null_upper_99 = nan(1, n_markers);
    pi_null_upper_95 = nan(1, n_markers);

    for i = 1:n_markers
        md_H = mean(baselineH(:, i), 'omitnan'); sd_H = std(baselineH(:, i), 'omitnan');
        md_D = mean(baselineD(:, i), 'omitnan'); sd_D = std(baselineD(:, i), 'omitnan');
        LLR = @(z) log(normpdf(z, md_H, sd_H)) - log(normpdf(z, md_D, sd_D));

        LLR_bA = LLR(evalPreA(:, i)); LLR_aA = LLR(evalPostA(:, i));
        LLR_bP = LLR(evalPreP(:, i)); LLR_aP = LLR(evalPostP(:, i));

        dA = sigmoid(LLR_aA) - sigmoid(LLR_bA);
        dP = sigmoid(LLR_aP) - sigmoid(LLR_bP);

        cA_cross = (LLR_bA < 0) & (LLR_aA > 0);
        cP_cross = (LLR_bP < 0) & (LLR_aP > 0);

        validA = ~isnan(dA); validP = ~isnan(dP);

        scatter_A{i} = dA(validA); cross_A{i} = cA_cross(validA);
        scatter_P{i} = dP(validP); cross_P{i} = cP_cross(validP);

        mu_A(i) = mean(scatter_A{i});  mu_P(i) = mean(scatter_P{i});

        if length(scatter_A{i}) > 2 && length(scatter_P{i}) > 2
            fprintf('Evaluating CV Excess and permutations for marker %d/%d...\n', i, n_markers);
            [cv_pi, ~, null_pi, rep_pi] = excess_restore_cv(evalPreA(:, i), evalPostA(:, i), evalPreP(:, i), evalPostP(:, i), ref_all(:, i), is_eval_D, 5, 50, n_perms);
            pi_scores(i) = cv_pi;
            pi_ci = prctile(rep_pi, [12.5 87.5]);
            pi_scores_ci_lo(i) = pi_ci(1); pi_scores_ci_hi(i) = pi_ci(2);
            [pi_null_floor(i), pi_null_upper_999(i), pi_null_upper_99(i), pi_null_upper_95(i)] = null_ci(null_pi);
        end
    end

    % --- Joint distribution for the greedy top-2 markers ---
    t_idx = top2_idx;
    base_D = baselineD(:, t_idx);
    clean_D = base_D(~any(isnan(base_D), 2), :);
    C_D_reg = cov(clean_D) + 1e-5 * eye(2);

    base_H = baselineH(:, t_idx);
    clean_H = base_H(~any(isnan(base_H), 2), :);
    C_H_reg = cov(clean_H) + 1e-5 * eye(2);

    LLR_2D = @(X) log(mvnpdf(X, mean(clean_H), C_H_reg)) - log(mvnpdf(X, mean(clean_D), C_D_reg));

    preA = evalPreA(:, t_idx); postA = evalPostA(:, t_idx);
    vA = ~any(isnan(preA), 2) & ~any(isnan(postA), 2);
    preA_LLR = LLR_2D(preA(vA, :)); postA_LLR = LLR_2D(postA(vA, :));
    dA_joint = sigmoid(postA_LLR) - sigmoid(preA_LLR);
    cross_A_joint = (preA_LLR < 0) & (postA_LLR > 0);

    preP = evalPreP(:, t_idx); postP = evalPostP(:, t_idx);
    vP = ~any(isnan(preP), 2) & ~any(isnan(postP), 2);
    preP_LLR = LLR_2D(preP(vP, :)); postP_LLR = LLR_2D(postP(vP, :));
    dP_joint = sigmoid(postP_LLR) - sigmoid(preP_LLR);
    cross_P_joint = (preP_LLR < 0) & (postP_LLR > 0);

    muA_j = mean(dA_joint); muP_j = mean(dP_joint);

    fprintf('Evaluating CV Excess and permutations for joint marker...\n');
    [pi_joint, ~, joint_null_pi, rep_pi_joint] = excess_restore_cv(preA(vA,:), postA(vA,:), preP(vP,:), postP(vP,:), ref_all(:, t_idx), is_eval_D, 5, 50, n_perms);
    joint_pi_ci = prctile(rep_pi_joint, [12.5 87.5]);
    joint_pi_ci_lo = joint_pi_ci(1); joint_pi_ci_hi = joint_pi_ci(2);
    [joint_pi_floor, joint_pi_upper_999, joint_pi_upper_99, joint_pi_upper_95] = null_ci(joint_null_pi);

    % --- Sort individual markers by CV score (descending) & prepend joint ---
    [~, all_idx] = sort(pi_scores(:), 'descend');
    all_idx = all_idx(:)';

    all_scat_A = [{dA_joint}, scatter_A(all_idx)];
    all_scat_P = [{dP_joint}, scatter_P(all_idx)];
    all_cross_A = [{cross_A_joint}, cross_A(all_idx)];
    all_cross_P = [{cross_P_joint}, cross_P(all_idx)];

    all_mu_A = [muA_j, mu_A(all_idx)];
    all_mu_P = [muP_j, mu_P(all_idx)];
    all_pi = [pi_joint, pi_scores(all_idx)];
    all_pi_ci_lo = [joint_pi_ci_lo, pi_scores_ci_lo(all_idx)];
    all_pi_ci_hi = [joint_pi_ci_hi, pi_scores_ci_hi(all_idx)];

    all_pi_floor     = [joint_pi_floor, pi_null_floor(all_idx)];
    all_pi_upper_999 = [joint_pi_upper_999, pi_null_upper_999(all_idx)];
    all_pi_upper_99  = [joint_pi_upper_99, pi_null_upper_99(all_idx)];
    all_pi_upper_95  = [joint_pi_upper_95, pi_null_upper_95(all_idx)];

    % Tick labels: black if significant (beats 95% null), grey otherwise;
    % joint and top-2 individual columns are bold only when significant
    N = length(all_pi);
    labels = cell(1, N);
    for k = 1:N
        if k == 1
            safe_name = 'top 2';
        else
            orig_idx = all_idx(k-1);
            safe_name = strrep(ind_markers{orig_idx}, '_', '\_');
        end
        if all_pi(k) > all_pi_upper_95(k)
            if k <= 3
                labels{k} = sprintf('\\bf{\\color{black}%s}', safe_name);
            else
                labels{k} = sprintf('\\color{black}%s', safe_name);
            end
        else
            labels{k} = sprintf('\\color[rgb]{0.7,0.7,0.7}%s', safe_name);
        end
    end

    % --- Unified plotting ---
    x_gap = 2.5;
    x = (1:N) * x_gap;
    figure('Color', '#DAF2FB');
    if is_eval_D, group_name = 'Depressed'; else, group_name = 'Healthy'; end
    sgtitle(sprintf('%s restoration: per-marker \\Delta P(H)', group_name));

    % TOP PLOT (raw full data)
    ax1 = axes('Position', [0.08, 0.40, 0.90, 0.55]); hold on;
    plot([0, x(end) + x_gap], [0, 0], 'k-', 'LineWidth', 0.1);
    for i = 1:N
        xA_cen = x(i) - 0.4;
        xP_cen = x(i) + 0.4;

        std_A = std(all_scat_A{i});
        std_P = std(all_scat_P{i});
        box_w = 0.35;

        patch([xA_cen-box_w, xA_cen+box_w, xA_cen+box_w, xA_cen-box_w], ...
            [all_mu_A(i)-std_A/2, all_mu_A(i)-std_A/2, all_mu_A(i)+std_A/2, all_mu_A(i)+std_A/2], ...
            cA, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

        patch([xP_cen-box_w, xP_cen+box_w, xP_cen+box_w, xP_cen-box_w], ...
            [all_mu_P(i)-std_P/2, all_mu_P(i)-std_P/2, all_mu_P(i)+std_P/2, all_mu_P(i)+std_P/2], ...
            cP, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

        idxA_cross = all_cross_A{i};
        plot(xA_cen * ones(1, sum(idxA_cross)), all_scat_A{i}(idxA_cross), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cA, 'MarkerSize', 4);
        plot(xA_cen * ones(1, sum(~idxA_cross)), all_scat_A{i}(~idxA_cross), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cA, 'MarkerSize', 2);

        idxP_cross = all_cross_P{i};
        plot(xP_cen * ones(1, sum(idxP_cross)), all_scat_P{i}(idxP_cross), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cP, 'MarkerSize', 4);
        plot(xP_cen * ones(1, sum(~idxP_cross)), all_scat_P{i}(~idxP_cross), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cP, 'MarkerSize', 2);
    end

    set(ax1, 'XTick', x, 'XTickLabel', labels, 'TickLabelInterpreter', 'tex', ...
        'xlim', [0, x(end) + x_gap], 'ylim', [-.5 1], 'ytick', [0 1], 'TickDir', 'out', 'fontsize', 13);
    xtickangle(ax1, 90);
    ylabel('restoration');

    % BOTTOM PLOT (CV excess & 1-tailed CIs)
    ax2 = axes('Position', [0.08, 0.08, 0.90, 0.15]); hold on;
    plot([0, x(end) + x_gap], [0, 0], '-', 'color', [1 .85 .6], 'LineWidth', 0.1);
    rect_w_bot = 0.8;
    for k = 1:N
        draw_ci_bands(x(k), rect_w_bot, all_pi_floor(k), ...
            all_pi_upper_999(k), all_pi_upper_99(k), all_pi_upper_95(k));
    end
    % Black line and error bars; markers filled black where significant (beats
    % the 95% null), hollow (white face) otherwise
    sig = all_pi > all_pi_upper_95;
    plot(x, all_pi, '-k', 'LineWidth', 1);
    errorbar(x(sig), all_pi(sig), all_pi(sig) - all_pi_ci_lo(sig), all_pi_ci_hi(sig) - all_pi(sig), ...
        'ok', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'CapSize', 0, 'Color', 'k');
    errorbar(x(~sig), all_pi(~sig), all_pi(~sig) - all_pi_ci_lo(~sig), all_pi_ci_hi(~sig) - all_pi(~sig), ...
        'ok', 'LineStyle', 'none', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 4, 'CapSize', 0, 'Color', 'k');

    % ylim sized to the data actually plotted, asymmetrically: the top must
    % clear the null-band ceiling (u999) and the CI/point estimates, but the
    % bottom is allowed to crop the (noisy, single-min) null floor rather
    % than waste space accommodating it.
    bot_hi = max([all_pi_upper_999, all_pi_ci_hi, all_pi], [], 'omitnan');
    bot_lo = min([all_pi_ci_lo, all_pi], [], 'omitnan');
    bot_pad = 0.1 * (bot_hi - bot_lo);
    set(ax2, 'XTick', [], 'ylim', [bot_lo - bot_pad, bot_hi + bot_pad], 'YTick', [0 round(bot_hi, 2)], ...
        'xlim', [0, x(end) + x_gap], 'fontsize', 13);
    ylabel('excess restoration');
    linkaxes([ax1, ax2], 'x');
end

% Class-balanced K-fold CV classification error (QDA, or linear if requested),
% averaged over n_reps repetitions. Rows with any NaN are dropped first.
function [mean_err, sd_err, rep_errs] = cv_classify_error(dist_1, dist_2, K_folds, n_reps, varargin)
% Parse optional 'linear' flag compactly
idx = find(strcmpi(varargin, 'linear'));
use_linear = ~isempty(idx) && varargin{idx+1};

if nargin < 4 || isempty(n_reps), n_reps = 10; end
if nargin < 3 || isempty(K_folds), K_folds = 5; end

% Drop rows with any NaNs
dist_1 = dist_1(~any(isnan(dist_1), 2), :);
dist_2 = dist_2(~any(isnan(dist_2), 2), :);

all_fold_errs = nan(n_reps, K_folds);

for r = 1:n_reps
    cv_1 = cvpartition(size(dist_1, 1), 'KFold', K_folds);
    cv_2 = cvpartition(size(dist_2, 1), 'KFold', K_folds);

    for k = 1:K_folds
        train_1 = dist_1(cv_1.training(k), :); test_1 = dist_1(cv_1.test(k), :);
        train_2 = dist_2(cv_2.training(k), :); test_2 = dist_2(cv_2.test(k), :);

        try
            reg = 1e-6 * eye(size(train_1, 2));
            m1 = mean(train_1); cov1 = cov(train_1) + reg;
            m2 = mean(train_2); cov2 = cov(train_2) + reg;
            if use_linear
                norm_linear_bd=best_linear_classifier(m1, cov1, m2, cov2);
                q1 = norm_linear_bd.q1(:); q0 = norm_linear_bd.q0;

                err_1 = sum((test_1 * q1 + q0) <= 0) / size(test_1, 1);
                err_2 = sum((test_2 * q1 + q0) >= 0) / size(test_2, 1);
            else
                err_1 = sum((log(mvnpdf(test_1, m1, cov1)) - log(mvnpdf(test_1, m2, cov2))) <= 0) / size(test_1, 1);
                err_2 = sum((log(mvnpdf(test_2, m1, cov1)) - log(mvnpdf(test_2, m2, cov2))) >= 0) / size(test_2, 1);
            end

            all_fold_errs(r, k) = (err_1 + err_2) / 2;
        catch
            all_fold_errs(r, k) = 1.0;
        end
    end
end

rep_errs = mean(all_fold_errs, 2);
mean_err = mean(rep_errs);
sd_err = std(rep_errs);
end

% Cross-validated Delta P(H) "excess" of Ayahuasca over Placebo, with an
% optional permutation null. is_eval_D = true evaluates Depressed recovery
% (ref_data = full Healthy baseline); false evaluates Healthy (ref = Depressed).
function [mean_score, sd_score, null_scores, rep_scores] = excess_restore_cv(preA, postA, preP, postP, ref_data, is_eval_D, K, n_reps, n_perms)
    function res = get_LLR(X, md_H, C_H, md_D, C_D)
        % Safely evaluate LLR comparing Healthy vs Depressed
        res = nan(size(X, 1), 1);
        valid = ~any(isnan(X), 2);
        if size(X, 2) == 1
            res(valid) = log(normpdf(X(valid,:), md_H, C_H)) - log(normpdf(X(valid,:), md_D, C_D));
        else
            res(valid) = log(mvnpdf(X(valid,:), md_H, C_H)) - log(mvnpdf(X(valid,:), md_D, C_D));
        end
    end

    % 1. Learn the STATIC parameters from the Reference Data once
    clean_ref = ref_data(~any(isnan(ref_data), 2), :);
    if size(clean_ref, 2) == 1
        md_ref = mean(clean_ref); C_ref = std(clean_ref);
    else
        md_ref = mean(clean_ref);
        C_ref = cov(clean_ref);
        C_ref = C_ref + 1e-5 * eye(size(C_ref, 1)); % Ridge Regularization
    end

    function score = eval_folds(fk_preA, fk_postA, fk_preP, fk_postP, idxA, idxP)
        fold_scores = zeros(1, K);
        for k_fold = 1:K
            trA = training(idxA, k_fold); teA = test(idxA, k_fold);
            trP = training(idxP, k_fold); teP = test(idxP, k_fold);

            % 2. Learn the DYNAMIC parameters from the CV Training Fold
            train_base = [fk_preA(trA, :); fk_preP(trP, :)];
            train_base = train_base(~any(isnan(train_base), 2), :);

            if size(train_base, 2) == 1
                md_fold = mean(train_base); C_fold = std(train_base);
            else
                md_fold = mean(train_base);
                C_fold = cov(train_base);
                C_fold = C_fold + 1e-5 * eye(size(C_fold, 1)); % Ridge Regularization
            end

            % 3. Assign H and D parameters correctly based on who we are evaluating
            if is_eval_D
                md_H = md_ref;  C_H = C_ref;
                md_D = md_fold; C_D = C_fold;
            else
                md_H = md_fold; C_H = C_fold;
                md_D = md_ref;  C_D = C_ref;
            end

            % Evaluate LLR and Sigmoid on TEST data only
            dA_test = sigmoid(get_LLR(fk_postA(teA, :), md_H, C_H, md_D, C_D)) - sigmoid(get_LLR(fk_preA(teA, :), md_H, C_H, md_D, C_D));
            dP_test = sigmoid(get_LLR(fk_postP(teP, :), md_H, C_H, md_D, C_D)) - sigmoid(get_LLR(fk_preP(teP, :), md_H, C_H, md_D, C_D));

            % Metric: Excess of Aya shift over max(0, Placebo shift)
            fold_scores(k_fold) = mean(dA_test, 'omitnan') - max(0, mean(dP_test, 'omitnan'));
        end
        score = mean(fold_scores, 'omitnan');
    end

    N_A = size(preA, 1); N_P = size(preP, 1);

    % Actual Cross-Validation
    scores = zeros(1, n_reps);
    for r = 1:n_reps
        idxA = cvpartition(N_A, 'KFold', K);
        idxP = cvpartition(N_P, 'KFold', K);
        scores(r) = eval_folds(preA, postA, preP, postP, idxA, idxP);
    end
    mean_score = mean(scores);
    sd_score = std(scores);
    rep_scores = scores;

    % Permutation Test
    null_scores = [];
    if n_perms > 0
        null_scores = zeros(1, n_perms);
        all_pre = [preA; preP];
        all_post = [postA; postP];
        for p = 1:n_perms
            perm_idx = randperm(size(all_pre, 1));
            fk_preA = all_pre(perm_idx(1:N_A), :); fk_postA = all_post(perm_idx(1:N_A), :);
            fk_preP = all_pre(perm_idx(N_A+1:end), :); fk_postP = all_post(perm_idx(N_A+1:end), :);

            idxA = cvpartition(N_A, 'KFold', K);
            idxP = cvpartition(N_P, 'KFold', K);
            null_scores(p) = eval_folds(fk_preA, fk_postA, fk_preP, fk_postP, idxA, idxP);
        end
    end
end