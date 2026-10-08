clear;
clc;
close all;

params = global_parameters_balanced();

T = readtable('posterior0915.csv');
theta = table2array(T);

rng(1, 'twister');

n_keep = 10;
idx_keep = randperm(size(theta,1), n_keep);
theta = theta(idx_keep,:);

tspan = [0, 60*60*24];

options = odeset( ...
    'RelTol', 1e-4, ...
    'AbsTol', 1e-6);


%% Generic initial condition for baseline calculation

y0 = [
    params.pyruvate_initial;
    params.oxaloacetate_initial;
    params.acetyl_coa_initial;
    params.citrate_initial;
    params.cis_aconitate_initial;
    params.isocitrate_initial;
    params.aketo_initial;
    params.succinyl_coa_initial;
    params.succinate_initial;
    params.fumarate_initial;
    params.acetyl_p_initial;
    params.nad_pool / (1 + params.nadh_ratio);
    params.nad_pool * params.nadh_ratio / (1 + params.nadh_ratio)
];


%% Start parallel pool

if isempty(gcp('nocreate'))
    parpool;
end


%% NADH-response settings

% Direct fold changes in NADH concentration relative to baseline
% An odd number of points ensures fold = 1 is included exactly
nadh_fold_list = logspace(-1, 1, 11);

beta_list = [0, 0.25, 0.50, 0.75, 1.00];

n_grid = numel(nadh_fold_list);
n_beta = numel(beta_list);

% I_beta / I_base
inhibition_relative = nan(n_keep, n_beta, n_grid);

baseline_nadh       = nan(n_keep,1);
baseline_inhibition = nan(n_keep,1);
baseline_valid      = false(n_keep,1);


%% Calculate baseline steady state and inhibition curves

parfor i = 1:n_keep

    p = params;

    p.enzyme_pdh       = theta(i,1);
    p.enzyme_aco       = theta(i,2);
    p.enzyme_akd       = theta(i,3);
    p.enzyme_scs       = theta(i,4);
    p.enzyme_sd        = theta(i,5);
    p.enzyme_pk        = theta(i,6);
    p.enzyme_cs        = theta(i,7);
    p.enzyme_id        = theta(i,8);
    p.enzyme_pta       = theta(i,9);
    p.acetate_constant = theta(i,10);
    p.enzyme_md        = theta(i,11);
    p.enzyme_f         = theta(i,12);
    p.enzyme_ack       = theta(i,13);

    p.k = theta(i,14);
    alpha0 = theta(i,15);

    inhibition_i = nan(n_beta, n_grid);

    I_base = NaN;
    nadh_base = NaN;
    baseline_ok = false;


    %% 1. Calculate this posterior draw's baseline steady state

    p_base = p;
    p_base.alpha = alpha0;

    % beta = 1 reproduces the original inhibition function
    p_base.beta = 1;

    % At beta = 1, the reference value cancels
    p_base.baseline_pdh = 1;

    % Initialise NADH/NAD+ consistently with alpha0
    y0_base = y0;

    y0_base(12) = ...
        p_base.nad_pool / (1 + alpha0);

    y0_base(13) = ...
        p_base.nad_pool * alpha0 / (1 + alpha0);

    try
        [~, baseline_solution] = ode15s( ...
            @(t,y) function_fsolve(y, p_base), ...
            tspan, ...
            y0_base, ...
            options);

        y_base = baseline_solution(end,:)';
        dydt_base = function_fsolve(y_base, p_base);

        if norm(dydt_base, inf) <= 1e-6

            nadh_base = y_base(13);

            I_base = nadh_inhibition( ...
                nadh_base, ...
                p_base.ki_nadh_pdh);

            baseline_ok = ...
                isfinite(nadh_base) && ...
                nadh_base > 0 && ...
                isfinite(I_base) && ...
                I_base > 0;
        end

    catch
        baseline_ok = false;
    end


    %% 2. Evaluate the inhibition function directly over NADH

    if baseline_ok

        for b = 1:n_beta

            beta = beta_list(b);

            for j = 1:n_grid

                % Directly vary NADH concentration relative to baseline
                nadh_test = ...
                    nadh_fold_list(j) * nadh_base;

                I_raw = nadh_inhibition( ...
                    nadh_test, ...
                    p_base.ki_nadh_pdh);

                % I_beta / I_base
                %
                % I_beta = I_base * (I_raw/I_base)^beta
                %
                % Therefore:
                %
                % I_beta/I_base = (I_raw/I_base)^beta
                inhibition_i(b,j) = ...
                    (I_raw / I_base).^beta;
            end
        end
    end

    inhibition_relative(i,:,:) = reshape( ...
        inhibition_i, ...
        1, n_beta, n_grid);

    baseline_nadh(i) = nadh_base;
    baseline_inhibition(i) = I_base;
    baseline_valid(i) = baseline_ok;
end


%% Report baseline convergence

fprintf( ...
    'Valid baseline steady states: %d/%d (%.1f%%)\n', ...
    sum(baseline_valid), ...
    n_keep, ...
    100 * mean(baseline_valid));

valid_nadh = baseline_nadh( ...
    baseline_valid & ...
    isfinite(baseline_nadh) & ...
    baseline_nadh > 0);

if ~isempty(valid_nadh)

    fprintf( ...
        'Median baseline NADH: %.5g mM\n', ...
        median(valid_nadh));

    q_nadh = prctile(valid_nadh, [25, 75]);

    fprintf( ...
        'Baseline NADH IQR: [%.5g, %.5g] mM\n', ...
        q_nadh(1), ...
        q_nadh(2));
end


%% Calculate posterior summaries

inhibition_median = nan(n_beta, n_grid);
inhibition_q25    = nan(n_beta, n_grid);
inhibition_q75    = nan(n_beta, n_grid);

for b = 1:n_beta

    Yb = squeeze(inhibition_relative(:,b,:));

    for j = 1:n_grid

        values = Yb(:,j);

        values = values( ...
            isfinite(values) & ...
            values > 0);

        if numel(values) >= 8

            inhibition_median(b,j) = median(values);

            quartiles = prctile(values, [25, 75]);

            inhibition_q25(b,j) = quartiles(1);
            inhibition_q75(b,j) = quartiles(2);
        end
    end
end


%% Create figure

[fig1, ax1] = make_journal_figure(15, 14);

hold(ax1, 'on');

% Grey: no NADH-dependent response
% Dark blue: original inhibition response
colors = [
    0.65 0.65 0.65;
    0.58 0.72 0.82;
    0.32 0.56 0.73;
    0.13 0.38 0.61;
    0.02 0.18 0.34
];

beta_labels = {
    '\beta = 0'
    '\beta = 0.25'
    '\beta = 0.50'
    '\beta = 0.75'
    '\beta = 1'
};


%% Draw IQR bands first

for b = n_beta:-1:1

    keep = ...
        isfinite(inhibition_median(b,:)) & ...
        isfinite(inhibition_q25(b,:)) & ...
        isfinite(inhibition_q75(b,:)) & ...
        inhibition_q25(b,:) > 0;

    x = nadh_fold_list(keep);
    y25 = inhibition_q25(b,keep);
    y75 = inhibition_q75(b,keep);

    if ~isempty(x)

        fill( ...
            ax1, ...
            [x, fliplr(x)], ...
            [y25, fliplr(y75)], ...
            colors(b,:), ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.05, ...
            'HandleVisibility', 'off');
    end
end


%% Baseline reference lines

xline( ...
    ax1, ...
    1, ...
    '--', ...
    'Color', [0.72 0.72 0.72], ...
    'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

yline( ...
    ax1, ...
    1, ...
    '--', ...
    'Color', [0.72 0.72 0.72], ...
    'LineWidth', 1.2, ...
    'HandleVisibility', 'off');


%% Draw posterior median curves

line_handles = gobjects(n_beta,1);

for b = 1:n_beta

    keep = isfinite(inhibition_median(b,:));

    line_handles(b) = plot( ...
        ax1, ...
        nadh_fold_list(keep), ...
        inhibition_median(b,keep), ...
        '-', ...
        'Color', colors(b,:), ...
        'LineWidth', 2, ...
        'DisplayName', beta_labels{b});
end


%% Axis formatting

set( ...
    ax1, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'Layer', 'top');

xlim(ax1, [0.1, 10]);

xticks(ax1, [0.1, 1, 10]);
xticklabels(ax1, {'0.1', '1', '10'});


ylim(ax1, [0.1, 10]);
yticks(ax1, [0.1, 1, 10]);
yticklabels(ax1, {'0.1', '1', '10'});


xlabel( ...
    ax1, ...
    ['Normalised NADH concentration, ' ...
     'NADH / NADH_0'], ...
    'Interpreter', 'tex');

ylabel( ...
    ax1, ...
    'Relative PDH activity, I_{\beta} / I_0', ...
    'Interpreter', 'tex');

legend( ...
    ax1, ...
    line_handles, ...
    beta_labels, ...
    'Location', 'northeast', ...
    'Box', 'off', ...
    'Interpreter', 'tex');

apply_journal_style( ...
    ax1, ...
    [], ...
    'AxisFontSize', 13, ...
    'LabelFontSize', 13);

ax1.TickLabelInterpreter = 'tex';
ax1.XAxis.MinorTick = 'off';
ax1.YAxis.MinorTick = 'off';


%% Export

% fig1.Color = 'white';
% 
% exportgraphics( ...
%     fig1, ...
%     'beta_inhibition_nadh_response0914.pdf', ...
%     'ContentType', 'vector', ...
%     'BackgroundColor', 'white');