%% ============================================================
%  Robust Profile Likelihood Analysis
% =============================================================

clear;
clc;
close all;


%% ============================================================
%  1. Parameters and initial conditions
% =============================================================

params = global_parameters_balanced();

tspan = [0 60*60*240];

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
    params.nad_pool * (1 - params.nadh_ratio);
    params.nad_pool * params.nadh_ratio
];


%% ============================================================
%  2. Initial guesses and parameter bounds
%
%  All parameters are optimised in log10 space.
% =============================================================

param0 = [
    -2, -2, -2, -2, -2, -2, -2, ...
    -2, -2, -2, -2, -2, 3, -2
];

lb = [
    -3, -3, -3, -3, -3, -3, -3, ...
    -3, -3, -3, -3, -4, 1, -3
];

ub = [
     2,  2,  3,  2,  2,  2,  2, ...
     2,  2,  2,  2,  1,  4, -1
];


%% ============================================================
%  3. Synthetic data
% =============================================================

data.values = [
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
    params.nad_pool * (1 - params.nadh_ratio);
    params.nad_pool * params.nadh_ratio;
    1.0833;
    0.3
];

data.obs_idx = 1:13;

sigma_s = 0.2;
sigma_f = 0.2;


%% ============================================================
%  4. GLOBAL MULTI-START FIT
% =============================================================

opts_fit = optimoptions('fmincon', ...
    'Display', 'off', ...
    'Algorithm', 'interior-point', ...
    'FiniteDifferenceType', 'forward', ...
    'FiniteDifferenceStepSize', 1e-5, ...
    'FunctionTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6, ...
    'MaxIterations', 2000);


n_starts = 10;
n_params = numel(param0);

all_param = NaN(n_starts, n_params);
all_fval = NaN(n_starts, 1);
all_exitflag = NaN(n_starts, 1);

start_points = NaN(n_starts, n_params);

% First start = original initial guess
start_points(1,:) = param0;

% Reproducible random starts
rng(1);

for s = 2:n_starts

    start_points(s,:) = ...
        lb + rand(size(lb)) .* (ub - lb);

end


fprintf('\n');
fprintf('============================================\n');
fprintf('STARTING MULTI-START GLOBAL FIT\n');
fprintf('Number of starts = %d\n', n_starts);
fprintf('============================================\n');


for s = 1:n_starts

    fprintf('Start %d/%d ...\n', s, n_starts);

    p0 = start_points(s,:);

    try

        [p_hat_s, fval_s, exitflag_s] = ...
            fmincon( ...
                @(p) objective_function( ...
                    p, params, data, tspan, y0, sigma_s, sigma_f), ...
                p0, ...
                [], [], [], [], ...
                lb, ub, ...
                [], ...
                opts_fit);

        all_param(s,:) = p_hat_s;
        all_fval(s) = fval_s;
        all_exitflag(s) = exitflag_s;

        fprintf( ...
            '  DONE: fval = %.8g, exitflag = %d\n', ...
            fval_s, exitflag_s);

    catch ME

        fprintf('  FAILED: %s\n', ME.message);

    end

end


%% ------------------------------------------------------------
% Select best valid global solution
% ------------------------------------------------------------

valid_fit = ...
    isfinite(all_fval) & ...
    all_exitflag > 0 & ...
    all_fval < 1e11;

if ~any(valid_fit)

    error('No valid multi-start fit converged.');

end

valid_idx = find(valid_fit);

[best_fval, local_best_idx] = ...
    min(all_fval(valid_fit));

best_idx = valid_idx(local_best_idx);

param_hat = all_param(best_idx,:);
fval_min = all_fval(best_idx);
exitflag_fit = all_exitflag(best_idx);


fprintf('\n');
fprintf('============================================\n');
fprintf('MULTI-START GLOBAL FIT FINISHED\n');
fprintf('Best start = %d/%d\n', best_idx, n_starts);
fprintf('Best objective value = %.8g\n', fval_min);
fprintf('Best exitflag = %d\n', exitflag_fit);
fprintf('============================================\n');

disp('Best estimated log10 parameters:');
disp(param_hat);


result_table = table( ...
    (1:n_starts)', ...
    all_fval, ...
    all_exitflag, ...
    'VariableNames', ...
    {'Start', 'Objective', 'Exitflag'});

result_table = sortrows(result_table, 'Objective');

disp(' ');
disp('Multi-start results sorted by objective:');
disp(result_table);


%% ============================================================
%  5. Update parameter structure using global MLE
% =============================================================

params.enzyme_pk       = 10^param_hat(1);
params.enzyme_pdh      = 10^param_hat(2);
params.enzyme_cs       = 10^param_hat(3);
params.enzyme_aco      = 10^param_hat(4);
params.enzyme_id       = 10^param_hat(5);
params.enzyme_akd      = 10^param_hat(6);
params.enzyme_scs      = 10^param_hat(7);
params.enzyme_sd       = 10^param_hat(8);
params.enzyme_fumarase = 10^param_hat(9);
params.enzyme_md       = 10^param_hat(10);
params.enzyme_pta      = 10^param_hat(11);
params.enzyme_ack      = 10^param_hat(12);

params.k               = 10^param_hat(13);
params.alpha           = 10^param_hat(14);


%% ============================================================
%  6. Check simulation at global MLE
% =============================================================

options_ode = odeset( ...
    'RelTol', 1e-4, ...
    'AbsTol', 1e-6);

[t, metabolite_concentration] = ...
    ode15s( ...
        @(t,y) function_fsolve(y, params), ...
        tspan, ...
        y0, ...
        options_ode);


%% ============================================================
%  7. PROFILE LIKELIHOOD RANGE
%
%  Default:
%
%       MLE - 2 log10 units
%             to
%       MLE + 2 log10 units
%
%  Equivalent to approximately:
%
%       MLE / 100   to   MLE * 100
%
%  but never outside lb / ub.
% =============================================================

n_profile_points = 80;
profile_half_width = 2;

theta_range_all = cell(n_params,1);

for i = 1:n_params

    lower_profile = ...
        max(lb(i), param_hat(i) - profile_half_width);

    upper_profile = ...
        min(ub(i), param_hat(i) + profile_half_width);

    theta_range_all{i} = ...
        linspace( ...
            lower_profile, ...
            upper_profile, ...
            n_profile_points);

end


%% ============================================================
%  8. Allocate profile results
% =============================================================

PL = cell(n_params,1);
PL_exitflag = cell(n_params,1);
PL_opt_params = cell(n_params,1);


%% ============================================================
%  9. Parallel pool
% =============================================================

% delete(gcp('nocreate'));
% 
% parpool('Processes', 4);


%% ============================================================
%  10. ROBUST PROFILE LIKELIHOOD
%
%  Strategy:
%
%  First sweep:
%      centre -> right
%      centre -> left
%
%  Each point compares:
%      1. neighbouring profile optimum
%      2. global MLE
%
%  Then:
%      reverse refinement
%
%  This allows a better likelihood branch found later
%  to propagate backwards.
% =============================================================

for i = 2 %1:n_params

    fprintf('\n');
    fprintf('##################################################\n');
    fprintf('START PARAMETER i = %d\n', i);
    fprintf('MLE log10 value = %.8g\n', param_hat(i));
    fprintf('##################################################\n');

    theta_range = theta_range_all{i};
    n_theta = numel(theta_range);

    idx_free = setdiff(1:n_params, i);

    PL_i = NaN(1, n_theta);
    exitflag_i = NaN(1, n_theta);
    opt_params_i = NaN(n_theta, n_params);


    %% --------------------------------------------------------
    % Profile optimisation options
    % ------------------------------------------------------------

    opts_profile = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'interior-point', ...
        'FiniteDifferenceType', 'forward', ...
        'FiniteDifferenceStepSize', 1e-4, ...
        'FunctionTolerance', 1e-6, ...
        'StepTolerance', 1e-6, ...
        'OptimalityTolerance', 1e-6, ...
        'MaxIterations', 2000);


    %% --------------------------------------------------------
    % Find point closest to MLE
    % ------------------------------------------------------------

    [~, j_center] = ...
        min(abs(theta_range - param_hat(i)));

    fprintf( ...
        'Parameter %d centre: j=%d/%d, theta=%.8g\n', ...
        i, j_center, n_theta, theta_range(j_center));


    %% ========================================================
    %  10A. Centre
    % =============================================================

    fixed_param = theta_range(j_center);

    param_template = param_hat;
    param_template(i) = fixed_param;

    fun_profile = @(p_free) ...
        objective_function( ...
            insert_params( ...
                param_template, ...
                p_free, ...
                idx_free), ...
            params, ...
            data, ...
            tspan, ...
            y0, ...
            sigma_s, ...
            sigma_f);


    starts = {
        param_hat(idx_free)
    };


    [best_free, best_fval, best_flag] = ...
        run_profile_multistart( ...
            fun_profile, ...
            starts, ...
            lb(idx_free), ...
            ub(idx_free), ...
            opts_profile);


    if isfinite(best_fval)

        p_full = insert_params( ...
            param_template, ...
            best_free, ...
            idx_free);

        PL_i(j_center) = best_fval;
        exitflag_i(j_center) = best_flag;
        opt_params_i(j_center,:) = p_full;

    end


    %% ========================================================
    %  10B. Sweep RIGHT
    % =============================================================

    if all(isfinite(opt_params_i(j_center,idx_free)))

        previous_free = ...
            opt_params_i(j_center,idx_free);

    else

        previous_free = param_hat(idx_free);

    end


    for j = (j_center + 1):n_theta

        fixed_param = theta_range(j);

        param_template = param_hat;
        param_template(i) = fixed_param;

        fun_profile = @(p_free) ...
            objective_function( ...
                insert_params( ...
                    param_template, ...
                    p_free, ...
                    idx_free), ...
                params, ...
                data, ...
                tspan, ...
                y0, ...
                sigma_s, ...
                sigma_f);


        starts = {
            previous_free
            param_hat(idx_free)
        };


        [best_free, best_fval, best_flag] = ...
            run_profile_multistart( ...
                fun_profile, ...
                starts, ...
                lb(idx_free), ...
                ub(idx_free), ...
                opts_profile);


        if isfinite(best_fval)

            p_full = insert_params( ...
                param_template, ...
                best_free, ...
                idx_free);

            PL_i(j) = best_fval;
            exitflag_i(j) = best_flag;
            opt_params_i(j,:) = p_full;

            previous_free = best_free;

        end


        if mod(j,10) == 0 || j == n_theta

            fprintf( ...
                'i=%d RIGHT progress %d/%d, f=%.6g\n', ...
                i, j, n_theta, best_fval);

        end

    end


    %% ========================================================
    %  10C. Sweep LEFT
    % =============================================================

    if all(isfinite(opt_params_i(j_center,idx_free)))

        previous_free = ...
            opt_params_i(j_center,idx_free);

    else

        previous_free = param_hat(idx_free);

    end


    for j = (j_center - 1):-1:1

        fixed_param = theta_range(j);

        param_template = param_hat;
        param_template(i) = fixed_param;

        fun_profile = @(p_free) ...
            objective_function( ...
                insert_params( ...
                    param_template, ...
                    p_free, ...
                    idx_free), ...
                params, ...
                data, ...
                tspan, ...
                y0, ...
                sigma_s, ...
                sigma_f);


        starts = {
            previous_free
            param_hat(idx_free)
        };


        [best_free, best_fval, best_flag] = ...
            run_profile_multistart( ...
                fun_profile, ...
                starts, ...
                lb(idx_free), ...
                ub(idx_free), ...
                opts_profile);


        if isfinite(best_fval)

            p_full = insert_params( ...
                param_template, ...
                best_free, ...
                idx_free);

            PL_i(j) = best_fval;
            exitflag_i(j) = best_flag;
            opt_params_i(j,:) = p_full;

            previous_free = best_free;

        end


        if mod(j,10) == 0 || j == 1

            fprintf( ...
                'i=%d LEFT progress %d/%d, f=%.6g\n', ...
                i, j, n_theta, best_fval);

        end

    end


    %% ========================================================
    %  10D. Reverse refinement: RIGHT -> centre
    % =============================================================

    fprintf( ...
        'i=%d: starting RIGHT reverse refinement\n', ...
        i);


    for j = (n_theta - 1):-1:(j_center + 1)

        fixed_param = theta_range(j);

        param_template = param_hat;
        param_template(i) = fixed_param;

        fun_profile = @(p_free) ...
            objective_function( ...
                insert_params( ...
                    param_template, ...
                    p_free, ...
                    idx_free), ...
                params, ...
                data, ...
                tspan, ...
                y0, ...
                sigma_s, ...
                sigma_f);


        starts = {};


        % Existing solution
        if all(isfinite(opt_params_i(j,idx_free)))

            starts{end+1,1} = ...
                opt_params_i(j,idx_free);

        end


        % Solution from right neighbour
        if all(isfinite(opt_params_i(j+1,idx_free)))

            starts{end+1,1} = ...
                opt_params_i(j+1,idx_free);

        end


        % Global MLE
        starts{end+1,1} = ...
            param_hat(idx_free);


        [best_free, best_fval, best_flag] = ...
            run_profile_multistart( ...
                fun_profile, ...
                starts, ...
                lb(idx_free), ...
                ub(idx_free), ...
                opts_profile);


        if isfinite(best_fval) && ...
                (~isfinite(PL_i(j)) || ...
                 best_fval < PL_i(j) - 1e-6)


            old_fval = PL_i(j);

            p_full = insert_params( ...
                param_template, ...
                best_free, ...
                idx_free);

            PL_i(j) = best_fval;
            exitflag_i(j) = best_flag;
            opt_params_i(j,:) = p_full;


            fprintf( ...
                ['i=%d j=%d improved: ' ...
                 '%.6g -> %.6g\n'], ...
                i, j, old_fval, best_fval);

        end

    end


    %% ========================================================
    %  10E. Reverse refinement: LEFT -> centre
    % =============================================================

    fprintf( ...
        'i=%d: starting LEFT reverse refinement\n', ...
        i);


    for j = 2:(j_center - 1)

        fixed_param = theta_range(j);

        param_template = param_hat;
        param_template(i) = fixed_param;

        fun_profile = @(p_free) ...
            objective_function( ...
                insert_params( ...
                    param_template, ...
                    p_free, ...
                    idx_free), ...
                params, ...
                data, ...
                tspan, ...
                y0, ...
                sigma_s, ...
                sigma_f);


        starts = {};


        % Existing solution
        if all(isfinite(opt_params_i(j,idx_free)))

            starts{end+1,1} = ...
                opt_params_i(j,idx_free);

        end


        % Solution from left neighbour
        if all(isfinite(opt_params_i(j-1,idx_free)))

            starts{end+1,1} = ...
                opt_params_i(j-1,idx_free);

        end


        % Global MLE
        starts{end+1,1} = ...
            param_hat(idx_free);


        [best_free, best_fval, best_flag] = ...
            run_profile_multistart( ...
                fun_profile, ...
                starts, ...
                lb(idx_free), ...
                ub(idx_free), ...
                opts_profile);


        if isfinite(best_fval) && ...
                (~isfinite(PL_i(j)) || ...
                 best_fval < PL_i(j) - 1e-6)


            old_fval = PL_i(j);

            p_full = insert_params( ...
                param_template, ...
                best_free, ...
                idx_free);

            PL_i(j) = best_fval;
            exitflag_i(j) = best_flag;
            opt_params_i(j,:) = p_full;


            fprintf( ...
                ['i=%d j=%d improved: ' ...
                 '%.6g -> %.6g\n'], ...
                i, j, old_fval, best_fval);

        end

    end


    %% ========================================================
    % Store profile
    % =============================================================

    PL{i} = PL_i;
    PL_exitflag{i} = exitflag_i;
    PL_opt_params{i} = opt_params_i;


    fprintf('\n');
    fprintf('##################################################\n');
    fprintf('DONE PARAMETER i = %d\n', i);
    fprintf('Successful points = %d / %d\n', ...
        sum(isfinite(PL_i)), n_theta);
    fprintf('##################################################\n');

end


fprintf('\n');
fprintf('==================================================\n');
fprintf('ALL PROFILE LIKELIHOODS FINISHED\n');
fprintf('==================================================\n');


%% ============================================================
%  11. Save results
% =============================================================

save( ...
    'profile_likelihood_results.mat', ...
    'PL', ...
    'PL_exitflag', ...
    'PL_opt_params', ...
    'theta_range_all', ...
    'param_hat', ...
    'fval_min', ...
    'all_param', ...
    'all_fval', ...
    'all_exitflag', ...
    'lb', ...
    'ub');


%% ============================================================
%  12. Parameter labels
% =============================================================
% Replace your plotting sections 12-13 with this script.
% Uses your existing make_journal_figure, apply_journal_style and export_figure.
% PL and theta_range_all must come from the same run.
% theta_range_all already contains log10(parameter value).

parameter_list = { ...
    'E_{PK}', 'E_{PDH}', 'E_{CS}', 'E_{ACO}', ...
    'E_{IDH}', 'E_{AKGDH}', 'E_{SCS}', 'E_{SDH}', ...
    'E_{FUM}', 'E_{MDH}', 'E_{PTA}', 'E_{ACKA}', ...
    'k', '\alpha'};

%% ============================================================
%  13. Plot profile likelihood against log10(parameter value)
% =============================================================

plot_indices = 1:n_params;
[fig, ~] = make_journal_figure(18, 12);
clf(fig);

tl = tiledlayout(fig, 3, 5, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for ii = 1:length(plot_indices)
    i = plot_indices(ii);
    ax = nexttile(tl, ii);

    % Existing grid is already in log10 space: no exponentiation or log10.
    x = theta_range_all{i}(:);
    y = fval_min - PL{i}(:);
    assert(numel(x) == numel(y), ...
        'Profile %d and its grid have different lengths.', i);
    valid = isfinite(x) & isfinite(y); % Negative and zero x values are valid.

    % Preserve missing-point gaps rather than joining across them.

    set(ax, 'XScale', 'linear');

    plot(ax, x(valid), y(valid), '-', ...
    'Color', [0.2 0.2 0.2], ...
    'LineWidth', 1.8);

    hold(ax, 'on');

    yline(ax, -1.92, '--', ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    title(ax, parameter_list{i}, ...
        'Interpreter', 'tex', 'FontWeight', 'normal');

    % First -10 crossing on each side of the peak.
    % x is now log10(parameter value), so interpolate directly in x.
    valid_indices = find(valid);
    xv = x(valid);
    yv = y(valid);
    trimmed_limits = [];

    if numel(xv) >= 2
        assert(all(diff(xv) > 0), 'Profile grid %d must be increasing.', i);
        [y_peak, j_peak] = max(yv);
        x_left = xv(1);
        x_right = xv(end);

        if y_peak > -10
            j_left = find(yv(1:j_peak) <= -10, 1, 'last');
            if ~isempty(j_left)
                a = j_left;
                b = j_left + 1;
                if valid_indices(b) == valid_indices(a) + 1
                    x_left = xv(a) + ...
                        (-10-yv(a))/(yv(b)-yv(a))*(xv(b)-xv(a));
                else
                    x_left = xv(a); % No interpolation across missing points.
                end
            end

            offset = find(yv(j_peak:end) <= -10, 1, 'first');
            if ~isempty(offset)
                b = j_peak + offset - 1;
                a = b - 1;
                if valid_indices(b) == valid_indices(a) + 1
                    x_right = xv(a) + ...
                        (-10-yv(a))/(yv(b)-yv(a))*(xv(b)-xv(a));
                else
                    x_right = xv(b);
                end
            end
        end

        if x_left < x_right
            trimmed_limits = [x_left x_right];
        end
    end

    % Preserve your journal styling helpers and their existing settings.
    apply_journal_style(ax, [], ...
        'AxisFontSize', 16, ...
        'LabelFontSize', 16, ...
        'LineWidth', 2, ...
        'TitleFontWeight', 'normal');

    % Apply scale, limits and ticks AFTER the style helper.
    ax.XScale = 'linear';
    ylim(ax, [-10 0.5]);
    if ~isempty(trimmed_limits)
        xlim(ax, trimmed_limits); % No 10.^ here.
    end
    ax.XMinorTick = 'off';
    ax.XMinorGrid = 'off';
    ax.XTickLabelRotation = 0;

    % Exactly two major x ticks in log10 units.
    xl = xlim(ax);
    width = diff(xl);
    targets = xl(1) + [0.2 0.8]*width;
    tick_step = 10^floor(log10(width/3));
    tick_values = round(targets/tick_step)*tick_step;
    if tick_values(1) >= tick_values(2) || ...
            any(tick_values <= xl(1) | tick_values >= xl(2))
        tick_values = targets;
    end
    xticks(ax, tick_values);
    xticklabels(ax, arrayfun(@(v) sprintf('%.4g', v), ...
        tick_values, 'UniformOutput', false));
end

ylabel(tl, '\Delta log L', 'FontSize', 20);
xlabel(tl, 'log_{10}(Parameter value)', ...
    'Interpreter', 'tex', 'FontSize', 20);

drawnow;
export_figure(fig, 'tcaferm_profile_likelihood_real_log10');

%% ============================================================
%  14. Objective function
% =============================================================

function negLL = objective_function( ...
    p, ...
    params, ...
    data, ...
    tspan, ...
    y0, ...
    sigma_s, ...
    sigma_f)


    penalty = 1e12;


    try

        %% Parameter update

        params.enzyme_pk       = 10^p(1);
        params.enzyme_pdh      = 10^p(2);
        params.enzyme_cs       = 10^p(3);
        params.enzyme_aco      = 10^p(4);
        params.enzyme_id       = 10^p(5);
        params.enzyme_akd      = 10^p(6);
        params.enzyme_scs      = 10^p(7);
        params.enzyme_sd       = 10^p(8);
        params.enzyme_fumarase = 10^p(9);
        params.enzyme_md       = 10^p(10);
        params.enzyme_pta      = 10^p(11);
        params.enzyme_ack      = 10^p(12);

        params.k               = 10^p(13);
        params.alpha           = 10^p(14);


        %% ODE

        options = odeset( ...
            'RelTol', 1e-4, ...
            'AbsTol', 1e-6);


        [~, y_sim] = ...
            ode15s( ...
                @(t,y) function_fsolve(y, params), ...
                tspan, ...
                y0, ...
                options);


        %% Check simulation

        if isempty(y_sim)

            negLL = penalty;
            return;

        end


        model_output = ...
            y_sim(end, data.obs_idx)';


        if ...
                any(~isfinite(model_output)) || ...
                any(~isreal(model_output)) || ...
                any(model_output <= 0)

            negLL = penalty;
            return;

        end


        %% Flux

        flux = compute_flux( ...
            model_output', ...
            params);


        if ...
                isempty(flux) || ...
                any(~isfinite(flux)) || ...
                any(~isreal(flux))

            negLL = penalty;
            return;

        end


        flux_group_1 = ...
            flux(end-2:end-1);

        flux_group_2 = ...
            flux(3:end-3);


        if ...
                any(flux_group_1 <= 0) || ...
                any(flux_group_2 <= 0)

            negLL = penalty;
            return;

        end


        %% Residuals

        residuals = ...
            ( ...
                log10(data.values(1:end-2)) ...
                - ...
                log10(model_output(1:13)) ...
            ) ...
            ./ sigma_s;


        flux_residual1 = ...
            ( ...
                log10(data.values(end-1)) ...
                - ...
                log10(flux_group_1) ...
            ) ...
            ./ sigma_f;


        flux_residual2 = ...
            ( ...
                log10(data.values(end)) ...
                - ...
                log10(flux_group_2) ...
            ) ...
            ./ sigma_f;


        negLL = ...
            0.5 * ( ...
                sum(residuals(:).^2) ...
                + ...
                sum(flux_residual1(:).^2) ...
                + ...
                sum(flux_residual2(:).^2) ...
            );


        if ...
                ~isfinite(negLL) || ...
                ~isreal(negLL)

            negLL = penalty;

        end


    catch

        negLL = penalty;

    end

end


%% ============================================================
%  15. Insert free parameters
% =============================================================

function p_full = ...
    insert_params( ...
        p_fixed, ...
        p_free, ...
        idx_free)

    p_full = p_fixed;

    p_full(idx_free) = p_free;

end


%% ============================================================
%  16. Multi-start optimisation for one profile point
% =============================================================

function [best_p, best_fval, best_flag] = ...
    run_profile_multistart( ...
        fun_profile, ...
        starts, ...
        lb_free, ...
        ub_free, ...
        opts_profile)


    best_p = [];
    best_fval = Inf;
    best_flag = NaN;


    for s = 1:numel(starts)

        p0 = starts{s};


        if ...
                isempty(p0) || ...
                any(~isfinite(p0))

            continue;

        end


        p0 = max(p0, lb_free);
        p0 = min(p0, ub_free);


        try

            [p_opt, fval, exitflag] = ...
                fmincon( ...
                    fun_profile, ...
                    p0, ...
                    [], [], [], [], ...
                    lb_free, ...
                    ub_free, ...
                    [], ...
                    opts_profile);


            if ...
                    exitflag > 0 && ...
                    isfinite(fval) && ...
                    isreal(fval)

                if fval < best_fval

                    best_fval = fval;
                    best_p = p_opt;
                    best_flag = exitflag;

                end

            end


        catch

            % Ignore failed start

        end

    end


    if isempty(best_p)

        best_fval = NaN;
        best_flag = NaN;

    end

end