function negLL = objective_function( ...
    p, ...
    params, ...
    data, ...
    tspan, ...
    y0, ...
    sigma_s, ...
    sigma_f)


    % ---------------------------------------------------------
    % Large objective value returned when simulation enters
    % an invalid numerical region.
    % ---------------------------------------------------------

    penalty = 1e12;


    try

        %% ----------------------------------------------------
        % Update parameters
        % -----------------------------------------------------

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


        %% ----------------------------------------------------
        % ODE solution
        % -----------------------------------------------------

        options = odeset( ...
            'RelTol', 1e-4, ...
            'AbsTol', 1e-6);


        [~, y_sim] = ...
            ode15s( ...
                @(t,y) function_fsolve(y, params), ...
                tspan, ...
                y0, ...
                options);


        %% ----------------------------------------------------
        % Check ODE result
        % -----------------------------------------------------

        if isempty(y_sim)

            negLL = penalty;
            return;

        end


        model_output = ...
            y_sim(end, data.obs_idx)';


        % log10 concentrations require:
        %
        % concentration > 0
        %
        if ...
                any(~isfinite(model_output)) || ...
                any(~isreal(model_output)) || ...
                any(model_output <= 0)

            negLL = penalty;
            return;

        end


        %% ----------------------------------------------------
        % Flux calculation
        % -----------------------------------------------------

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


        %% ----------------------------------------------------
        % These are exactly the same flux subsets that were
        % used in your original objective function.
        % -----------------------------------------------------

        flux_group_1 = ...
            flux(end-2:end-1);

        flux_group_2 = ...
            flux(3:end-3);


        % Because you take log10 below,
        % all selected fluxes must be positive.

        if ...
                any(flux_group_1 <= 0) || ...
                any(flux_group_2 <= 0)

            negLL = penalty;
            return;

        end


        %% ----------------------------------------------------
        % Concentration residuals
        % -----------------------------------------------------

        residuals = ...
            ( ...
                log10(data.values(1:end-2)) ...
                - ...
                log10(model_output(1:13)) ...
            ) ...
            ./ sigma_s;


        %% ----------------------------------------------------
        % Flux residuals
        %
        % IMPORTANT:
        % I preserved your original definition here.
        % -----------------------------------------------------

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


        %% ----------------------------------------------------
        % Negative log likelihood
        % -----------------------------------------------------

        negLL = ...
            0.5 * ( ...
                sum(residuals(:).^2) ...
                + ...
                sum(flux_residual1(:).^2) ...
                + ...
                sum(flux_residual2(:).^2) ...
            );


        %% ----------------------------------------------------
        % Final numerical safety check
        % -----------------------------------------------------

        if ...
                ~isfinite(negLL) || ...
                ~isreal(negLL)

            negLL = penalty;

        end


    catch

        % Any ODE / flux / numerical failure becomes
        % a bad objective value rather than killing fmincon.

        negLL = penalty;

    end

end