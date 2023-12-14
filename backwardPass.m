function gains = backwardPass(states, inputs, derivatives, options)
    hover_thrust = 0.5 * 9.81 / 4;
    gains = repmat(struct('K', ones(options.m, options.n), 'k', ones(options.m, 1), 'optimal_control', hover_thrust * ones(options.m, 1)), options.horizon, 1);
   
    [V_x, V_xx] = derivatives.cost_final(states(:, end));
    backwardsUpdateTerms.V_x = V_x;
    backwardsUpdateTerms.V_xx = V_xx;
    
    for i = options.horizon:-1:1
        q_terms = calcQTerms(derivatives, states(:, i), inputs(:, i), backwardsUpdateTerms, options);
        % if any(isnan(q_terms.q_uu), 'all')
        %     error('q_terms.q_uu has NaN values at iteration %d', i);
        % end
        backwardsUpdateTerms = backwardsUpdate(q_terms);
        deltaState = states(:,i+1) - states(:,i);   % ADDED BY ME
        gains(i) = calcGains(q_terms, deltaState, inputs(:, i));  % replaced states(:, i) with deltaState
    end
end

function q_terms = calcQTerms(derivatives, state, input, backwardsUpdateTerms, options)
    derivativeValues = calcDerivatives(state, input, derivatives); % Replaces struct w real values, might be ok

    q_terms.q_x = derivativeValues.l_x + derivativeValues.f_x.' * backwardsUpdateTerms.V_x; % OK
    q_terms.q_u = derivativeValues.l_u + derivativeValues.f_u.' * backwardsUpdateTerms.V_x; % OK
    q_terms.q_xx = derivativeValues.l_xx + derivativeValues.f_x.' * backwardsUpdateTerms.V_xx * derivativeValues.f_x; % OK
    q_terms.q_uu = derivativeValues.l_uu + derivativeValues.f_u.' * backwardsUpdateTerms.V_xx * derivativeValues.f_u; % OK
    q_terms.q_xu = derivativeValues.l_xu + derivativeValues.f_x.' * backwardsUpdateTerms.V_xx * derivativeValues.f_u; % OK

    for i = 1:options.n
        q_terms.q_xx = q_terms.q_xx + backwardsUpdateTerms.V_x(i) * reshape(derivativeValues.f_xx(i, :, :), options.n, options.n);
        q_terms.q_uu = q_terms.q_uu + backwardsUpdateTerms.V_x(i) * reshape(derivativeValues.f_uu(i, :, :), options.m, options.m);
        q_terms.q_xu = q_terms.q_xu + backwardsUpdateTerms.V_x(i) * reshape(derivativeValues.f_xu(i, :, :), options.n, options.m);
    end

    % Add regularization
    %q_terms.q_uu = q_terms.q_uu + regularizationFactor * eye(m);
    regularizationFactor= 0.0001;
    eig_tol = 1e3;
    while min(eig(q_terms.q_uu + regularizationFactor * eye(options.m))) <= eig_tol
        regularizationFactor = regularizationFactor * 2;
    end

    q_terms.q_uu = q_terms.q_uu + regularizationFactor * eye(options.m);

    q_terms.q_ux = q_terms.q_xu'; % OK
end

function backwardsUpdateTerms = backwardsUpdate(q_terms)
    backwardsUpdateTerms.V_x = q_terms.q_x -  q_terms.q_xu * (q_terms.q_uu \ q_terms.q_u);
    backwardsUpdateTerms.V_xx = q_terms.q_xx -  q_terms.q_xu * (q_terms.q_uu \ q_terms.q_ux);
end

function gains = calcGains(q_terms, state, control)
    gains.K = -pinv(q_terms.q_uu) * q_terms.q_ux;
    gains.k = -pinv(q_terms.q_uu) * q_terms.q_u;
    % gains.optimal_control = control + gains.k - gains.K * state;    %
    % changed to plus K from minus K, not sure why this is recalculated
    gains.optimal_control = 0.5 * 9.81 / 4;     % hardcoded whoops
end