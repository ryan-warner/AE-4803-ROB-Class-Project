function gains = backwardPass(states, inputs, derivatives, options)
    hover_thrust = 0.5 * 9.81 / 4;
    gains = repmat(struct('K', randn(options.m, options.n), 'k', randn(options.m, 1), 'optimal_control', hover_thrust * randn(options.m, 1)), options.horizon, 1);
   
    [V_x, V_xx] = derivatives.cost_final(states(:, end));
    backwardsUpdateTerms.V_x = V_x;
    backwardsUpdateTerms.V_xx = V_xx;
    
    for i = 800:-1:1
        q_terms = calcQTerms(derivatives, states(:, i), inputs(:, i), backwardsUpdateTerms);
        % if any(isnan(q_terms.q_uu), 'all')
        %     error('q_terms.q_uu has NaN values at iteration %d', i);
        % end
        backwardsUpdateTerms = backwardsUpdate(q_terms);
        gains(i) = calcGains(q_terms, states(:, i), inputs(:, i));
    end
end

function q_terms = calcQTerms(derivatives, state, input, backwardsUpdateTerms)
    derivatives = calcDerivatives(state, input, derivatives); % Replaces struct w real values, might be ok

    q_terms.q_x = derivatives.l_x + derivatives.f_x.' * backwardsUpdateTerms.V_x; % OK
    q_terms.q_u = derivatives.l_u + derivatives.f_u.' * backwardsUpdateTerms.V_x; % OK
    q_terms.q_xx = derivatives.l_xx + derivatives.f_x.' * backwardsUpdateTerms.V_xx * derivatives.f_x; % OK
    q_terms.q_uu = derivatives.l_uu + derivatives.f_u.' * backwardsUpdateTerms.V_xx * derivatives.f_u; % OK
    q_terms.q_xu = derivatives.l_xu + derivatives.f_x.' * backwardsUpdateTerms.V_xx * derivatives.f_u; % OK

    n = 12;
    m = 4;
    for i = 1:n
        q_terms.q_xx = q_terms.q_xx + backwardsUpdateTerms.V_x(i) * reshape(derivatives.f_xx(i, :, :), n, n);
        q_terms.q_uu = q_terms.q_uu + backwardsUpdateTerms.V_x(i) * reshape(derivatives.f_uu(i, :, :), m, m);
        q_terms.q_xu = q_terms.q_xu + backwardsUpdateTerms.V_x(i) * reshape(derivatives.f_xu(i, :, :), n, m);
    end

    % Add regularization
    %q_terms.q_uu = q_terms.q_uu + regularizationFactor * eye(m);
    regularizationFactor= 0.0001;
    eig_tol = 4;
    while min(eig(q_terms.q_uu + regularizationFactor * eye(m))) <= eig_tol
        regularizationFactor = regularizationFactor * 2;
    end

    q_terms.q_uu = q_terms.q_uu + regularizationFactor * eye(m);

    q_terms.q_ux = q_terms.q_xu'; % OK
end

function backwardsUpdateTerms = backwardsUpdate(q_terms)
    backwardsUpdateTerms.V_x = q_terms.q_x -  q_terms.q_xu * (q_terms.q_uu \ q_terms.q_u);
    backwardsUpdateTerms.V_xx = q_terms.q_xx -  q_terms.q_xu * (q_terms.q_uu \ q_terms.q_ux);
end

function gains = calcGains(q_terms, state, control)
    gains.K = -pinv(q_terms.q_uu) * q_terms.q_ux;
    gains.k = -pinv(q_terms.q_uu) * q_terms.q_u;
    gains.optimal_control = control + gains.k - gains.K * state;
end