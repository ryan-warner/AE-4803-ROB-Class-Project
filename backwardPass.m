function gains = backwardPass(states, inputs, derivatives)
    gains = repmat(struct('K', zeros(m, n), 'k', zeros(m, 1)), horizon, 1);
    
    [Vx, Vxx] = term_cost_derivs(states(:, end));
    
    for i = 800:-1:1
        q_terms = calcQTerms()
        backwardsUpdateTerms = backwardsUpdate(q_terms, gains(i))
        gains(i) = calcGains(q_terms)
    end
end

%% 

% My derivatives might need to be evaluated... think I left them in the HW
% as already evaluated, so this is a bit problematic. Can rework tho, not a
% huge deal at all.

%%

function q_terms = calcQTerms(derivatives, backwardsUpdateTerms)
    q_terms.q_x = derivatives.l_x + derivatives.f_x.' * backwardsUpdateTerms.V_x; % OK
    q_terms.q_u = derivatives.l_u + derivatives.f_u.' * backwardsUpdateTerms.V_x; % OK
    q_terms.q_xx = derivatives.l_xx + derivatives.f_x.' * backwardsUpdateTerms.V_xx * derivatives.f_x; % OK
    q_terms.q_uu = derivatives.l_uu + derivatives.f_u.' * backwardsUpdateTerms.V_xx * derivatives.f_u; % OK
    q_terms.q_xu = derivatives.l_xu + derivatives.f_x.' * backwardsUpdateTerms.V_xx * derivatives.f_u; % OK

    % He does other things to sort out the other matrices...?

    q_terms.q_ux = q_terms.q_xu; % OK
end

function backwardsUpdateTerms = backwardsUpdate(q_terms, gains)
        backwardsUpdateTerms.V_x = q_terms.q_x + gains.K.' * q_terms.q_uu * gains.k + gains.K.' * q_terms.q_u + q_terms.q_ux.' * gains.k;
        backwardsUpdateTerms.V_xx = q_terms.q_xx + gains.K.' * q_terms.q_uu * gains.K + gains.K.' * q_terms.q_ux + q_terms.q_ux.' * gains.K;
    end

function gains = calcGains(q_terms)
    gains.K = -inv(q_terms.q_uu) * q_terms.q_ux;
    gains.k = -inv(q_terms.q_uu) * q_terms.q_u;
end