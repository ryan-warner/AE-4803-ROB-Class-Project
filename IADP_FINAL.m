function [currentGains, iterationCosts] = IADP_FINAL(initial, dynamics, costs, derivatives)
    maxIterations = 20;
    iterationCosts = zeros([1, maxIterations]);
    n = 12;
    m = 4;
    horizon = 800;
    hover_thrust = 0.5 * 9.81 / 4;
    currentGains = repmat(struct('K', randn(m, n), 'k', randn(m, 1), 'optimal_control', hover_thrust * randn(m, 1)), horizon, 1);

    for i = 1:maxIterations
        [iterationCosts(i), currentGains] = step(initial, dynamics, costs, currentGains, derivatives);
        
        fprintf('Iteration: %i | Cost: %.02f\n', i, iterationCosts(i));
    end

    % Hardcoding a bunch of simulation vars for rn
    function [cost, gains, states] = step(initial, dynamics, costs, currentGains, derivatives)
        [states, inputs, cost] = forwardPass(initial, dynamics, costs, currentGains);
        gains = backwardPass(states, inputs, derivatives);
    end    
end