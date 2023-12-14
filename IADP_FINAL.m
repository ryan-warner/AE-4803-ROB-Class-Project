function [currentGains, iterationCosts] = IADP_FINAL(initial, dynamics, costs, derivatives, options)
    iterationCosts = zeros([1, options.maxIterations]);
    hover_thrust = 0.5 * 9.81 / 4;
    currentGains = repmat(struct('K', ones(options.m, options.n), 'k', ones(options.m, 1), 'optimal_control', hover_thrust * ones(options.m, 1)), options.horizon, 1);

    for i = 1:options.maxIterations
        [iterationCosts(i), currentGains] = step(initial, dynamics, costs, currentGains, derivatives, options);
        
        fprintf('Iteration: %i | Cost: %.02f\n', i, iterationCosts(i));
    end

    function [cost, gains, states] = step(initial, dynamics, costs, currentGains, derivatives, options)
        [states, inputs, cost] = forwardPass(initial, dynamics, costs, currentGains, options);
        gains = backwardPass(states, inputs, derivatives, options);
    end    
end