function [currentGains, iterationCosts] = IADP_FINAL(initial, dynamics, costs, derivatives)
    maxIterations = 3;
    iterationCosts = zeros([1, maxIterations]);
    n = 12;
    m = 4;
    horizon = 800;
    currentGains = repmat(struct('K', zeros(m, n), 'k', zeros(m, 1)), horizon, 1);

    for i = 1:maxIterations
        [iterationCosts(i), currentGains] = step(initial, dynamics, costs, currentGains, derivatives);
    end

    % Hardcoding a bunch of simulation vars for rn
    function [cost, gains] = step(initial, dynamics, costs, currentGains, derivatives)
        [states, inputs, cost] = forwardPass(initial, dynamics, costs, currentGains);
        gains = backwardPass(states, inputs, derivatives);
    end    
end