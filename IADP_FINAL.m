function [gains, iterationCosts] = IADP_FINAL(initial, dynamics, costFn, finalCost)
    iterationCosts = zeros([1, maxIterations]);
    currentGains.feedbackGains = zeros([12, 800]);
    currentGains.feedbackGains = zeros(4, 800);

    for i = 1:maxIterations
        [iterationCosts(i), currentGains] = step(initial, dynamics, costFn, finalCost, currentGains);
    end

    % Hardcoding a bunch of simulation vars for rn
    function [cost, gains] = step(initial, dynamics, costFn, finalCost, currentGains)
        [states, inputs, cost] = forwardPass(initial, dynamics, costFn, finalCost, currentGains)
        gains = backwardPass(states, inputs, );


        %gains.feedbackGains = feedbackGains;
        %gains.feedforwardGains = feedforwardGains;
    end    
end