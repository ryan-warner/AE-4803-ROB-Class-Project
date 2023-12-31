function [currentGains, iterationCosts] = IADP_FINAL(initial, dynamics, costs, derivatives, options)
    iterationCosts = zeros([1, options.maxIterations]);
    hover_thrust = 0.5 * 9.81 / 4;
    % currentGains = repmat(struct('K', randn(options.m, options.n), 'k', randn(options.m, 1), 'optimal_control', hover_thrust * ones(options.m, 1)), options.horizon, 1);
    currentGains = options.currentGains;
    
    for i = 1:options.maxIterations
        [iterationCosts(i), tempGains] = step(initial, dynamics, costs, currentGains, derivatives, options);
        noiseLevel = 100;        
        if i > 1 && iterationCosts(i) > iterationCosts(i - 1)
            fprintf('Cost increased, taking last iteration\n');
            iterationCosts(i) = iterationCosts(i - 1);
            currentGains = addNoiseToGains(tempGains, noiseLevel);

        elseif iterationCosts(i) == 0 && i > 1
            fprintf('Cost is 0, taking last iteration\n');
            iterationCosts(i) = iterationCosts(i - 1);
            currentGains = addNoiseToGains(tempGains, noiseLevel);
        else
            currentGains = tempGains;
        end

        fprintf('Iteration: %i | Cost: %.08f\n', i, iterationCosts(i));
    end

    function [cost, gains, states] = step(initial, dynamics, costs, currentGains, derivatives, options)
        [states, inputs, cost] = forwardPass(initial, dynamics, costs, currentGains, options);
        gains = backwardPass(states, inputs, derivatives, options);
    end    
end

function gains = addNoiseToGains(gains, noiseLevel)
    for j = 1:length(gains)
        gains(j).K = gains(j).K + noiseLevel * randn(size(gains(j).K));
        gains(j).k = gains(j).k + noiseLevel * randn(size(gains(j).k));
        gains(j).optimal_control = gains(j).optimal_control + noiseLevel * randn(size(gains(j).optimal_control));
    end
end