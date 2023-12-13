function [states, inputs, trajectoryCosts] = forwardPass(initial, dynamics, costs, gains, options)
    states = [initial zeros(options.n, options.horizon)];
    inputs = zeros(options.m, options.horizon);
    trajectoryCosts = zeros(1, options.horizon + 1);

    for i = 1:options.horizon
        iterGains = gains(i);
        inputs(:, i) =  iterGains.K * states(:, i) + iterGains.optimal_control;
        trajectoryCosts(i) = costs.cost(states(:, i), inputs(:, i));
        states(:, i + 1) = dynamics(states(:, i), inputs(:, i));
    end

    trajectoryCosts(end) = costs.finalCost(states(:, end));
    trajectoryCosts = sum(trajectoryCosts);
end