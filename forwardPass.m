function [states, inputs, trajectoryCostsTotal] = forwardPass(initial, dynamics, costs, gains, options)
    states = [initial zeros(options.n, options.horizon)];
    inputs = zeros(options.m, options.horizon);
    trajectoryCosts = zeros(1, options.horizon + 1);

    for i = 1:options.horizon
        iterGains = gains(i);
        inputs(:, i) =  max(iterGains.K * states(:, i) + iterGains.optimal_control, [0;0;0;0]);
        trajectoryCosts(i) = costs.cost(states(:, i), inputs(:, i));
        states(:, i + 1) = dynamics(states(:, i), inputs(:, i)) * options.timestep + states(:, i);
    end

    trajectoryCosts(end) = costs.finalCost(states(:, end));
    trajectoryCostsTotal = sum(trajectoryCosts);
    if isnan(trajectoryCostsTotal)
        disp('Final cost is NaN');
    end
end