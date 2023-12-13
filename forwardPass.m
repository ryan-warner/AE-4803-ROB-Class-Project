function [states, inputs, trajectoryCosts] = forwardPass(initial, dynamics, costs, gains)
    states = [initial zeros(12, 800)];
    inputs = zeros(4, 800);
    trajectoryCosts = zeros(1, 800 + 1);

    for i = 1:800
        iterGains = gains(i);
        inputs(:, i) =  iterGains.K * states(:, i) + iterGains.optimal_control;
        trajectoryCosts(i) = costs.cost(states(:, i), inputs(:, i));
        states(:, i + 1) = dynamics(states(:, i), inputs(:, i));
    end

    trajectoryCosts(end) = costs.finalCost(states(:, end));
    trajectoryCosts = sum(trajectoryCosts);
end