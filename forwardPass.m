function [states, inputs, costs] = forwardPass(initial, dynamics, cost, finalCost, gains)
    states = [initial zeros(12, 800)];
    inputs = zeros(4, 800);
    costs = zeros(1, 800 + 1);

    for i = 1:800
        iterGains = gains(:, t);
        inputs(:, i) =  iterGains.K * states(:, t) + iterGains.k;
        costs(i) = cost(states(:, i), inputs(:, i));
        states(:, i + 1) = dynamics(states(:, i), inputs(:, i));
    end

    costs(end) = finalCost(states(:, end));
    costs = sum(costs);
end