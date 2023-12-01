function [states, inputs, traj_cost] = forward_pass(ic, dyn, cost, term_cost, ctl_params)
    horizon = length(ctl_params);
    n = length(ic);
    m = length(controller(ic, ctl_params(1)));
    states = [ic zeros(n, horizon)];
    inputs = zeros(m, horizon);
    costs = zeros(1, horizon + 1);

    for t = 1:length(ctl_params)
        inputs(:, t) = controller(states(:, t), ctl_params(t));
        costs(t) = cost(states(:, t), inputs(:, t));
        states(:, t + 1) = dyn(states(:, t), inputs(:, t));
    end

    costs(end) = 
    traj_cost = sum(costs);
end