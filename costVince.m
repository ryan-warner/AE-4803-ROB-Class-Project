function [cost, cost_derivs, term_cost, term_cost_derivs] = quadcost(Q, R, Qf, goal_state)
    [n, ~] = size(Q);
    [m, ~] = size(R);

    cost = @(state, control) 0.5 * (state - goal_state)' * Q * (state - goal_state) + 0.5 * control' * R * control;
    cost_derivs = @(state, control) deal(Q * (state - goal_state), R * control, Q, R, zeros(n, m));
    
    term_cost = @(state) 0.5 * (state - goal_state)' * Qf * (state - goal_state);
    term_cost_derivs = @(state) deal(Qf * (state - goal_state), Qf);
end