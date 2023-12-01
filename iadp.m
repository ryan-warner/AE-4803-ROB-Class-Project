function [ctl_params, traj_costs] = iadp(ic, init_ctl_params, dyn, dyn_derivs, cost, cost_derivs, term_cost, term_cost_derivs, options)
    arguments
        ic (:, 1) {mustBeNumeric}
        init_ctl_params (:, 1) struct
        dyn (1, 1) function_handle
        dyn_derivs (1, 1) function_handle
        cost (1, 1) function_handle
        cost_derivs (1, 1) function_handle
        term_cost (1, 1) function_handle
        term_cost_derivs (1, 1) function_handle

        options.StopTol (1, 1) {mustBeNonnegative} = 0
        options.MaxIters (1, 1) {mustBePositive, mustBeInteger} = 100
        % options.RegDamping (1, 1) {mustBeNonnegative} = 0.0;
        % options.RegFactor (1, 1) {mustBeNonnegative} = 2;
        options.Mode (1, 1) string {mustBeMember(options.Mode, ["iLQR", "DDP"])} = "DDP"
    end

    traj_costs = zeros(1, options.MaxIters);
    ctl_params = init_ctl_params;

    create_quad = @ddp_quad;

    for i = 1:options.MaxIters
        [states, controls, traj_costs(i)] = forward_pass(ic, dyn, cost, term_cost, ctl_params);
        ctl_params = backward_pass(states, controls, dyn_derivs, cost_derivs, term_cost_derivs, create_quad, 0.0);

        % if options.RegDamping > 0
        %     % Not the most efficient but works. Could reuse trajectory
        %     % rollout to make it faster.
        %     [~, ~, tc] = forward_pass(ic, dyn, cost, term_cost, ctl_params);
        %     reg_iter = 1;
        % 
        %     while tc > traj_costs(i)
        %         ctl_params = backward_pass(states, controls, dyn_derivs, cost_derivs, term_cost_derivs, create_quad, options.RegDamping * (options.RegFactor ^ reg_iter));
        %         [~, ~, tc] = forward_pass(ic, dyn, cost, term_cost, ctl_params);
        %         reg_iter = reg_iter + 1;
        %     end
        % end

        % if i > 1 && options.StopTol > 0 && stop_cond(traj_costs(i - 1), traj_costs(i), options.StopTol)
        %     traj_costs = traj_costs(1:i);
        %     break;
        % end
    end
end

function ctl_params = solve_quad(state, control, q_params)
    K = -q_params.Quu \ q_params.Qux;
    k = -q_params.Quu \ q_params.Qu;
    ctl_params = struct('K', K, 'd', control + k - K * state);
end

function q_params = ddp_quad(state, control, next_value_params, dyn_derivs, cost_derivs)
    [fx, fu, fxx, fuu, fxu] = dyn_derivs(state, control);
    [cx, cu, cxx, cuu, cxu] = cost_derivs(state, control);
    cux = cxu';
    n = 12;
    m = 4;

    % This is not a fast way to do this but it is less error prone when
    % encoding the equations. You can speed it up a lot using vectorized
    % operations

    Qxx = cxx + fx' * next_value_params.Vxx * fx;
    Quu = cuu + fu' * next_value_params.Vxx * fu;
    Qxu = cxu + fx' * next_value_params.Vxx * fu;

    for i = 1:n
        Qxx = Qxx + next_value_params.Vx(i) * reshape(fxx(i, :, :), n, n);
        Quu = Quu + next_value_params.Vx(i) * reshape(fuu(i, :, :), m, m);
        Qxu = Qxu + next_value_params.Vx(i) * reshape(fxu(i, :, :), n, m);
    end

    q_params = struct('Qx', cx + fx' * next_value_params.Vx, ...
                      'Qu', cu + fu' * next_value_params.Vx, ...
                      'Qxx', Qxx, ...
                      'Quu', Quu, ...
                      'Qxu', Qxu);

    q_params.Qux = q_params.Qxu';
end

function params = value_params(q_params)
    params = struct('Vx', q_params.Qx -  q_params.Qxu * (q_params.Quu \ q_params.Qu), ...
                    'Vxx', q_params.Qxx -  q_params.Qxu * (q_params.Quu \ q_params.Qux));
end

function ctl_params = backward_pass(states, controls, dyn_derivs, cost_derivs, term_cost_derivs, create_quad, regularizer)
    horizon = size(controls, 2);
    n = size(states, 1);
    m = size(controls, 1);
    ctl_params = repmat(struct('K', zeros(m, n), 'd', zeros(m, 1)), horizon, 1);
    
    [Vx, Vxx] = term_cost_derivs(states(:, end));
    valfun_params = repmat(struct('Vx', Vx, 'Vxx', Vxx), horizon + 1, 1);
    qfun_params = repmat(struct('Qx', zeros(n, 1), 'Qu', zeros(m, 1), 'Qxx', zeros(n, n), ...
        'Quu', zeros(m, m), 'Qux', zeros(m, n), 'Qxu', zeros(n, m)), horizon, 1);
    
    for t = horizon:-1:1
        qfun_params(t) = create_quad(states(:, t), controls(:, t), valfun_params(t + 1), dyn_derivs, cost_derivs);
        %qfun_params(t).Quu = qfun_params(t).Quu; %+ regularizer * eye(m);
        valfun_params(t) = value_params(qfun_params(t));
        ctl_params(t) = solve_quad(states(:, t), controls(:, t), qfun_params(t));
    end
    1;
end

function stop = stop_cond(prev_traj_cost, current_traj_cost, tol)
    stop = (prev_traj_cost - current_traj_cost) < tol;
end
