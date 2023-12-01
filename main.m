%% Initialization
clear all, close all, clc

pos_init = [-3 -2 -1];
pos_final = [5 3 2];

tf = 8;
time_step = 0.01;

horizon = tf / time_step;

%% Quadrotor DDP

% Think we'll generally have pretty high Q_tf and Q Matrices, low R Matrix.

quadVince;

Q = 1.0 * diag([1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
R = 10.0 * diag([1.0, 1.0, 1.0, 1.0]);
Qf = 100.0 * diag([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
goal_state = [5; 3; 2; 0; 0; 0; 0; 0; 0; 0; 0; 0];

cost = @(state, control) 0.5 * (state - goal_state)' * Q * (state - goal_state) + 0.5 * control' * R * control;
cost_derivs = @(state, control) deal(Q * (state - goal_state), R * control, Q, R, 0);

cost_final = @(state) 0.5 * (state - goal_state)' * Qf * (state - goal_state);
cost_final_derivs = @(state) deal(Qf * (state - goal_state), Qf);

derivatives.dynamics = dyn_derivs;
derivatives.cost = cost_derivs;
derivatives.cost_final = cost_final_derivs;

costs.cost = cost;
costs.finalCost = cost_final;

ic = [-3; -2; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0];

horizon = 800;
init_ctl_params = repmat(struct('K', zeros(m, n), 'k', zeros(m, 1)), horizon, 1);

[gains, iadp_costs] = IADP_FINAL(ic, dyn, costs, derivatives)


% dyn_derivs = function handle returning f_x, f_u, f_xx, f_uu, f_ux
% cost_derivs = funciton handle returning l_x. l_u, l_xx, l_uu, l_ux
% term_cost_derivs = function handle returning l_x_f, l_xx_f

% Maybe we matlab function all the derivatives and return a struct...?
% Could

%[ctl_params, ddp_traj_costs] = iadp(ic, init_ctl_params, dyn, dyn_derivs, cost, cost_derivs, term_cost, term_cost_derivs,  MaxIters=5, Mode="DDP");
[ddp_states, ~, ~] = forwardPass(ic, dyn, costs, gains);

disp("Done :)")

plot3(ddp_states(1, :), ddp_states(2, :), ddp_states(3, :))

%% Quadcopter Recursive Model Predictive Control (rMPC)

%% Obstacle Avoidance with Discrete Barrier States