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

Q = 0.0 * diag([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
R = 1.0 * diag([1.0, 1.0, 1.0, 1.0]);
Qf = 10.0 * diag([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
goal_state = [5; 3; 2; 0; 0; 0; 0; 0; 0; 0; 0; 0];

cost = @(state, control) 0.5 * (state - goal_state)' * Q * (state - goal_state) + 0.5 * control' * R * control;
cost_derivs = @(state, control) deal(Q * (state - goal_state), R * control, Q, R, 0);

term_cost = @(state) 0.5 * (state - goal_state)' * Qf * (state - goal_state);
term_cost_derivs = @(state) deal(Qf * (state - goal_state), Qf);

ic = [-3; -2; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0];
horizon = 800;
init_ctl_params = repmat(struct('K', zeros(m, n), 'd', zeros(m, 1)), horizon, 1);

[ctl_params, ddp_traj_costs] = iadp(ic, init_ctl_params, dyn, dyn_derivs, cost, cost_derivs, term_cost, term_cost_derivs,  MaxIters=10, Mode="DDP");

[ddp_states, ~, ~] = forward_pass(ic, dyn, cost, term_cost, ctl_params);

disp("Done :)")

%% Quadcopter Recursive Model Predictive Control (rMPC)

%% Obstacle Avoidance with Discrete Barrier States