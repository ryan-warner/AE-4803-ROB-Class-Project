%% Initialization
clear all, close all, clc

pos_init = [-3 -2 -1];
pos_final = [5 3 2];

tf = 8;
time_step = 0.01;

horizon = tf / time_step;


warning('off','MATLAB:illConditionedMatrix')

%% Quadrotor DDP

% Think we'll generally have pretty high Q_tf and Q Matrices, low R Matrix.

quadVince;
hover_thrust = mass * gravity / 4;

Q = 1.0 * diag([0.1, 0.1, 100, 1.0, 1.0, 350, 1.0, 1.0, 1.0, 70, 70, 70]);
R = .01 * diag([1.0, 1.0, 1.0, 1.0]);
Qf = 1.0 * diag([100, 100, 100, 100, 100, 1000, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]);

%Q = 1.0 * diag([1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
%R = 10000000.0 * diag([1, 1, 1, 1]);
%Qf = 1000.0 * diag([1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]);

ic = [-3; -2; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0];
goal_state = [5; 3; 2; 0; 0; 0; 0; 0; 0; 0; 0; 0];

cost = @(state, control) 0.5 * (state - goal_state)' * Q * (state - goal_state) + 0.5 * (control - hover_thrust)' * R * (control - hover_thrust);
cost_derivs = @(state, control) deal(Q * (state - goal_state), R * (control - hover_thrust), Q, R, 0);

cost_final = @(state) 0.5 * (state - goal_state)' * Qf * (state - goal_state);
cost_final_derivs = @(state) deal(Qf * (state - goal_state), Qf);

derivatives.dynamics = dynamics_derivatives;
derivatives.cost = cost_derivs;
derivatives.cost_final = cost_final_derivs;

costs.cost = cost;
costs.finalCost = cost_final;

% Simulation Options
options.timestep = 0.01;
options.tf = 8;
options.horizon = options.tf / options.timestep;
options.n = size(x, 1);
options.m = size(u, 1);
options.maxIterations = 20;

[gains, iadp_costs, prevStates] = IADP_FINAL(ic, dynamics, costs, derivatives, options);

[ddp_states, ddp_inputs, ~] = forwardPass(ic, dynamics, costs, gains, options, prevStates);

disp("Done :)")
time(1)=0;
for i= 2:options.horizon
    time(i) =time(i-1) + options.timestep;  
end

[test, test2] = quadrotor_visualize(ddp_states, ddp_inputs, time);

%plot3(ddp_states(1, :), ddp_states(2, :), ddp_states(3, :))

%% Quadcopter Recursive Model Predictive Control (rMPC)

%% Obstacle Avoidance with Discrete Barrier States