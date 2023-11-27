%% General Setup
close all, clear all, clc
box on

%Simulation Parameters
simulation_vars.max_iterations = 20;
simulation_vars.time_step = 0.05;
simulation_vars.time_horizon = 200;

% Symbolics
syms x [3 1]
syms x_g [3 1]
syms u [1 1]
syms dt
vehicle_syms.position = x;
vehicle_syms.position_goal = x_g;
vehicle_syms.input = u;
vehicle_syms.dt = dt;

syms r
syms Q [3 3]
syms Q_tf [3 3]
syms R [1 1]
syms p [2 1]
cost_syms.position = Q;
cost_syms.input = R;
cost_syms.position_final = Q_tf;
cost_syms.obstacle_position = p;
cost_syms.car_radius = r;

%% Question 1
% Case 1
% Initial Location 0 0 0
% Goal Location 4 5 pi/2
% Q 0.1 * diag([1 1 0])
% R 1.0
% Q_tf = 0

% Setup
initial_location = [0 0 0].';
target_location = [4 5 pi/2].';
vehicle_velocity = 1;
has_obstacles = false;
Q = 0.1 * diag([1 1 0]);
R = 1.0;
Q_tf = 0.0 * diag([1 1 0]);

% Car Instantiation
car = DubinsCar(vehicle_velocity, initial_location, vehicle_syms, cost_syms, simulation_vars, has_obstacles, [], 0);
car = car.setCostMatrices(Q, R, Q_tf, target_location);

% IADP
[feedback_gains_ddp, feedforward_gains_ddp, iteration_costs_ddp, x_trajectory_ddp, ~] = car.IADP(true);
[feedback_gains_ilqr, feedforward_gains_ilqr, iteration_costs_ilqr, x_trajectory_ilqr, ~] = car.IADP(false);

% Simulation
%[x_trajectory_ddp, ~, ~] = car.simulate(feedback_gains_ddp, feedforward_gains_ddp);
%[x_trajectory_ilqr, ~, ~] = car.simulate(feedback_gains_ilqr, feedforward_gains_ilqr);

fig1 = figure(Name="Case 1 - Cost");
hold on, box on
plot(1:car.simulation_vars.max_iterations + 1, iteration_costs_ddp)
plot(1:car.simulation_vars.max_iterations + 1, iteration_costs_ilqr)
legend('DDP Cost', 'iLQR Cost', Location='northwest')

fig2 = figure(Name="Case 1 - Trajectory");
hold on, box on
plot(target_location(1), target_location(2), Marker="*", Color='r')
plot(x_trajectory_ddp(1, :), x_trajectory_ddp(2, :))
plot(x_trajectory_ilqr(1, :), x_trajectory_ilqr(2, :))
axis equal
legend("Target Location", 'DDP Trajectory', 'iLQR Trajectory', Location='northwest')

saveImage(fig1)
saveImage(fig2)

% Case 2
% Initial Location 0 0 0
% Goal Location 4 5 pi/2
% Q 0.0
% R 1.0
% Q_tf = 10.0 * diag([1 1 0])

% Setup
initial_location = [0 0 0].';
target_location = [4 5 pi/2].';
vehicle_velocity = 1;
Q = 0.0 * diag([1 1 0]);
R = 1.0;
Q_tf = 10.0 * diag([1 1 0]);

% Car Instantiation
car = DubinsCar(vehicle_velocity, initial_location, vehicle_syms, cost_syms, simulation_vars, has_obstacles, [], 0);
car = car.setCostMatrices(Q, R, Q_tf, target_location);

% IADP
[feedback_gains_ddp, feedforward_gains_ddp, iteration_costs_ddp, x_trajectory_ddp, ~] = car.IADP(true);
[feedback_gains_ilqr, feedforward_gains_ilqr, iteration_costs_ilqr, x_trajectory_ilqr, ~] = car.IADP(false);

% Simulation
%[x_trajectory_ddp, ~, ~] = car.forwardsPass(feedback_gains_ddp, feedforward_gains_ddp);
%[x_trajectory_ilqr, ~, ~] = car.forwardsPass(feedback_gains_ilqr, feedforward_gains_ilqr);

fig1 = figure(Name="Case 2 - Cost");
hold on, box on
plot(1:car.simulation_vars.max_iterations + 1, iteration_costs_ddp)
plot(1:car.simulation_vars.max_iterations + 1, iteration_costs_ilqr)
legend('DDP Cost', 'iLQR Cost', Location='northwest')

fig2 = figure(Name="Case 2 - Trajectory");
hold on, box on
plot(target_location(1), target_location(2), Marker="*", Color='r')
plot(x_trajectory_ddp(1, :), x_trajectory_ddp(2, :))
plot(x_trajectory_ilqr(1, :), x_trajectory_ilqr(2, :))
axis equal
legend("Target Location", 'DDP Trajectory', 'iLQR Trajectory', Location='northwest')

saveImage(fig1)
saveImage(fig2)

%% Question 2
% Case 1
% Initial Location 0 0 -pi/4
% Goal Location 4 5 pi/2
% Obtacles (4, 1) (2, 6)
% r 1
% Q 0.01 * diag([1 1 0])
% R 3.0
% Q_tf = diag([1 1 0])

% Setup
initial_location = [0 0 -pi/4].';
target_location = [4 5 pi/2].';
vehicle_velocity = 1;
has_obstacles = true;
Q = 0.01 * diag([1 1 0]);
R = 3.0;
Q_tf = 1.0 * diag([1 1 0]);
obstacles = [[4 1]; [2 6]];
r = 1;

% Car Instantiation
car = DubinsCar(vehicle_velocity, initial_location, vehicle_syms, cost_syms, simulation_vars, has_obstacles, obstacles, r);
car = car.setCostMatrices(Q, R, Q_tf, target_location);

% IADP
[feedback_gains_ddp, feedforward_gains_ddp, iteration_costs_ddp, x_trajectory_ddp, ~] = car.IADP(true);
[feedback_gains_ilqr, feedforward_gains_ilqr, iteration_costs_ilqr, x_trajectory_ilqr, ~] = car.IADP(false);

% Simulation
%[x_trajectory_ddp, ~, ~] = car.forwardsPass(feedback_gains_ddp, feedforward_gains_ddp);
%[x_trajectory_ilqr, ~, ~] = car.forwardsPass(feedback_gains_ilqr, feedforward_gains_ilqr);

fig1 = figure(Name="Q2 Case 1 - Cost");
hold on, box on
plot(1:car.simulation_vars.max_iterations + 1, iteration_costs_ddp)
plot(1:car.simulation_vars.max_iterations + 1, iteration_costs_ilqr)
legend('DDP Cost', 'iLQR Cost', Location='northwest')

fig2 = figure(Name="Q2 Case 1 - Trajectory");
hold on, box on
plot(target_location(1), target_location(2), Marker="*", Color='r')
plot(x_trajectory_ddp(1, :), x_trajectory_ddp(2, :))
plot(x_trajectory_ilqr(1, :), x_trajectory_ilqr(2, :))
ylim([-1, 7])
axis equal
for obstacle = obstacles
    viscircles(obstacle', r)
end
legend("Target Location", 'DDP Trajectory', 'iLQR Trajectory', Location='northwest')


saveImage(fig1)
saveImage(fig2)

% Case 2
% Initial Location 0 0 pi/4
% Goal Location 4 5 pi/2
% Obtacles (4, 1) (2, 6)
% r 1
% Q 0.01 * diag([1 1 0])
% R 3.0
% Q_tf = diag([1 1 0])

% Setup
initial_location = [0 0 pi/4].';
target_location = [4 5 pi/2].';
vehicle_velocity = 1;
has_obstacles = true;
Q = 0.01 * diag([1 1 0]);
R = 3.0;
Q_tf = 1.0 * diag([1 1 0]);
obstacles = [[4 1]; [2 6]];
r = 1;

% Car Instantiation
car = DubinsCar(vehicle_velocity, initial_location, vehicle_syms, cost_syms, simulation_vars, has_obstacles, obstacles, r);
car = car.setCostMatrices(Q, R, Q_tf, target_location);

% IADP
[feedback_gains_ddp, feedforward_gains_ddp, iteration_costs_ddp, x_trajectory_ddp, ~] = car.IADP(true);
[feedback_gains_ilqr, feedforward_gains_ilqr, iteration_costs_ilqr, x_trajectory_ilqr, ~] = car.IADP(false);

% Simulation
%[x_trajectory_ddp, ~, ~] = car.forwardsPass(feedback_gains_ddp, feedforward_gains_ddp);
%[x_trajectory_ilqr, ~, ~] = car.forwardsPass(feedback_gains_ilqr, feedforward_gains_ilqr);

fig1 = figure(Name="Q2 Case 2 - Cost");
hold on, box on
plot(1:car.simulation_vars.max_iterations + 1, iteration_costs_ddp)
plot(1:car.simulation_vars.max_iterations + 1, iteration_costs_ilqr)
legend('DDP Cost', 'iLQR Cost', Location='northwest')

fig2 = figure(Name="Q2 Case 2 - Trajectory");
hold on, box on
plot(target_location(1), target_location(2), Marker="*", Color='r')
plot(x_trajectory_ddp(1, :), x_trajectory_ddp(2, :))
plot(x_trajectory_ilqr(1, :), x_trajectory_ilqr(2, :))
ylim([-1, 7])
axis equal
for obstacle = obstacles
    viscircles(obstacle', r)
end
legend("Target Location", 'DDP Trajectory', 'iLQR Trajectory', Location='northwest')


saveImage(fig1)
saveImage(fig2)

% Case 3
% Initial Location 0 0 pi/2
% Goal Location 4 5 pi/2
% Obtacles (4, 1) (2, 6)
% r 1
% Q 0.01 * diag([1 1 0])
% R 3.0
% Q_tf = diag([1 1 0])

% Setup
initial_location = [0 0 pi/2].';
target_location = [4 5 pi/2].';
vehicle_velocity = 1;
has_obstacles = true;
Q = 0.01 * diag([1 1 0]);
R = 3.0;
Q_tf = 1.0 * diag([1 1 0]);
obstacles = [[4 1]; [2 6]];
r = 1;

% Car Instantiation
car = DubinsCar(vehicle_velocity, initial_location, vehicle_syms, cost_syms, simulation_vars, has_obstacles, obstacles, r);
car = car.setCostMatrices(Q, R, Q_tf, target_location);

% IADP
[feedback_gains_ddp, feedforward_gains_ddp, iteration_costs_ddp, x_trajectory_ddp, ~] = car.IADP(true);
[feedback_gains_ilqr, feedforward_gains_ilqr, iteration_costs_ilqr, x_trajectory_ilqr, ~] = car.IADP(false);

% Simulation
%[x_trajectory_ddp, ~, ~] = car.forwardsPass(feedback_gains_ddp, feedforward_gains_ddp);
%[x_trajectory_ilqr, ~, ~] = car.forwardsPass(feedback_gains_ilqr, feedforward_gains_ilqr);

fig1 = figure(Name="Q2 Case 3 - Cost");
hold on, box on
plot(1:car.simulation_vars.max_iterations + 1, iteration_costs_ddp)
plot(1:car.simulation_vars.max_iterations + 1, iteration_costs_ilqr)
legend('DDP Cost', 'iLQR Cost', Location='northwest')

fig2 = figure(Name="Q2 Case 3 - Trajectory");
hold on, box on
plot(target_location(1), target_location(2), Marker="*", Color='r')
plot(x_trajectory_ddp(1, :), x_trajectory_ddp(2, :))
plot(x_trajectory_ilqr(1, :), x_trajectory_ilqr(2, :))
ylim([-1, 7])
axis equal
for obstacle = obstacles
    viscircles(obstacle', r)
end
legend("Target Location", 'DDP Trajectory', 'iLQR Trajectory', Location='northwest')


saveImage(fig1)
saveImage(fig2)

close all