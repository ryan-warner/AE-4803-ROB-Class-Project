%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Differential Dynamic Programming  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Course: Robotics and Autonomy     %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  AE8803  Fall  2018                %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Author: Evangelos Theodorou       %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Modified by: Ashkar I. Awal       %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  AE 4803  Fall  2022               %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Quadcopter DDP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

numStates = 12;
numInputs = 4;
trim = 1.2263; % N

% Discretization
dt = 0.01;  % s
% total time allowed
t = 8; % s

% Horizon
Horizon = t/dt; % 800
% Number of Iterations
num_iter = 100;

% Running Cost Weights
Q = eye(numStates);
Q(1,1)   = 6;   % x
Q(2,2)   = 6;   % y
Q(3,3)   = 150; % z
Q(4,4)   = 2;   % x_dot
Q(5,5)   = 2;   % y_dot
Q(6,6)   = 10;  % z_dot
Q(7,7)   = 1;   % phi
Q(8,8)   = 1;   % theta
Q(9,9)   = 1;   % psi
Q(10,10) = 70;  % p
Q(11,11) = 70;  % q
Q(12,12) = 70;  % r


% Weight in Final State:
Q_f = 1*eye(numStates);
Q_f(1,1) = 500; % x
Q_f(2,2) = 500; % y
Q_f(3,3) = 100; % z
Q(4,4)   = 10;   % x_dot
Q(5,5)   = 10;   % y_dot
Q(6,6)   = 50;  % z_dot


% Weight in the Control:
R = 0.01*eye(numInputs);


% Initial Configuration:
x0 = zeros(numStates,1);
x0(1,1) = -3; % x
x0(2,1) = -2; % y
x0(3,1) = -1; % z

% Initial Control:
u_k = trim*ones(numInputs,Horizon-1);
du_k = zeros(numInputs,Horizon-1);


% Initial trajectory:
x_traj = zeros(numStates,Horizon);


% Target: 
p_target = zeros(numStates,1);
p_target(1,1) = 5;  % x
p_target(2,1) = 3;  % y
p_target(3,1) = 2;  % z


% Learning Rate:
gamma = 0.2;

%
current_state = x0;
stateArray = [];
iterDDP = 200;
Horizon = 75;
iterMPC = 400;

    for kk = 1:iterMPC
        % Run DDP algorithm
        [x_traj, u_k, L_k, l_k, Cost] = fnDDP(iterDDP, numStates, numInputs, Horizon, dt, Q, Q_f, R, current_state(:,kk), u_k, du_k, x_traj, p_target, gamma);
        
        % calcualte next step using u(1)
        [next_state] = fnDynamics(current_state(:,kk), u_k(:,1), Horizon,dt,0);
        
        % set current state = next state and loop again
        x_traj = next_state.*ones(numStates,Horizon);
        current_state = [current_state,next_state];
        
        % slide inputs
        u_k = u_k(:,2:end);
        u_k(:,end+1) = u_k(:,end);
    
    end

      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPC plotting vars
x_traj = current_state; x_traj(:,1) = [];
Horizon = iterMPC;

% time array   
time(1) = 0;
for i = 2:Horizon
    time(i) = time(i-1) + dt;
end


% Positions
f = figure();
  subplot(3,1,1)
plot(time,x_traj(1,:),'linewidth',4); hold on;  
plot(time,p_target(1,1)*ones(1,Horizon),'red','linewidth',4)
title('Linear displacement, $x$','fontsize',20,'interpreter','latex'); 
xlabel('Time (s)','fontsize',20)
ylabel('x-Displacement (m)')
grid; hold off;

  subplot(3,1,2); 
plot(time,x_traj(2,:),'linewidth',4); hold on;
plot(time,p_target(2,1)*ones(1,Horizon),'red','linewidth',4)
title('Linear displacement, $y$','fontsize',20,'interpreter','latex');
xlabel('Time (s)','fontsize',20)
ylabel('y-displacement (m)')
grid; hold off;

  subplot(3,1,3); 
plot(time,x_traj(3,:),'linewidth',4); hold on
plot(time,p_target(3,1)*ones(1,Horizon),'red','linewidth',4)
title('Linear displacement, $z$','fontsize',20,'interpreter','latex')
xlabel('Time (s)','fontsize',20)
ylabel('z-displacement (m)')
grid; hold off;
f.Position=[0,0,1000,1500];
exportgraphics(f,'PosMPC.jpg')

% Velocities
f = figure();
  subplot(3,1,1)
plot(time,x_traj(4,:),'linewidth',4); hold on;  
plot(time,p_target(4,1)*ones(1,Horizon),'red','linewidth',4)
title('Velocity, $\dot{x}$','fontsize',20,'interpreter','latex'); 
xlabel('Time (s)','fontsize',20)
ylabel('x-Velocity (m/s)')
grid; hold off;

  subplot(3,1,2); 
plot(time,x_traj(5,:),'linewidth',4); hold on;
plot(time,p_target(5,1)*ones(1,Horizon),'red','linewidth',4)
title('Velocity, $\dot{y}$','fontsize',20,'interpreter','latex');
xlabel('Time (s)','fontsize',20)
ylabel('y-Velocity (m/s)')
grid; hold off;

  subplot(3,1,3); 
plot(time,x_traj(6,:),'linewidth',4); hold on
plot(time,p_target(6,1)*ones(1,Horizon),'red','linewidth',4)
title('Velocity, $\dot{z}$','fontsize',20,'interpreter','latex')
xlabel('Time (s)','fontsize',20)
ylabel('z-Velocity (m/s)')
grid; hold off;
f.Position=[0,0,1000,1500];
exportgraphics(f,'VelsMPC.jpg')


% Euler angles
f = figure();
  subplot(3,1,1)
plot(time,x_traj(7,:),'linewidth',4); hold on;
plot(time,p_target(7,1)*ones(1,Horizon),'red','linewidth',4)
title('Roll, $\phi$','fontsize',20,'interpreter','latex'); 
xlabel('Time (s)','fontsize',20)
ylabel('Roll (rad)')
grid; hold off;

  subplot(3,1,2); 
plot(time,x_traj(8,:),'linewidth',4); hold on;
plot(time,p_target(8,1)*ones(1,Horizon),'red','linewidth',4)
title('Pitch, $\theta$','fontsize',20,'interpreter','latex');
xlabel('Time (s)','fontsize',20)
ylabel('Pitch (rad)')
grid; hold off;

  subplot(3,1,3); 
plot(time,x_traj(9,:),'linewidth',4); hold on
plot(time,p_target(9,1)*ones(1,Horizon),'red','linewidth',4)
title('Yaw, $\psi$','fontsize',20,'interpreter','latex')
xlabel('Time (s)','fontsize',20)
ylabel('Yaw (rad)')
grid; hold off;
f.Position=[0,0,1000,1500];
exportgraphics(f,'EulerAnglesMPC.jpg')


% Body rates
f = figure();
  subplot(3,1,1)
plot(time,x_traj(10,:),'linewidth',4); hold on;
plot(time,p_target(10,1)*ones(1,Horizon),'red','linewidth',4)
title('Roll Rate, $r$','fontsize',20,'interpreter','latex'); 
xlabel('Time (s)','fontsize',20)
ylabel('Roll Rate (rad/s)')
grid; hold off;

  subplot(3,1,2); 
plot(time,x_traj(11,:),'linewidth',4); hold on;
plot(time,p_target(11,1)*ones(1,Horizon),'red','linewidth',4)
title('Pitch rate, $q$','fontsize',20,'interpreter','latex');
xlabel('Time (s)','fontsize',20)
ylabel('Pitch rate (rad/s)')
grid; hold off;

  subplot(3,1,3); 
plot(time,x_traj(12,:),'linewidth',4); hold on
plot(time,p_target(12,1)*ones(1,Horizon),'red','linewidth',4)
title('Yaw rate, $r$','fontsize',20,'interpreter','latex')
xlabel('Time (s)','fontsize',20)
ylabel('Yaw rate (rad/s)')
grid; hold off;
f.Position=[0,0,1000,1500];
exportgraphics(f,'BodyRatesMPC.jpg')

% Cost
f = figure()
plot(Cost,'r','linewidth',4); 
xlabel('Iterations','fontsize',20)
ylabel('Cost magnitude','fontsize',20)
title('Cost','fontsize',20);
f.Position=[0,0,1000,1500];
exportgraphics(f,'Cost.jpg')

% Control Trajectories
f = figure()
plot(u_k','linewidth',4); 
xlabel('Time (s)','fontsize',20)
ylabel('Thrust input (N)','fontsize',20)
title('Control','fontsize',20);
legend('$f_1$','$f_2$','$f_3$','$f_4$','interpreter','latex','fontsize',20)
f.Position=[0,0,1000,1500];
exportgraphics(f,'ControlMPC.jpg')