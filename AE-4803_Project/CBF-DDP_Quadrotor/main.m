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

numStates = 13;
numInputs = 4;

% constants
trim = 1.2263; % N
m = 0.5; I_xx = 0.0032; I_yy = 0.0032; I_zz = 0.0055;
k_t = 0.01691; l = 0.17; grav = 9.81;


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
Q(1,1)   = .1;  % x
Q(2,2)   = .1;  % y
Q(3,3)   = 100; % z
Q(4,4)   = 1;   % x_dot
Q(5,5)   = 1;   % y_dot
Q(6,6)   = 350;  % z_dot
Q(7,7)   = 1;   % phi
Q(8,8)   = 1;   % theta
Q(9,9)   = 1;   % psi
Q(10,10) = 70;  % p
Q(11,11) = 70;  % q
Q(12,12) = 70;  % r
Q(13,13) = 1000;   % barrier state


% Weight in Final State:
Q_f = 10*eye(numStates);
Q_f(1,1) = 2000; % x
Q_f(2,2) = 2000; % y
Q_f(3,3) = 2000; % z
Q_f(4,4)   = 100;   % x_dot
Q_f(5,5)   = 110;   % y_dot
Q_f(6,6)   = 1000;  % z_dot
Q_f(13,13) = 100;  % barrier state


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
gamma = 0.3;

    %%%% Barrier State $$$$
    syms x y z 
    syms dummy [9 1] % dummy vars
    X = [x;y;z; dummy];

    h1(x,y,z) = (x-2.2)^2 + (y-2.2)^2 + (z-1)^2 - 1;
    h2(x,y,z) = (x)^2 + (y+0.2)^2 + (z)^2 - 1;
    h3(x,y,z) = (x-3)^2 + (y)^2 + (z-0.5)^2 - 1;
    
    h1_x = jacobian(h1,X);
    h2_x = jacobian(h2,X);
    h3_x = jacobian(h3,X);

    h10 = double(h1(0,0,0));
    h20 = double(h2(0,0,0));
    h30 = double(h3(0,0,0));

    h1_x_0 = double(h1_x(0,0,0));
    h2_x_0 = double(h2_x(0,0,0));
    h3_x_0 = double(h3_x(0,0,0));
    
    % "noise" for system
    gamma1 = 0.7;

    % Run DDP algorithm
    for k = 1:num_iter
    
    %%%%%%%%%%%%%%%%%%%% Linearization of the dynamics %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Quadratic Approximations of the cost function %%%%%%%%%%%%
        for  j = 1:(Horizon-1)
            
            [l0,l_x,l_xx,l_u,l_uu,l_ux] = fnCost(x_traj(:,j), p_target, u_k(:,j), numStates, numInputs, Q, R);
            q0(j)      = dt * l0;
            q_k(:,j)   = dt * l_x;
            Q_k(:,:,j) = dt * l_xx;
            r_k(:,j)   = dt * l_u;
            R_k(:,:,j) = dt * l_uu;
            P_k(:,:,j) = dt * l_ux; 
            
            [dfx,dfu] = fnState_And_Control_Transition_Matrices(h10,h20,h30,h1_x_0,h2_x_0,h3_x_0,gamma1, x_traj(:,j),u_k(:,j));
           
            A(:,:,j) = eye(numStates) + dfx*dt;
            B(:,:,j) = dfu*dt;  

        end
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Find the controls %%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vxx(:,:,Horizon) = Q_f;
    Vx(:,Horizon)    = Q_f * (x_traj(:,Horizon) - p_target); 
    V(Horizon)       = 0.5 * (x_traj(:,Horizon) - p_target)' * Q_f * (x_traj(:,Horizon) - p_target); 
    
    
    %%%%%%%%%%%%%%%% Backpropagation of the Value Function %%%%%%%%%%%%%%%%
        for j = (Horizon-1):-1:1
             
            H = R_k(:,:,j) + B(:,:,j)' * Vxx(:,:,j+1) * B(:,:,j);
            G = P_k(:,:,j) + B(:,:,j)' * Vxx(:,:,j+1) * A(:,:,j);   
            g = r_k(:,j)   + B(:,:,j)' * Vx(:,j+1);

            % Feedback
            L_k(:,:,j) = - H \ G;
            % Feedforward
            l_k (:,j)  = - H \ g;  
        
            % Gradients
            Vxx(:,:,j) = Q_k(:,:,j) + A(:,:,j)'*Vxx(:,:,j+1)*A(:,:,j) + L_k(:,:,j)'*H*L_k(:,:,j) + L_k(:,:,j)'*G + G'*L_k(:,:,j);
            Vx(:,j)    = q_k(:,j) + A(:,:,j)'*Vx(:,j+1) + L_k(:,:,j)'*g + G'*l_k(:,j) + L_k(:,:,j)'*H*l_k(:,j);
            V(:,j)     = q0(j) + V(j+1) + 0.5*l_k (:,j)'*H*l_k(:,j) + l_k(:,j)'*g;
        end
    
    %%%%%%%%%%%%%%%%%%%%%% Forward Propagation: %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%% Find the controls/ forward %%%%%%%%%%%%%%%%%%%%%%
    dx = zeros(numStates,1);
        for i = 1:(Horizon-1)
           du         = l_k(:,i)    + L_k(:,:,i)*dx;
           dx         = A(:,:,i)*dx + B(:,:,i)*du;
           u_new(:,i) = u_k(:,i)    + gamma*du;

           % prevent control from going negative
           u_new(:,i) = max(u_new(:,i), [0;0;0;0]);
        end
    
    u_k = u_new;
    


    %%%%%%%%%%%%%%%%% Simulation of the Nonlinear System %%%%%%%%%%%%%%%%%%
    [x_traj]    = fnsimulate(gamma1,x0,u_k,Horizon,dt);
    [Cost(:,k)] = fnCostComputation(x_traj,u_k,p_target,dt,Q_f,R);
    
    %fprintf('DDP Iteration %d,  Current Cost = %e \n',k,Cost(1,k));
    
    end %% end iterating over the DDP algorithm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
exportgraphics(f,'Displacements.jpg')


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
exportgraphics(f,'Velocities.jpg')


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
exportgraphics(f,'Euler Angles.jpg')


% Body rates
f=figure();
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
exportgraphics(f,'Body Rates.jpg')

% Cost
f=figure();
plot(Cost,'r','linewidth',4); 
xlabel('Iterations','fontsize',20)
ylabel('Cost magnitude','fontsize',20)
title('Cost','fontsize',20);
f.Position=[0,0,1000,1500];
exportgraphics(f,'Cost.jpg')

% DBaS
f = figure()
plot(time,x_traj(13,:),'linewidth',4); hold on
title('Discrete Barrier State, $w_k$','fontsize',20,'interpreter','latex')
xlabel('Time (s)','fontsize',20)
ylabel('Magnitude')
grid; hold off;
f.Position=[0,0,1000,1500];
exportgraphics(f,'Control Trajectory.jpg')

% Control Trajectories
f=figure()
plot(u_k','linewidth',4); 
xlabel('Time (s)','fontsize',20)
ylabel('Thrust input (N)','fontsize',20)
title('Control','fontsize',20);
legend('$f_1$','$f_2$','$f_3$','$f_4$','interpreter','latex','fontsize',20)
f.Position=[0,0,1000,1500];
exportgraphics(f,'Control Trajectory.jpg')
   