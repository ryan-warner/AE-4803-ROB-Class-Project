clear

%% dynamics
% xdot = Ax+Bu
%A = [1 -5; 0 -1]; B = [0; 1]; % unstable
A = [0.5 1; -1 2]; B = [0; 1];
% satisfies f(0,0) = 0

%% constraints
% h = (x1-2)^2 + (x2-1.5)^2 - 1^2;  
% hx = [2(x1-2), 2(x2-1.5)];

%% unsafe control
Q = eye(2); R = 1;
Q(2,2) = 1;
K = lqr(A,B,Q,R);


%% simulate system

tspan = 0:0.01:10;
x0 = [4;4];
[T,X] = ode45(@(t,x) unsafe_ode(t,x,A,B,K),tspan, x0);
circ = [2-1 1.5-1.5 2*1 2*1];


figure(1)
plot(T,X,'LineWidth',2);
title('States vs time')

figure(2)
plot(X(:,1),X(:,2),'LineWidth',2)
title('X1 vs X2')
hold on
plot(x0(1),x0(2),'r*')
rectangle('Position', circ, 'Curvature', [1 1]);
hold off


%% Linear Embedded System
h0 = (0-2)^2 + (0-1.5)^2 - 1^2;
hx0 = [2*(0-2) 2*(0-1.5)];
gamma = 1;

Az1 = -1/h0^2*hx0*A - gamma * (1/h0^2*hx0);
Az2 = -gamma;

Bz = -1/h0^2*hx0*B;

Abar = [A, zeros(2,1); Az1 Az2];
Bbar = [B; Bz];

Qbas = 10;
Qbar = [Q zeros(2,1); zeros(1,2) Qbas];

K_safe = lqr(Abar,Bbar,Qbar,R);


%% simulate safe control

tspan = 0:0.01:10;
z0 = 1/(x0(1)-2)^2 + (x0(2)-1.5)^2 - 1 - 1/h0;
xbar0 = [x0;z0];

[T_safe,X] = ode45(@(t,x) safe_ode(t,x,A,B,K_safe,gamma,h0), tspan, xbar0);
circ = [2-1 1.5-1.5 2*1 2*1];

figure(3)
plot(T_safe,X,'LineWidth',2);
title('States vs time')

figure(4)
plot(X(:,1),X(:,2),'LineWidth',2)
title('X1 vs X2')
hold on
rectangle('Position', circ, 'Curvature', [1 1]);
hold off

%% ode's
function dx = unsafe_ode(t,x,A,B,K)
    u = -K*x;    
    dx = A*x + B*u;
end

function dx = safe_ode(t,x,A,B,K,gamma,h0)
    u = -K*x; % LQR
    dx = zeros(3,1);
    dx(1:2) = A*x(1:2) + B*u;
    
    h = (x(1)-2)^2 + (x(2) -1.5)^2 - 1^2;
    hx = [2*(x(1)-2) 2*(x(2)-1.5)];

    zdot = -1/h^2* hx * dx(1:2) - gamma*(x(3) + 1/h0 - 1/h);
    dx(3) = zdot;
end




