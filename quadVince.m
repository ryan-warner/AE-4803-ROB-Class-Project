n = 12;
m = 4;
dt = 0.01; % Idk

mass = 0.5                                                  % kg
inertiaMatrix = diag([0.0032, 0.0032, 0.0055])              % kgm^2
rotorTorqueConstant = 0.01691                               % 1/m
rotorMomentArm = 0.17                                       % m
gravity = 9.81                                              % m/s^2

x = sym('x', [n, 1]);
u = sym('u', [m, 1]);

assume(x, 'real');
assume(u, 'real');

% Phi       7
% Theta     8
% Psi       9

R = [cos(x(8)) * cos(x(9)), cos(x(8)) * sin(x(9)), -sin(x(8)); ...
    sin(x(8)) * sin(x(7)) * cos(x(9)) - cos(x(7)) * sin(x(9)), sin(x(8)) * sin(x(7)) * sin(x(9)) + cos(x(7)) * cos(x(9)), sin(x(7)) * cos(x(8)); ...
    sin(x(7)) * sin(x(9)) + cos(x(7)) * sin(x(8)) * cos(x(9)), sin(x(8)) * sin(x(9)) * cos(x(7)) - sin(x(7)) * sin(x(9)), cos(x(8)) * cos(x(7))];

x_ddot = (R * [0; 0; sum(u)] - [0; 0; mass*gravity]) / mass;

euler_matrix = [1, tan(x(8)) * sin(x(7)), tan(x(8)) * cos(x(7)); ...
    0, cos(x(7)), -sin(x(7)); ...
    0, sin(x(7)) / cos(x(8)), cos(x(7)) / cos(x(8))];

euler_dot = euler_matrix * [x(10); x(11); x(12)];

tempInertiaMatrix = num2cell(diag(inertiaMatrix));
[I_xx, I_yy, I_zz] = tempInertiaMatrix{:};

p_dot = (sqrt(2) / 2 * (u(1) + u(3) - u(2) - u(4)) * rotorMomentArm - (I_zz - I_yy) * x(11) * x(12)) / I_xx;
q_dot = (sqrt(2) / 2 * (u(3) + u(4) - u(1) - u(2)) * rotorMomentArm - (I_zz - I_xx) * x(10) * x(12)) / I_yy;
r_dot =  (rotorTorqueConstant * (u(1) + u(4) - u(2) - u(3))) / I_zz;

sym_dyn = x + dt * [x(4); x(5); x(6); x_ddot(1); x_ddot(2); x_ddot(3); euler_dot(1); euler_dot(2); euler_dot(3); p_dot; q_dot; r_dot];
dyn = matlabFunction(sym_dyn, Vars={x, u});

dyn_derivs = matlabFunction(jacobian(sym_dyn, x), ...
                            jacobian(sym_dyn, u), ...
                            reshape(jacobian(reshape(jacobian(sym_dyn, x), n * n, 1), x), n, n, n), ...
                            reshape(jacobian(reshape(jacobian(sym_dyn, u), n * m, 1), u), n, m, m), ...
                            reshape(jacobian(reshape(jacobian(sym_dyn, u), n * m, 1), x), n, m, n), ...
                            Vars={x, u}, Outputs={'fx', 'fu', 'fxx', 'fuu', 'fxu'});