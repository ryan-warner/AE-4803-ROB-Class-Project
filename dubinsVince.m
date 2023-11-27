n = 3;
m = 1;
dt = 0.05;

x = sym('x', [n, 1]);
u = sym('u', [m, 1]);

assume(x, 'real');
assume(u, 'real');

sym_dyn = x + dt * [cos(x(3, :)); sin(x(3, :)); u(1, :)];
dyn = matlabFunction(sym_dyn, Vars={x, u});

dyn_derivs = matlabFunction(jacobian(sym_dyn, x), ...
                            jacobian(sym_dyn, u), ...
                            reshape(jacobian(reshape(jacobian(sym_dyn, x), n * n, 1), x), n, n, n), ...
                            reshape(jacobian(reshape(jacobian(sym_dyn, u), n * m, 1), u), n, m, m), ...
                            reshape(jacobian(reshape(jacobian(sym_dyn, u), n * m, 1), x), n, m, n), ...
                            Vars={x, u}, Outputs={'fx', 'fu', 'fxx', 'fuu', 'fxu'});