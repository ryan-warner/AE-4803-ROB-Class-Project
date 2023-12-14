function [A,B] = fnState_And_Control_Transition_Matrices(h10,h20,h30,h1_x_0,h2_x_0,h3_x_0, gamma1, X,u)

    m = 0.5; I_xx = 0.0032; I_yy = 0.0032; I_zz = 0.0055;
    k_t = 0.01691; l = 0.17; g = 9.81;
    
    x = X(1,1); y = X(2,1); z = X(3,1);
    x_dot = X(4,1); y_dot = X(5,1); z_dot = X(6,1);
    phi = X(7,1); theta = X(8,1); psi = X(9,1);
    p = X(10,1); q = X(11,1); r = X(12,1);
    
    f1 = u(1,1);
    f2 = u(2,1);
    f3 = u(3,1);
    f4 = u(4,1);  
    
    A = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0; 
         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 
         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0; 
         0, 0, 0, 0, 0, 0, 0, -(cos(theta)*(f1 + f2 + f3 + f4))/m, 0, 0, 0, 0; 
         0, 0, 0, 0, 0, 0, (cos(phi)*cos(theta)*(f1 + f2 + f3 + f4))/m, -(sin(phi)*sin(theta)*(f1 + f2 + f3 + f4))/m, 0, 0, 0, 0; 
         0, 0, 0, 0, 0, 0, -(f1*cos(theta)*sin(phi) + f2*cos(theta)*sin(phi) + f3*cos(theta)*sin(phi) + f4*cos(theta)*sin(phi))/m, -(f1*cos(phi)*sin(theta) + f2*cos(phi)*sin(theta) + f3*cos(phi)*sin(theta) + f4*cos(phi)*sin(theta))/m, 0, 0, 0, 0; 
         0, 0, 0, 0, 0, 0, q*cos(phi)*tan(theta) - r*sin(phi)*tan(theta), r*cos(phi)*(tan(theta)^2 + 1) + q*sin(phi)*(tan(theta)^2 + 1), 0, 1, sin(phi)*tan(theta), cos(phi)*tan(theta); 
         0, 0, 0, 0, 0, 0, - r*cos(phi) - q*sin(phi), 0, 0, 0, cos(phi), -sin(phi); 
         0, 0, 0, 0, 0, 0, (q*cos(phi))/cos(theta) - (r*sin(phi))/cos(theta), (r*cos(phi)*sin(theta))/cos(theta)^2 + (q*sin(phi)*sin(theta))/cos(theta)^2, 0, 0, sin(phi)/cos(theta), cos(phi)/cos(theta); 
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (r*(I_yy - I_zz))/I_xx, (q*(I_yy - I_zz))/I_xx; 
         0, 0, 0, 0, 0, 0, 0, 0, 0, -(r*(I_xx - I_zz))/I_yy, 0, -(p*(I_xx - I_zz))/I_yy; 
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];   
    
    B = [0, 0, 0, 0; 
         0, 0, 0, 0; 
         0, 0, 0, 0; 
         -sin(theta)/m, -sin(theta)/m, -sin(theta)/m, -sin(theta)/m; 
         (cos(theta)*sin(phi))/m, (cos(theta)*sin(phi))/m, (cos(theta)*sin(phi))/m, (cos(theta)*sin(phi))/m; 
         (cos(phi)*cos(theta))/m, (cos(phi)*cos(theta))/m, (cos(phi)*cos(theta))/m, (cos(phi)*cos(theta))/m; 
         0, 0, 0, 0; 
         0, 0, 0, 0; 
         0, 0, 0, 0; 
         (sqrt(2)*l)/(2*I_xx),  -(sqrt(2)*l)/(2*I_xx), (sqrt(2)*l)/(2*I_xx), -(sqrt(2)*l)/(2*I_xx); 
         -(sqrt(2)*l)/(2*I_yy), -(sqrt(2)*l)/(2*I_yy), (sqrt(2)*l)/(2*I_yy), (sqrt(2)*l)/(2*I_yy); 
         k_t/I_zz, -k_t/I_zz, -k_t/I_zz, k_t/I_zz];

    % barrier functions at current step
    % h1_k = (x_k-2.2)^2 + (y_k-2.2)^2 + (z_k-1)^2   - 1;
    % h2_k = (x_k)^2     + (y_k+0.2)^2 + (z_k)^2     - 1;
    % h3_k = (x_k-3)^2   + (y_k)^2     + (z_k-0.5)^2 - 1;
    
    XXX = (1/h10^2)*h1_x_0 + (1/h20^2)*h2_x_0 + (1/h30^2)*h3_x_0;
    Az1 = -(XXX)*A - gamma1*(XXX);
    Az2 = -gamma1;

    Bz = -XXX*B;
    
    Abar = double([A, zeros(12,1); Az1 Az2]);
    Bbar = double([B; Bz]);

    A = Abar;
    B = Bbar;

% doesn't work
%     h = [(x-2.2)^2 + (y-2.2)^2 + (z-1)^2   - 1;
%          (x)^2     + (y+0.2)^2 + (z)^2     - 1;
%          (x-3)^2   + (y)^2     + (z-0.5)^2 - 1];
% 
%     h0 = [(0-2.2)^2 + (0-2.2)^2 + (0-1)^2   - 1;
%           (0)^2     + (0+0.2)^2 + (0)^2     - 1;
%           (0-3)^2   + (0)^2     + (0-0.5)^2 - 1];
%     
%     h0sq = h0.'*h0;
%     h_x0 = [h1_x_0; h2_x_0; h3_x_0];
% 
%     Az1 = -1/h0sq * h_x0 * A - gamma1 * 1/h0sq*h_x0;
%     Az2 = -gamma1 * eye(3);
% 
%     Bz = -1/h0sq * h_x0 * B;
% 
%     Abar = double([A, zeros(12,3); Az1 Az2]);
%     Bbar = double([B; Bz]); 
%     
%     A = Abar;
%     B = Bbar;



end