function [X] = fnsimulate(xo,u_new,Horizon,dt,sigma)

X = xo;

m = 0.5; I_xx = 0.0032; I_yy = 0.0032; I_zz = 0.0055;
k_t = 0.01691; l = 0.17; g = 9.81;

    for k = 1:(Horizon-1)
        
        % states
        x = X(1,k); y = X(2,k); z = X(3,k);
        x_dot = X(4,k); y_dot = X(5,k); z_dot = X(6,k);
        phi = X(7,k); theta = X(8,k); psi = X(9,k);
        p = X(10,k); q = X(11,k); r = X(12,k);
        
        % inputs
        f1 = u_new(1,k);
        f2 = u_new(2,k);
        f3 = u_new(3,k);
        f4 = u_new(4,k);

        Fx = [x_dot; 
              y_dot; 
              z_dot; 
              -(sin(theta)*(f1 + f2 + f3 + f4))/m; 
              (cos(theta)*sin(phi)*(f1 + f2 + f3 + f4))/m; 
              (f1*cos(phi)*cos(theta) - g*m + f2*cos(phi)*cos(theta) + f3*cos(phi)*cos(theta) + f4*cos(phi)*cos(theta))/m; 
              p + r*cos(phi)*tan(theta) + q*sin(phi)*tan(theta); 
              q*cos(phi) - r*sin(phi); 
              (r*cos(phi))/cos(theta) + (q*sin(phi))/cos(theta); 
              (q*r*(I_yy - I_zz) + (sqrt(2)*l*(f1 - f2 + f3 - f4))/2)/I_xx; 
              -(p*r*(I_xx - I_zz) + (sqrt(2)*l*(f1 + f2 - f3 - f4))/2)/I_yy; 
              (k_t*(f1 - f2 - f3 + f4))/I_zz];
     

        X(:,k+1) = X(:,k) + Fx*dt;% + G_x*u_new(:,k)*dt + G_x*u_new(:,k)*sqrt(dt)*sigma*randn;
    
    end

end