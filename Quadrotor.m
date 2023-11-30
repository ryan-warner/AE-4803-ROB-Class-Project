classdef Quadrotor
    % Instance variables, should probably do some property validation
    % Did not know that was a thing ngl
    properties
        bodyRates (1, 3) = [0, 0, 0]                                % rad/s
        eulerRates (1, 3) = [0, 0, 0]                               % rad/s
        thrust (1,4) = [0, 0, 0, 0]                                 % N?
        rotationMatrix (3, 3) = [0, 0, 0; 0, 0, 0; 0, 0, 0]         % Body to inertial
        mass = 0.5                                                  % kg
        inertiaMatrix (3, 3) = diag([0.0032, 0.0032, 0.0055])       % kgm^2
        rotorTorqueConstant = 0.01691                               % 1/m
        rotorMomentArm = 0.17                                       % m
        gravity = 9.81                                              % m/s^2

        % This representation of state may be wrong
        state (1, 12) = [0, 0, 0, 0, 0, 0] + eulerRates + bodyRates
    end
    
    % f(x, u) = x_t+1 = x_t + d_t * x_dot
    % f'(x, u) = x_ddot * dt


    % Methods (instance and static)
    methods
        function [forces, accelerations] = dynamics(obj, inputs)
            forces = obj.rotationMatrix * [0; 0; sum(inputs)] + [0; 0; -obj.mass * obj.gravity];
            accelerations = forces / obj.mass;
        end

        function [moments, accelerations] = rotationalDynamics(obj, inputs)
            tempBodyRates = num2cell(obj.bodyRates);
            [p, q, r] = tempBodyRates{:};

            tempInertiaMatrix = num2cell(diag(obj.inertiaMatrix));
            [I_xx, I_yy, I_zz] = tempInertiaMatrix{:};

            moments = [];
            
            moments(1) = sqrt(2) / 2 * (inputs(1) + inputs(3) - inputs(2) - inputs(4)) * obj.rotorMomentArm - (I_zz - I_yy) * q * r;
            moments(2) = sqrt(2) / 2 * (inputs(3) + inputs(4) - inputs(1) - inputs(2)) * obj.rotorMomentArm - (I_zz - I_xx) * p * r;
            moments(3) = obj.rotorTorqueConstant * (inputs(1) + inputs(4) - inputs(2) - inputs(3));

            accelerations = moments ./ [I_xx, I_yy, I_zz];
        end
        
        function obj = calcRotationMatrix(obj)
            tempEulerRates = num2cell(obj.eulerRates);
            [phi, theta, psi] = tempEulerRates{:};

            obj.rotationMatrix = [cos(theta) * cos(psi), cos(theta) * sin(psi), -sin(theta); ...
                sin(theta) * sin(phi) * cos(psi) - cos(phi) * sin(psi), sin(theta) * sin(phi) * sin(psi) + cos(phi) * cos(psi), sin(phi) * cos(theta); ...
                sin(phi) * sin(psi) + cos(phi) * sin(theta) * cos(psi), sin(theta) * sin(psi) * cos(phi) - sin(phi) * sin(psi), cos(theta) * cos(phi)];
        end

        function accelerations = calcEulerRates(obj)
            tempBodyRates = num2cell(obj.bodyRates);
            [p, q, r] = tempBodyRates{:};

            tempEulerRates = num2cell(obj.eulerRates);
            [phi, theta, ~] = tempEulerRates{:};

            accelerations = [1, tan(theta) * sin(phi), tan(theta) * cos(phi);...
                0, cos(phi), -sin(phi);...
                0, sin(phi) / cos(theta), cos(phi) / cos(theta)] * [p; q; r];
        end

    end
end