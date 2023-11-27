classdef Quadrotor
    % Instance variables, should probably do some property validation
    % Did not know that was a thing ngl
    properties
        bodyRates (1, 3) = [0, 0, 0]                                % rad/s
        eulerRates (1, 3) = [0, 0, 0]                               % rad/s
        thrust (1,4) = [0, 0, 0, 0]                                 % N?
        rotationMatrix (3, 3) = [0, 0, 0; 0, 0, 0; 0, 0, 0]         % Body to inertial
        mass = 0.5                                                  % kg
        inertiaMatrix (3, 3) = diag(0.0032, 0.0032, 0.0055)         % kgm^2
        rotorTorqueConstant = 0.01691                               % 1/m
        rotorMomentArm = 0.17                                       % m
        gravity = 9.81                                              % m/s^2

        % This representation of state may be wrong
        state (1, 12) = [0, 0, 0, 0, 0, 0] + eulerRates + bodyRates
    end
    
    % Methods (instance and static)
    methods

    end
end