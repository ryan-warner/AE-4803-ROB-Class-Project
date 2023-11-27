classdef DubinsCar
    properties
        velocity                    % Vehicle total velocity
        state                       % Current vehicle state (y, z, theta)
        initial_state               % Initial vehicle state (y, z, theta)
        state_vars                  % Symbolic references to state
        cost_vars                   % Symbolic references to cost
        cost                        % Symbolic representation of cost fn
        dynamics                    % Symbolic representation of dynamics
        true_symbolic_derivatives   % Derivatives of cost fn - symbolic
        symbolic_derivatives        % Derivatives of cost fn - Q/R Subbed
        derivatives                 % Derivatives of cost fn - Numeric
        q_terms                     % Bellman eqn. Q terms
        gains                       % Control gains (K - Feedback, k - Feedforward)
        backwards_update_terms      % Backwards update cost terms (H_x, H_xx)
        simulation_vars             % Simulation variables (time horizon/step)
    end
    methods
        % Constructor
        function obj = DubinsCar(velocity, state, vehicle_syms, cost_syms, simulation_vars, has_obstacles, obstacles, r)
            obj.velocity = velocity;            % scalar
            obj.state = state;                  % 3x1 (Column vector)
            obj.initial_state = state;          % 3x1 (Column vector)
            
            obj.state_vars.position = vehicle_syms.position;  % 3x1 (Column vector)
            obj.state_vars.input = vehicle_syms.input;        % scalar 
            obj.state_vars.dt = vehicle_syms.dt;              % scalar
            obj.state_vars.position_goal = vehicle_syms.position_goal;  % 3x1 (Column vector)

            obj.cost_vars.position = cost_syms.position;                % 3x3 Matrix
            obj.cost_vars.position_final = cost_syms.position_final;    % 3x3 Matrix
            obj.cost_vars.input = cost_syms.input;                      % 1x1 Matrix -> Scalar

            % Set up dynamics
            obj.dynamics = obj.state_vars.position + [cos(obj.state_vars.position(3)) * obj.velocity; ...
                sin(obj.state_vars.position(3)) * obj.velocity; ...
                obj.state_vars.input] .* obj.state_vars.dt;

            % Set up costs - ALL SCALAR OUTPUT
            obj.cost.state = 0.5 * (obj.state_vars.position - obj.state_vars.position_goal)' * cost_syms.position * (obj.state_vars.position - obj.state_vars.position_goal);
            obj.cost.state_final = 0.5 * (obj.state_vars.position - obj.state_vars.position_goal)' * cost_syms.position_final * (obj.state_vars.position - obj.state_vars.position_goal);
            obj.cost.input = 0.5 * obj.state_vars.input' * cost_syms.input * obj.state_vars.input;
            obj.cost.obstacles = exp(-((obj.state_vars.position(1) - cost_syms.obstacle_position(1))^2 + (obj.state_vars.position(2) - cost_syms.obstacle_position(2))^2)/(2*cost_syms.car_radius));
        
            obj.cost.total = obj.cost.state + obj.cost.input;
            
            % Add if there are obstacles!
            if has_obstacles
                obj.cost.obstacles = subs(obj.cost.obstacles, cost_syms.car_radius, r);

                for obstacle = obstacles
                    obj.cost.total = obj.cost.total + subs(obj.cost.obstacles, cost_syms.obstacle_position, obstacle);
                end
            end
            
            % First Derivatives
            obj.true_symbolic_derivatives.l_u = jacobian(obj.cost.total, obj.state_vars.input).';
            obj.true_symbolic_derivatives.l_x = jacobian(obj.cost.total, obj.state_vars.position).';

            % Second Derivatives - Diagonals
            obj.true_symbolic_derivatives.l_uu = jacobian(obj.true_symbolic_derivatives.l_u, obj.state_vars.input);
            obj.true_symbolic_derivatives.l_xx = jacobian(obj.true_symbolic_derivatives.l_x, obj.state_vars.position);

            % Second Derivatives - Off Diagonal
            obj.true_symbolic_derivatives.l_ux = jacobian(obj.true_symbolic_derivatives.l_u, obj.state_vars.position);
            obj.true_symbolic_derivatives.l_xu = diff(obj.true_symbolic_derivatives.l_x, obj.state_vars.input);
        
            % Final Derivatives
            % Need to be 100% on these!
            obj.true_symbolic_derivatives.l_final_x = jacobian(obj.cost.state_final, obj.state_vars.position).';
            obj.true_symbolic_derivatives.l_final_xx = jacobian(obj.true_symbolic_derivatives.l_final_x, obj.state_vars.position);

            % Linearized Dynamics
            obj.true_symbolic_derivatives.f_x = jacobian(obj.dynamics, obj.state_vars.position);
            obj.true_symbolic_derivatives.f_u = jacobian(obj.dynamics, obj.state_vars.input);

            obj.true_symbolic_derivatives.f_xx = diff(obj.true_symbolic_derivatives.f_x, obj.state_vars.position(3));

            obj.true_symbolic_derivatives.f_uu = 0;
            obj.true_symbolic_derivatives.f_ux = 0;
            obj.true_symbolic_derivatives.f_xu = 0;

            % Set up simulation
            obj.simulation_vars.time_step = simulation_vars.time_step;
            obj.simulation_vars.time_horizon = simulation_vars.time_horizon;
            obj.simulation_vars.max_iterations = simulation_vars.max_iterations;
        end

        % Substitute Cost Matrices
        function obj = setCostMatrices(obj, Q, R, Q_tf, x_g)
            fieldNames = fieldnames(obj.true_symbolic_derivatives);
            for i = 1:size(fieldNames)
                obj.symbolic_derivatives.(fieldNames{i}) = subs(obj.true_symbolic_derivatives.(fieldNames{i}), obj.cost_vars.position, Q);
                obj.symbolic_derivatives.(fieldNames{i}) = subs(obj.symbolic_derivatives.(fieldNames{i}), obj.cost_vars.input, R);
                obj.symbolic_derivatives.(fieldNames{i}) = subs(obj.symbolic_derivatives.(fieldNames{i}), obj.cost_vars.position_final, Q_tf);
                obj.symbolic_derivatives.(fieldNames{i}) = subs(obj.symbolic_derivatives.(fieldNames{i}), obj.state_vars.position_goal, x_g);
                obj.symbolic_derivatives.(fieldNames{i}) = subs(obj.symbolic_derivatives.(fieldNames{i}), obj.state_vars.dt, obj.simulation_vars.time_step);
                
                % Turn into function
                obj.symbolic_derivatives.(fieldNames{i}) = matlabFunction(obj.symbolic_derivatives.(fieldNames{i}), 'Vars',{obj.state_vars.position(1), obj.state_vars.position(2),obj.state_vars.position(3), obj.state_vars.input});
            end
            obj.cost.total = subs(obj.cost.total, obj.cost_vars.position, Q);
            obj.cost.total = subs(obj.cost.total, obj.cost_vars.input, R);
            obj.cost.total = subs(obj.cost.total, obj.state_vars.position_goal, x_g);
            obj.cost.total = matlabFunction(obj.cost.total, 'Vars',{obj.state_vars.position(1), obj.state_vars.position(2),obj.state_vars.position(3), obj.state_vars.input});

            obj.cost.state_final = subs(obj.cost.state_final, obj.cost_vars.position_final, Q_tf);
            obj.cost.state_final = subs(obj.cost.state_final, obj.state_vars.position_goal, x_g);
            obj.cost.state_final = matlabFunction(obj.cost.state_final, 'Vars',{obj.state_vars.position(1), obj.state_vars.position(2),obj.state_vars.position(3), obj.state_vars.input});
        end

        % Substitute x, u
        function obj = setState(obj, x, u)
            obj.derivatives.l_u = obj.symbolic_derivatives.l_u(x(1), x(2), x(3), u);
            obj.derivatives.l_x = obj.symbolic_derivatives.l_x(x(1), x(2), x(3), u);
            
            obj.derivatives.l_uu = obj.symbolic_derivatives.l_uu(x(1), x(2), x(3), u);
            obj.derivatives.l_xx = obj.symbolic_derivatives.l_xx(x(1), x(2), x(3), u);

            obj.derivatives.l_ux = obj.symbolic_derivatives.l_ux(x(1), x(2), x(3), u);
            obj.derivatives.l_xu = obj.symbolic_derivatives.l_xu(x(1), x(2), x(3), u);

            obj.derivatives.l_final_x = obj.symbolic_derivatives.l_final_x(x(1), x(2), x(3), u);
            obj.derivatives.l_final_xx = obj.symbolic_derivatives.l_final_xx(x(1), x(2), x(3), u);

            obj.derivatives.f_x = obj.symbolic_derivatives.f_x(x(1), x(2), x(3), u);
            obj.derivatives.f_u = obj.symbolic_derivatives.f_u(x(1), x(2), x(3), u);

            obj.derivatives.f_xx = obj.symbolic_derivatives.f_xx(x(1), x(2), x(3), u);
            obj.derivatives.f_uu = obj.symbolic_derivatives.f_uu(x(1), x(2), x(3), u);

            obj.derivatives.f_xu = obj.symbolic_derivatives.f_xu(x(1), x(2), x(3), u);
            obj.derivatives.f_ux = obj.symbolic_derivatives.f_ux(x(1), x(2), x(3), u);
        end

        % Q terms - iLQR
        function obj = setQTermsILQR(obj)
            obj.q_terms.q_x = obj.derivatives.l_x + obj.derivatives.f_x.' * obj.backwards_update_terms.V_x; % TRANSPOSE OK
            obj.q_terms.q_u = obj.derivatives.l_u + obj.derivatives.f_u.' * obj.backwards_update_terms.V_x; % TRANSPOSE OK
            obj.q_terms.q_xx = obj.derivatives.l_xx + obj.derivatives.f_x.' * obj.backwards_update_terms.V_xx * obj.derivatives.f_x;
            obj.q_terms.q_uu = obj.derivatives.l_uu + obj.derivatives.f_u.' * obj.backwards_update_terms.V_xx * obj.derivatives.f_u;
            obj.q_terms.q_xu = obj.derivatives.l_xu + obj.derivatives.f_x.' * obj.backwards_update_terms.V_xx * obj.derivatives.f_u;
            obj.q_terms.q_ux = obj.derivatives.l_ux + obj.derivatives.f_u.' * obj.backwards_update_terms.V_xx * obj.derivatives.f_x;
        end
        
        % Q terms - DDP
        function obj = setQTermsDDP(obj)
            obj = obj.setQTermsILQR();
            obj.q_terms.q_xx = obj.q_terms.q_xx + obj.backwards_update_terms.V_x .* obj.derivatives.f_xx;
            obj.q_terms.q_uu = obj.q_terms.q_uu; 
            obj.q_terms.q_xu = obj.q_terms.q_xu;
            obj.q_terms.q_ux = obj.q_terms.q_ux;
        end

        % Identical for DDP + iLQR
        function obj = setGains(obj)
            obj.gains.K = -inv(obj.q_terms.q_uu) * obj.q_terms.q_ux;
            obj.gains.k = -inv(obj.q_terms.q_uu) * obj.q_terms.q_u;
        end

        function obj = backwardsUpdate(obj)
            obj.backwards_update_terms.V_x = obj.q_terms.q_x + obj.gains.K.' * obj.q_terms.q_uu * obj.gains.k + obj.gains.K.' * obj.q_terms.q_u + obj.q_terms.q_ux.' * obj.gains.k;
            obj.backwards_update_terms.V_xx = obj.q_terms.q_xx + obj.gains.K.' * obj.q_terms.q_uu * obj.gains.K + obj.gains.K.' * obj.q_terms.q_ux + obj.q_terms.q_ux.' * obj.gains.K;
        end

        function [feedback_gains, feedforward_gains, iteration_costs, x_trajectory, u_trajectory] = IADP(obj, ddp)
            % K (1x3) (Will need to transpose!!)
            feedback_gains = zeros([3, obj.simulation_vars.time_horizon]);

            % k (1x1)
            feedforward_gains = zeros([1, obj.simulation_vars.time_horizon]);

            iteration_costs = zeros([1, obj.simulation_vars.max_iterations + 1]);
            iteration_cost_prev = inf;
            
            % x_trajectory (3x1)
            x_trajectory = zeros([3, obj.simulation_vars.time_horizon + 1]);

            % u_trajectory (1x1)
            u_trajectory = zeros([1, obj.simulation_vars.time_horizon]);
            for iteration = 0:obj.simulation_vars.max_iterations
                index = iteration + 1;
                [x_trajectory_new, u_trajectory_new, iteration_cost] = obj.forwardsPass(feedback_gains, feedforward_gains, x_trajectory, u_trajectory);
                

                if true %iteration_cost < iteration_cost_prev
                    x_trajectory = x_trajectory_new;
                    u_trajectory = u_trajectory_new;
                    iteration_cost_prev = iteration_cost;
                else
                    iteration_cost = iteration_cost_prev;
                    iteration_cost_prev = iteration_cost;
                end

                if ddp
                    [feedback_gains, feedforward_gains] = obj.backwardsPassDDP(x_trajectory, u_trajectory);
                else
                    [feedback_gains, feedforward_gains] = obj.backwardsPassiLQR(x_trajectory, u_trajectory);
                end
                
                iteration_costs(index) = iteration_cost;
                disp(["Done with iteration: ", iteration])
            end

            disp("Done with IADP! :)")
        end
           
        function [feedback_gains, feedforward_gains] = backwardsPassiLQR(obj, x_trajectory, u_trajectory)
            % K (1x3) (Will need to transpose!!)
            feedback_gains = zeros([3, obj.simulation_vars.time_horizon]);

            % k (1x1)
            feedforward_gains = zeros([1, obj.simulation_vars.time_horizon]);
            
            % Initialize backwards update terms based on final costs
            obj = obj.setState(x_trajectory(:, end), u_trajectory(:, end));
            obj.backwards_update_terms.V_x = obj.derivatives.l_final_x;
            obj.backwards_update_terms.V_xx = obj.derivatives.l_final_xx;
            obj = obj.setQTermsILQR();
            obj = obj.setGains();

            for step = obj.simulation_vars.time_horizon:-1:1
                obj = obj.backwardsPassStepILQR(x_trajectory(:, step), u_trajectory(:, step));
                feedback_gains(:, step) = obj.gains.K.';
                feedforward_gains(step) = obj.gains.k;
            end
        end

        function [feedback_gains, feedforward_gains] = backwardsPassDDP(obj, x_trajectory, u_trajectory)
            % K (1x3) (Will need to transpose!!)
            feedback_gains = zeros([3, obj.simulation_vars.time_horizon]);

            % k (1x1)
            feedforward_gains = zeros([1, obj.simulation_vars.time_horizon]);
            
            % Initialize backwards update terms based on final costs
            obj = obj.setState(x_trajectory(:, end), u_trajectory(:, end));
            obj.backwards_update_terms.V_x = obj.derivatives.l_final_x;
            obj.backwards_update_terms.V_xx = obj.derivatives.l_final_xx;
            obj = obj.setQTermsDDP();
            obj = obj.setGains();

            for step = obj.simulation_vars.time_horizon:-1:1
                obj = obj.backwardsPassStepDDP(x_trajectory(:, step), u_trajectory(:, step));
                feedback_gains(:, step) = obj.gains.K.';
                feedforward_gains(step) = obj.gains.k;
            end
        end

        function obj = backwardsPassStepILQR(obj, x, u)
            obj = obj.setState(x, u);
            obj = obj.setQTermsILQR();
            obj = obj.setGains();
            obj = obj.backwardsUpdate();
        end

        function obj = backwardsPassStepDDP(obj, x, u)
            obj = obj.setState(x, u);
            obj = obj.setQTermsDDP();
            obj = obj.setGains();
            obj = obj.backwardsUpdate();
        end

        function [x_trajectory, u_trajectory, iteration_cost] = forwardsPass(obj, feedback_gains, feedforward_gains, prev_trajectory_x, prev_trajectory_u)
            % x_trajectory (3x1)
            x_trajectory = zeros([3, obj.simulation_vars.time_horizon + 1]);

            % u_trajectory (1x1)
            u_trajectory = zeros([1, obj.simulation_vars.time_horizon]);

            % iteration cost (1x1)
            iteration_cost = 0;
            
            % Set initial location to initial state
            x_trajectory(:, 1) = obj.initial_state;

            iteration_cost = iteration_cost + obj.cost.total(obj.initial_state(1), obj.initial_state(2), obj.initial_state(3), 0);

            for step = 1:obj.simulation_vars.time_horizon
                next_pos = x_trajectory(:, step);
                delta_x = next_pos - prev_trajectory_x(:, step);
                delta_u = feedback_gains(:, step).' * delta_x + feedforward_gains(step);
                new_u = delta_u + prev_trajectory_u(step);
                obj.setState(next_pos, next_pos(3));
                obj = obj.dynamicsFn(new_u);
                new_x = obj.state;

                x_trajectory(:, step + 1) = new_x;
                iteration_cost = iteration_cost + obj.cost.total(new_x(1), new_x(2), new_x(3), new_u);
            end
            iteration_cost = iteration_cost + obj.cost.state_final(new_x(1), new_x(2), new_x(3), 0);
        end

        
        % Vehicle Dynamics Function
        function obj = dynamicsFn(obj, u)
            obj.state = obj.state + obj.simulation_vars.time_step * [obj.velocity * cos(obj.state(3)) obj.velocity * sin(obj.state(3)) u].';
        end
    end
end