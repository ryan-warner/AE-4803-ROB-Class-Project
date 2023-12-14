function [x_traj, u_k, L_k, l_k, Cost ] = fnDDP(num_iter, numStates, numInputs, Horizon, dt, Q, Q_f, R, x0, u_k, du_k, x_traj, p_target, gamma)

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
            
            [dfx,dfu] = fnState_And_Control_Transition_Matrices(x_traj(:,j),u_k(:,j));
           
            A(:,:,j) = eye(numStates,numStates) + dfx*dt;
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

%             Attempts to regularize H            
%             if any(eig(H) < 1e-3); or
%             if cond(H) > 1e3
%                 lambda = abs(min(eig(H))) + 1e-3;
%                 H = H + lambda*eye(length(H));
%             end
            
%             if any(isnan(1\H), 'all')
%                 hdiag = [0.05 0.05 0.05 0.05];
%                 inv_H = diag(h_diag);
%             else 
%                 inv_H = 1\H;
%             end

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
           du         = l_k(:,i) + L_k(:,:,i)*dx;
           dx         = A(:,:,i)*dx + B(:,:,i)*du;
           u_new(:,i) = u_k(:,i) + gamma*du;

           % prevent control from going negative
           u_new(:,i) = max(u_new(:,i), [0;0;0;0]);
        end
    
    u_k = u_new;
    
    
    %%%%%%%%%%%%%%%%% Simulation of the Nonlinear System %%%%%%%%%%%%%%%%%%
    [x_traj]    = fnsimulate(x0,u_k,Horizon,dt,0);
    [Cost(:,k)] = fnCostComputation(x_traj,u_k,p_target,dt,Q_f,R);
    
    %fprintf('DDP Iteration %d,  Current Cost = %e \n',k,Cost(1,k));
    
    end %% end iterating over the DDP algorithm


end