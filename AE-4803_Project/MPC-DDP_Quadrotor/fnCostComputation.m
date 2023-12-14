function [Cost] =  fnCostComputation(x_traj,u_new,Horizon,p_target,dt,Q_f,R)

 Cost = 0;
 
 %for j = 1:(Horizon-1)
 for j = 1:(Horizon-2)
     
    Cost = Cost + 0.5 * u_new(:,j)' * R * u_new(:,j) * dt;
     
 end
 
 TerminalCost = (x_traj(:,Horizon) - p_target)' * Q_f * (x_traj(:,Horizon) - p_target);
 
 Cost = Cost + TerminalCost;

end