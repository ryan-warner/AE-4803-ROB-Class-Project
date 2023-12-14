function  [l0,l_x,l_xx,l_u,l_uu,l_ux] = fnCost(x,p_target,u, numStates,numInputs, Q, R)

l0   = 0.5*(u - 1.2263)' * R * (u - 1.2263) + 0.5*(x-p_target)' * Q * (x-p_target);
l_x  = Q*(x-p_target);
l_xx = Q;
l_u  = R*(u - 1.2263);
l_uu = R;
l_ux = zeros(numInputs,numStates);