function result = calcDerivatives(state, input, derivatives)
    [f_x, f_u, f_xx, f_uu, f_xu] = derivatives.dynamics(state, input);
    [l_x, l_u, l_xx, l_uu, l_xu] = derivatives.cost(state, input);
    [l_final_x, l_final_xx] = derivatives.cost_final(state);

    result.f_x = f_x;
    result.f_u = f_u;
    result.f_xx = f_xx;
    result.f_uu = f_uu;
    result.f_xu = f_xu;
    %result.f_ux = f_xu';

    result.l_x = l_x;
    result.l_u = l_u;
    result.l_xx = l_xx;
    result.l_uu = l_uu;
    result.l_xu = l_xu;
    result.l_ux = l_xu';

    result.l_final_x = l_final_x;
    result.l_final_xx = l_final_xx;
end