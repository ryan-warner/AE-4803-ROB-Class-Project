function ctl = controller(state, params)
    ctl = params.K * state + params.d;
end