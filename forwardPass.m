% function [states, inputs, trajectoryCostsTotal] = forwardPass(initial, dynamics, costs, gains, options)
%     states = [initial zeros(options.n, options.horizon)];
%     inputs = zeros(options.m, options.horizon);
%     trajectoryCosts = zeros(1, options.horizon + 1);
% 
%     for i = 1:options.horizon
%         iterGains = gains(i);
%         inputs(:, i) =  max(iterGains.K * states(:, i) + iterGains.optimal_control, [0;0;0;0]);
%         trajectoryCosts(i) = costs.cost(states(:, i), inputs(:, i));
%         states(:, i + 1) = dynamics(states(:, i), inputs(:, i)) * options.timestep + states(:, i);
%     end
% 
%     trajectoryCosts(end) = costs.finalCost(states(:, end));
%     trajectoryCostsTotal = sum(trajectoryCosts);
%     if isnan(trajectoryCostsTotal)
%         disp('Final cost is NaN');
%     end
% end

% Attempt at adding line search
function [states, inputs, trajectoryCostsTotal] = forwardPass(initial, dynamics, costs, gains, options)
    stepSizes = [1, 0.5, 0.25, 0.1, 0.05, 0.01, 0.001, 0.0001];
    bestCost = inf;
    states = [initial zeros(12, 800)];
    inputs = zeros(4, 800);
    trajectoryCosts = zeros(1, 800 + 1);

    for stepSize = stepSizes
        tempStates = [initial zeros(12, 800)];
        tempInputs = zeros(4, 800);
        tempCosts = zeros(1, 800 + 1);

        for i = 1:800
            iterGains = gains(i);
            tempInputs(:, i) = iterGains.K * tempStates(:, i) + stepSize * iterGains.k + iterGains.optimal_control;
            tempCosts(i) = costs.cost(tempStates(:, i), tempInputs(:, i));
            tempStates(:, i + 1) = dynamics(tempStates(:, i), tempInputs(:, i));
        end

        tempCosts(end) = costs.finalCost(tempStates(:, end));
        totalCost = sum(tempCosts);

        if totalCost < bestCost
            bestCost = totalCost;
            states = tempStates;
            inputs = tempInputs;
            trajectoryCosts = tempCosts;
        end

        trajectoryCostsTotal = sum(trajectoryCosts);
    end
end