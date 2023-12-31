% Recursive MPC
tf = 8;     % seconds
deltaT = .01;
finalStep = tf / deltaT;
finalControl = finalStep - 1;   % final time step for control to be input, tf-1
xinit = 0;      % NEEDS TO BE ASSIGNED
inputInit = 0;  % NEEDS TO BE ASSIGNED

currState = xinit;
numIterations = 3;  % number of iterations of DDP to be ran
[warmDDP = IADP_Final(initial, dynamics, costs, derivatives, options);    % a warm start for the DDP code, giving full vector of control inputs
prevInputs = warmDDP;     % some initial guess for the control, maybe from a full iteration warm DDP start
stateVec = zeros(length(xinit), finalStep+1);   % initializing a state vector to track, not sure on dimensions
stateVec(:,1) = xinit;
inputVec = zeros(length(warmDDP(:,1)), finalControl+1);    % initializing control vector
% appInput =  prevInputs(:,1);
% inputVec(:,1) = appInput;    % pulling out initial value from warm DDP and saving as applied




% we have input and state at time 0
for t = 1:finalState        % vector 1:800, length 800, making it easier to index
    timeRem = finalControl - t; % number of control inputs remaining
 
    nextInputs = DDP(currState, prevInputs, numIterations, timeRem); % will run the warm start an additional # of times
    appInput = nextInputs(:,1);
    nextState = dynamics(prevState, appInput);
    inputVec(:,t) = appInput;
    stateVec(:,t+1) = nextState;
    prevState = nextState;

    
end

    
    
    
