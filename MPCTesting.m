% test file for implementing Recursive MPC
tf = 8;     % seconds
deltaT = .01;
finalStep = tf / deltaT;
finalControl = finalStep - 1;   % final time step for control to be input, tf-1
xinit = 0;      % NEEDS TO BE ASSIGNED
inputInit = 0;  % NEEDS TO BE ASSIGNED

stateCurr = xinit;
numIterations = 3;  % number of iterations of DDP to be ran
prevInputs = inputInit;     % some initial guess for the control, maybe from a full iteration warm DDP start
stateVec = zeros(length(xinit), finalStep+1);   % initializing a state vector to track, not sure on dimensions
stateVec(:,1) = xinit;
inputVec = zeros(length(inputInit), finalControl+1);    % initializing control vector
inputVec(:,1) = inputInit;

for t = 0:finalControl      % vector from 1:799, length 799
    
    timeRem = finalControl - t; % number of control inputs remaining
    nextInputs = DDP(stateCurr, prevInputs, numIterations, timeRem);    
    % DDP will need to take all of the above in, probably with different
    % syntax, prev inputs passed in as starting point for DDP, timeRem
    % becomes tf for DDP, stateCurr is xinit for DDP, numIterations is the
    % number of iterations the DDP will be ran
    nextState = dynamics(stateCurr, nextInputs(:,1));   % apply determined input from DDP
    appInput = nextInputs(:,1);     % saving the actual input that was applied
    
    stateVec(:,t+2) = nextState;
    stateCurr = nextState;  % could probably combine lines, here for clarity
    prevInputs = nextInputs;

end
