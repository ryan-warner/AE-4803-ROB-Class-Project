%% Gradient Calc Function
function[dxH, duH, dx2H, du2H, dudxH] = gradientCalc(cost, dynamics, dxVp1, dx2Vp1, state, control, DPDCalc, funcInput)
% trying to manually impose transposes
% dxc = jacobian(cost, state).';  % moved to DPD
% dxc = DPDCalc{1};
% % dxf = jacobian(dynamics, state);  % moved to DPD
% dxf = DPDCalc{2};
% dxH = dxc + dxf.' * dxVp1;  % OUTPUT, changing some transposes made this work
% % duc = jacobian(cost, control).';    % put in transposes here, but not necessary when control is 1x1, moved to DPD
% duc = DPDCalc{3};
% % duf = jacobian(dynamics, control);      % moved to DPD
% duf = DPDCalc{4};
% duH = duc + duf.' * dxVp1;  % OUTPUT

%% numeric section
% starting V plus 1 values are now being input as evaluated+

% pulling out the numeric functions from the cell array
num_dxc = DPDCalc{1}; 
num_dxf = DPDCalc{2};
num_duc = DPDCalc{3};
num_duf = DPDCalc{4};
num_dx2c  = DPDCalc{5};
num_du2c = DPDCalc{6};
num_dudxc = DPDCalc{7};
num_dx2f = DPDCalc{8};  % tricky
num_du2f = DPDCalc{9};
num_dudxf = DPDCalc{10};    % tricky

z = funcInput(1);
y = funcInput(2);
theta = funcInput(3);
u = funcInput(4);

eval_dxc = num_dxc(z, y, theta, u);
eval_dxf = num_dxf(z, y, theta, u);

dxH = eval_dxc + eval_dxf.' * dxVp1;  % OUTPUT

eval_duc = num_duc(z, y, theta, u);
eval_duf = num_duf(z, y, theta, u);

duH = eval_duc + eval_duf.' * dxVp1;  % OUTPUT

eval_dx2c = num_dx2c(z, y, theta, u);
eval_dx2f = num_dx2f(z, y, theta, u);


dx2H = zeros(length(state), length(state)); % OUTPUT
for i = 1:length(state)       % row number    
    for j = 1:length(state)     % column number
        dxidxjc = eval_dx2c(i,j);
        dxif = eval_dxf(:,i);
        dxjf = eval_dxf(:,j);
        dxidxjf = eval_dx2f(1+3*(i-1):3*i,j);
        dx2H(i,j) = dxidxjc + dxif.'*dx2Vp1 * dxjf + dxVp1.'*dxidxjf;   
    end
end

eval_du2c = num_du2c(z, y, theta, u);
eval_du2f = num_du2f(z, y, theta, u);

du2H = zeros(length(control), length(control)); % OUTPUT
for i = 1:length(control)       % row number    
    for j = 1:length(control)     % column number
        duidujc = eval_du2c(i,j);
        duif = eval_duf(:,i);
        dujf = eval_duf(:,j);
        duidujf = eval_du2f;
        du2H(i,j) = duidujc + duif.'*dx2Vp1 * dujf + dxVp1.'*duidujf;   % removed {} when changing from cell to sym array
    end
end

eval_dudxc = num_dudxc(z, y, theta, u);
eval_dudxf = num_dudxf(z, y, theta, u);

dudxH = zeros(length(control), length(state));  % OUTPUT
for i = 1:length(control)       % row number    
    for j = 1:length(state)     % column number
        % duidxjc = jacobian(jacobian(cost, state(j)), control(i));
        duidxjc = eval_dudxc(i,j);
        % duif = jacobian(dynamics, control(i));
        duif = eval_duf(:,i);
        % dxjf = jacobian(dynamics, state(j));
        dxjf = eval_dxf(:,j);
        % duidxjf = jacobian(jacobian(dynamics, state(j)), control(i));
        duidxjf = eval_dudxf(1+3*(i-1):3*i,j);
        dudxH(i,j) = duidxjc + duif.'*dx2Vp1 * dxjf + dxVp1.'*duidxjf;  % removed {} when changing from cell to sym array
    end
end


%% old work
% 
% % creating placeholder A to initialize symbolic matrices
% syms A B C
% 
% % section for dx2H, should be a 3 x 3 matrix, making as a cell array
% % dx2H = cell(length(state), length(state));
% dx2H = sym('A', [length(state) length(state)]); % OUTPUT
% dx2c = DPDCalc{5};
% dxf = DPDCalc{6};
% dx2f = DPDCalc{10};
% for i = 1:length(state)       % row number    
%     for j = 1:length(state)     % column number
%         % dxidxjc = jacobian(jacobian(cost, state(j)), state(i)); moved to
%         dxidxjc = dx2c(i,j);
%         % dxif = jacobian(dynamics, state(i)); moved to DPD
%         dxif = dxf(:,i);
%         % dxjf = jacobian(dynamics, state(j)); moved to DPD
%         dxjf = dxf(:,j);
%         % dxidxjf = jacobian(jacobian(dynamics, state(j)), state(i)); % had
%         % to use cells to remove
%         dxidxjf = dx2f{i,j};
%         dx2H(i,j) = dxidxjc + dxif.'*dx2Vp1 * dxjf + dxVp1.'*dxidxjf;   % removed {} when changing from cell to sym array
%         % adding a transpose to Vp1, like in the last loop, i think only
%         % way to make dimensions work
%     end
% end
% % seems to work but idrk
% % showdx2H = dx2H
% 
% % building du2H, doing robustly even though not necessary
% % du2H = cell(length(control), length(control));
% du2H = sym('B', [length(control), length(control)]);    % OUTPUT
% du2c = DPDCalc{7};
% duf = DPDCalc{8};
% du2f = DPDCalc{11};
% for i = 1:length(control)       % row number    
%     for j = 1:length(control)     % column number
%         % duidujc = jacobian(jacobian(cost, control(j)), control(i));
%         duidujc = du2c(i,j);
%         % duif = jacobian(dynamics, control(i));
%         duif = duf(:,i);
%         dujf = jacobian(dynamics, control(j));
%         dujf = duf(:,j);
%         % duidujf = jacobian(jacobian(dynamics, control(j)), control(i));
%         duidujf = du2f;
%         du2H(i,j) = duidujc + duif.'*dx2Vp1 * dujf + dxVp1.'*duidujf;   % removed {} when changing from cell to sym array
%     end
% end
% % showdu2H = du2H     % this was mostly copy paste, so may need to change
% 
% % building dudxH as cell array
% % dudxH = cell(length(control), length(state));
% dudxH = sym('C', [length(control), length(state)]); % OUTPUT
% dudxc = DPDCalc{9};
% dudxf = DPDCalc{12};
% for i = 1:length(control)       % row number    
%     for j = 1:length(state)     % column number
%         % duidxjc = jacobian(jacobian(cost, state(j)), control(i));
%         duidxjc = dudxc(i,j);
%         % duif = jacobian(dynamics, control(i));
%         duif = duf(:,i);
%         % dxjf = jacobian(dynamics, state(j));
%         dxjf = dxf(:,j);
%         % duidxjf = jacobian(jacobian(dynamics, state(j)), control(i));
%         duidxjf = dudxf{i,j};
%         dudxH(i,j) = duidxjc + duif.'*dx2Vp1 * dxjf + dxVp1.'*duidxjf;  % removed {} when changing from cell to sym array
%     end
% end
% % showdudxH = dudxH

end