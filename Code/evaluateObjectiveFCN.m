function J = evaluateObjectiveFCN(u,x,xref,Q,R,Ru,p)
%% Cost function of nonlinear MPC for SEIR model
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1 
%   x:      current state at time k
%   xref:   state references, constant from time k+1 to k+N
%   Q,R,Ru: MPC cost function weights 
%   p:      model parameters
%
% Output:
%   J:      objective function cost
%

%% Nonlinear MPC design parameters

%% Cost Calculation
% Set initial cost and input
N = size(x,2);
J = zeros(N,1);
u0 = 0;

% Loop through each time step
for ct=1:N
    
    % Accumulate state tracking cost from x(k+1) to x(k+N)
    J(ct) = (x(:,ct)-xref(:,ct))'*Q*(x(:,ct)-xref(:,ct));
    
    % Accumulate MV rate of change cost from u(k) to u(k+N-1)
    if ct==1
        J(ct) = J(ct) + (u(:,ct)-u0)'*R*(u(:,ct)-u0) + (p.beta0-u(:,ct))'*Ru*(p.beta0-u(:,ct));
    else
        J(ct) = J(ct) + (u(:,ct)-u(:,ct-1))'*R*(u(:,ct)-u(:,ct-1)) + (p.beta0-u(:,ct))'*Ru*(p.beta0-u(:,ct));
    end
end

