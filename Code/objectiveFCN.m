function J = ObjectiveFCN(u,x,N,xref,u0,p,Q,R,Ru, select_model)
%% Cost function of nonlinear MPC for Lotka-Volterra system
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1
%   x:      current state at time k
%   Ts:     controller sample time
%   N:      prediction horizon
%   xref:   state references, constant from time k+1 to k+N
%   u0:     previous controller output at time k-1
%
% Output:
%   J:      objective function cost
%

%% Integrate system
if strcmp(select_model,'DMD')
    [xk,~] = lsim(p.sys,[u' 0],[0:N].*p.dt,x);
    xk = xk(2:end,:); xk = xk';
elseif strcmp(select_model,'SINDy')
    Ns = size(x,1);
    xk = zeros(Ns,N+1); xk(:,1) = x;
    for ct=1:N
        % Obtain plant state at next prediction step.
        xk(:,ct+1) = rk4u(@sparseGalerkinControl_Discrete,xk(:,ct),u(ct),p.dt,1,[],p);
    end
    xk = xk(:,2:N+1);
end

%% Cost Calculation
% Set initial plant states, controller output and cost
uk = u(1);
J = 0;
% Loop through each prediction step
for ct=1:N
    % Obtain plant state at next prediction step
    xk1 = xk(:,ct);
    
    % Accumulate state tracking cost from x(k+1) to x(k+N)
    J = J + (xk1-xref)'*Q*(xk1-xref);
    % Accumulate MV rate of change cost from u(k) to u(k+N-1)
    if ct==1
        J = J + (uk-u0)'*R*(uk-u0) + (p.beta0-uk)'*Ru*(p.beta0-uk);
    else
        J = J + (uk-u(ct-1))'*R*(uk-u(ct-1)) + (p.beta0-uk)'*Ru*(p.beta0-uk);
    end
    % Update uk for the next prediction step
    if ct<N
        uk = u(ct+1);
    end
end

