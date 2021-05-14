function [c, ceq] = constraintFCN(u,x,N,p,select_model)
%% Constraint function of nonlinear MPC for Lotka-Volterra system
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1 
%   x:      current state at time k
%   N:      prediction horizon
%   p:      model parameters
%
% Output:
%   c:      inequality constraints applied across prediction horizon
%   ceq:    equality constraints (empty)
%

%% Nonlinear MPC design parameters
% Infectious population <= zMax
zMax = p.Imax;

%% Integrate system
if strcmp(select_model,'DMD')
    [xk,~] = lsim(p.sys,[u',0],[0:N].*p.Ts,x);
    xk = xk(2:end,:); xk = xk'; 
elseif strcmp(select_model,'SINDy')
    Ns = size(x,1);
    xk = zeros(Ns,N+1); xk(:,1) = x;
    for ct=1:N
        % Obtain plant state at next prediction step.
        xk(:,ct+1) = rk4u(@sparseGalerkinControl,xk(:,ct),u(ct),p.Ts,1,[],p);
    end
    xk = xk(:,2:N+1);
end


%% Inequality constraints calculation
c = zeros(N,1);
% Apply N population size constraints across prediction horizon, from time k+1 to k+N
uk = u(1);
for ct=1:N
    % -z + zMin < 0 % lower bound
%     c(ct) = -xk(2,ct)+zMin; %c(2*ct-1)
    % z - zMax < 0 % upper bound
    c(ct) = xk(3,ct)-zMax; % max infected people
    % update  input for next step
    if ct<N
        uk = u(ct+1);
    end
end
%% No equality constraints
ceq = [];

