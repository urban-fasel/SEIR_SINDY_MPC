function dx = SEIR(t,x,u,p)
%% SEIR nonlinear dynamical system
% Inputs:
%   t:      time
%   x:      current state
%   u:      control input: u = beta
%   p:      SEIRmodel parameters
%
% Output:
%   dx:     time derivative of states
%

S = x(1); E = x(2); I = x(3); R = x(4);

dx = [  -u*S*I;
        u*S*I - p.k*E;
        p.k*E - p.gamma*I;
        p.gamma*I];
    
end
