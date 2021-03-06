%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tutorial using SINDy or DMD to identify a SEIR model
% Identified model used in MPC infectious disease control
%
% Main parameters to initialize before running: 
% line 16: select model for system identification -> DMD or SINDy
% line 20: choose to run MPC with or without constraint
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Select DMD or SINDy model for system identification
% select_model = 'DMD';
select_model = 'SINDy';

% Choose to run MPC with or without constraint
constrON = 1; % Constraint on
% constrON = 0; % Constraint off


%% Initialize problem

% Parameters of SEIR model
p.beta0	= 0.5;      % Transmitting rate 
p.gamma	= 0.2;      % Recovery rate 
p.k     = 0.2;      % Incubation rate

% Initial conditions
S0      = 0.996;    % Initial susceptible population
E0      = 0.002;    % Initial exposed population 
I0      = 0.002;    % Initial infected population 
R0      = 0;        % Initial recovered population 
x0      = [S0, E0, I0, R0]; % Initial conditions
n       = 4;        % Number of states (compartments)

% Model specific parameters
if strcmp(select_model,'DMD')
    dt = 1.0; % Time step
elseif strcmp(select_model,'SINDy')
    dt        = 0.01; % Time step
    polyorder = 3;    % Library terms polynomial order
    usesine   = 0;    % Using sine functions in library
end

% Duration of differnt simulations
T_train = 100; % Training phase, number of days
T_test  = 300; % Testing phase, number of days
T_opt   = 300; % MPC phase, number of days


%% Define PRBS forcing function to excite system for model identification
Nic     = 2;         % Number of different PRBS signals: 1 testing and 1 training
taulim  = [14 28];   % Update training and testing intervention each 14 to 28 days
states  = 0.2:0.1:1; % Interventions vary between u(t) = [0.3, 1] with 0.1 steps
Nswitch = 1000;      % Number of switches in prbs function
rng(1,'twister');    % Control random number generator
seed    = randi(Nic,Nic,1); % Vary actuation in each trajectory
forcing = @(x,t,seedval) p.beta0*prbs(taulim, Nswitch, states, t,0,seedval);

SeedTrain = 1;  % PRBS signal 1 for training
SeedTest = 2;   % PRBS signal 2 for testing

% Define training and testing forcing u(t)
tspanTrain  = 0:dt:T_train; % Time span training
for i = 1:length(tspanTrain)
    uTrain(i,:) = forcing(0,tspanTrain(i),seed(SeedTrain));
end
tspanTest   = 0:dt:T_test;  % Time span testing
for i = 1:length(tspanTest)
    uTest(i,:) = forcing(0,tspanTest(i),seed(SeedTest));
end


%% Generate Data
% Integrate forced SEIR system
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n)); % Define ode options
[~,xTrain]=ode45(@(t,x) SEIR(t,x,forcing(x,t,seed(SeedTrain)),p),tspanTrain,x0,options); % Integrate SEIR model using ode45

% Plot SEIR training data
figure;
plot(tspanTrain,xTrain,'LineWidth',1.5)
set(gca,'FontSize',14)
title('Training data')
xlabel('Time')
ylabel('Population ratios (Compartments)')
legend('Suscepticle','Exposed','Infected','Recovered')


%% Compute true derivative
eps = 0.0; % Noise level
for i=1:length(xTrain)
    dx(i,:) = SEIR(0,xTrain(i,:),uTrain(i),p); % Run SEIR model to compute derivatives
end
dx(:,n+1) = 0*dx(:,n);          % Add additional column of zeros for derivatives of u, needed in SINDy algorithm
dx = dx + eps*randn(size(dx));  % Add measurement noise
xTrain = [xTrain uTrain];       % Stack state and input measurements


%% Model Identification
if strcmp(select_model,'DMD')
    [sysmodel_DMDc,U,Up] = DMDc(xTrain(:,1:4)',uTrain',dt); % Identify DMDc model
elseif strcmp(select_model,'SINDy')
    Theta = poolData(xTrain,(n+1),polyorder,usesine); % Generate library
    lambda = 1e-1;     % lambda is our sparsification hyperparameter.
    Xi = sparsifyDynamics(Theta,dx,lambda,(n+1)); % Run SINDy using sequential treshold least squares
    poolDataLIST({'S','E','I','R','u'},Xi,(n+1),polyorder,usesine); % List library terms
end


%% Run validation
% True dynamics: integrate forced true SEIR model using ode45
[~,xTrueTest]=ode45(@(t,x)SEIR(t,x,forcing(x,t,seed(SeedTest)),p),tspanTest,x0,options);

% DMD or SINDy dynamics
if strcmp(select_model,'DMD')
    % DMD dynamics: integrate DMD model using linear simulation tool lsim
    [xModelTrain,~] = lsim(sysmodel_DMDc,uTrain,tspanTrain,x0); % Training
    [xModelTest,~] = lsim(sysmodel_DMDc,uTest,tspanTest,x0); % Testing
elseif strcmp(select_model,'SINDy')
    % SINDy dynamics: integrate SINDy model using ode45
    % the SINDy ode with control is computed using the function sparseGalerkinControl
    p.ahat = Xi(:,1:n);      % SINDy model terms
    p.polyorder = polyorder; 
    p.usesine = usesine;
    [~,xModelTrain]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t,seed(SeedTrain)),p),tspanTrain,x0,options); % Training
    [~,xModelTest]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t,seed(SeedTest)),p),tspanTest,x0,options); % Testing
end

% Integrate unforced true SEIR dynamics for comparison
[~,xUnforced]=ode45(@(t,x)SEIR(t,x,p.beta0,p),tspanTest,x0,options);   

% Plot model comparison
clear ph
figure,box on,
ccolors = get(gca,'colororder');
plot(tspanTest,xTrueTest(:,1),'-','Color',ccolors(1,:),'LineWidth',1); hold on
plot(tspanTest,xTrueTest(:,2),'-','Color',ccolors(2,:),'LineWidth',1);
plot(tspanTest,xTrueTest(:,3),'-','Color',ccolors(3,:),'LineWidth',1);
plot(tspanTest,xTrueTest(:,4),'-','Color',ccolors(4,:),'LineWidth',1);
plot(tspanTest,xModelTest(:,1),'--','Color',ccolors(1,:)-[0 0.2 0.2],'LineWidth',2);
plot(tspanTest,xModelTest(:,2),'--','Color',ccolors(2,:)-[0.1 0.2 0.09],'LineWidth',2);
plot(tspanTest,xModelTest(:,3),'--','Color',ccolors(3,:)-[0 0.2 0.09],'LineWidth',2);
plot(tspanTest,xModelTest(:,4),'--','Color',ccolors(4,:)-[0.1 0.1 0.09],'LineWidth',2);
ylim([0 1])
xlabel('Time, days')
ylabel('Population ratios (Compartments)')
set(gca,'LineWidth',1, 'FontSize',14)
title(['Comparing identified ' select_model ' vs. true dynamics'])
legend('Suscepticle','Exposed','Infected','Recovered',['Suscepticle_{',select_model,'}'],['Exposed_{',select_model,'}'],['Infected_{',select_model,'}'],['Recovered_{',select_model,'}'])


%% Initialize MPC
Nvec        = 14;           % Choose prediction horizon over which the optimization is performed
Imax        = 0.05;         % Constraint: infectious population <= zMax  
Ts          = 1;            % Sampling time
Ton         = 7;            % Time when control starts
Tcontrol    = 7;            % Time when control updates: once a week
Duration    = T_opt;        % Run control for 100 time units
Q           = [0 0 1 0];	% State weights
R           = 0.1;          % Control variation du weights
Ru          = 0.1;          % Control weights
x0n         = x0';          % Initial condition
uopt0       = p.beta0;      % Set initial control input to zero

% Constraints on control optimization
LB          = 0.3*p.beta0*ones(Nvec,1);	% Lower bound of control input
UB          = p.beta0*ones(Nvec,1);     % Upper bound of control input

% Reference state, which shall be achieved: only deviation from xref1(3) (infectious cases) is penalized
xref1 = [0; 0; 0; 0];

% MPC parameters used in objective and constraint function
if strcmp(select_model,'DMD')
    pMPC.sys = sysmodel_DMDc;
elseif strcmp(select_model,'SINDy')
    pMPC.ahat = Xi(:,1:n);
    pMPC.polyorder = polyorder;
    pMPC.usesine = usesine;
end
pMPC.Ts = Ts;
pMPC.beta0 = p.beta0;
pMPC.Imax = Imax;   
    

%% run MPC    
% Options for optimization routine
options = optimoptions('fmincon','Display','none');%,'MaxIterations',100);

% Start simulation
x        = x0n;     % Initial state
uopt     = uopt0.*ones(Nvec,1); % Initial control input
xHistory = x;       % Stores state history
uHistory = uopt(1); % Stores control history
tHistory = 0;       % Stores time history
rHistory = xref1;   % Stores reference (could be trajectory and vary with time)


% For each iteration: take measurements & optimize control input & apply control input
for ct = 1:(Duration/Ts)   
    if ct*Ts>Ton            % Turn control on
        if mod(ct*Ts,Tcontrol) == 0  % Update (Optimize) control input: once a week
            % Set references
            xref = xref1;

            % NMPC with full-state feedback
            COSTFUN = @(u) objectiveFCN(u,x,Nvec,xref,uopt(1),pMPC,diag(Q),R,Ru,select_model); % Cost function
            if constrON
                CONSFUN = @(u) constraintFCN(u,x,Nvec,pMPC,select_model); % Constraint function
                uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,options); % Optimization with constraint
            else
                uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,[],options); % Optimization without constraint
            end
        end
    else    % If control is off
        uopt = uopt0.*ones(Nvec,1);
        xref = [nan; nan; nan; nan];
    end

    % Integrate system: Apply control & Step one timestep forward
    x = rk4u(@SEIR,x,uopt(1),Ts/1,1,[],p);      % Runge-Kutta scheme of order 4 for control system
    xHistory = [xHistory x];                    % Store current state
    uHistory = [uHistory uopt(1)];              % Store current control
    tHistory = [tHistory tHistory(end)+Ts/1];   % Store current time
    rHistory = [rHistory xref];                 % Store current reference
end

% Evaluate objective function
J = evaluateObjectiveFCN(uHistory,xHistory,rHistory,diag(Q),R,Ru,pMPC);


%% Plot results

% State: infectious cases I(t)
figure
plot(tspanTest,xUnforced(:,3),'r','LineWidth',2); hold on
plot(tHistory,xHistory(3,:),'b','LineWidth',2);
plot([tHistory(1),tHistory(end)],[Imax Imax],'k--','LineWidth',2);
xlabel('Time')
ylabel('Infectious population')
set(gca,'FontSize',14)
legend('No control','MPC','Constraint')

% Control input u(t)
figure
plot(tHistory,p.beta0-uHistory,'b','LineWidth',2);
xlabel('Time')
ylabel('Control')
set(gca,'FontSize',14)
legend('Control input u(t)')

% Cost J(t)
figure
plot(tHistory,J,'g','LineWidth',2)
xlabel('Time')
ylabel('Cost')
set(gca,'FontSize',14)
legend('Cost J(t)')



%% Save results
% Training
Results.tTrain = tspanTrain;
Results.xTrueTrain = xTrain;
Results.xModelTrain = xModelTrain;
Results.uTrain = uTrain;

% Testing
Results.tTest = tspanTest;
Results.xTrueTest = xTrueTest;
Results.xModelTest = xModelTest;
Results.uTest = uTest;

% MPC
Results.t = tHistory;
Results.x = xHistory;
Results.u = uHistory;
Results.J = J;

% Unforced
Results.tUnforced = tspanTest;
Results.xUnforced = xUnforced;

% Save
save(['../Results/Results_SEIR_MPC_' select_model '_C' num2str(constrON) '.mat'],'Results')


%% Plot figure 3 from tutorial paper
% First, run 4 cases (DMD and SINDy with and without constraint) 
plotFigures

