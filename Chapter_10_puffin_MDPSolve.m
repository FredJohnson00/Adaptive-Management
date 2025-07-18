% puffin example: when to expose a colony to harvest

clear all;
close all;

S = [1;2;3]; % system (i.e., colony) states (collapsed, vulnerable, robust)
A = [1;2]; % actions (no-harvest, harvest)
X = rectgrid(A,S); % all state-action pairs
Ix = getI(X,2); % column of X with system states

 % no-harvest transition matrix - null model 
A1 = [0.7 0.0 0.0
         0.3 0.5 0.4
         0.0 0.5 0.6];
  
% harvest transition matrix 
A2 = [1.0 0.6 0.0
         0.0 0.4 0.7
         0.0 0.0 0.3];

P = [A1 A2]; % array of action-dependent transitions

% rewards
R1 = [0; 0; 0]; % no-harvest rewards

R2 = [-99; 1; 1]; % harvest rewards

R = [R1 R2]; % array of action-dependent rewards

% FINITE TIME HORIZON, UNDISCOUNTED
clear model
model.name = 'backward iteration';
model.X = X;
model.P = P;
model.R = R;
model.T = 10; % 10 years
model.d = 1.0; % discount factor
model.Ix = Ix;
options = struct('print',2,'keepall',1); % show full summary of solution
results = mdpsolve(model,options);
actions = X(results.Ixopt,1);
policy = reshape(actions,3,10) % time and state-specific policy

% INFINITE TIME HORIZON, UNDISCOUNTED; RELATIVE (AVERAGE) VALUE
model.name = 'average value';
options=struct('algorithm','f','print',2,'tol',5e-6,'relval',1);
results = mdpsolve(model,options);
X(results.Ixopt,2:-1:1)

% INFINITE TIME HORIZON, DISCOUNTED; MODIFIED POLICY ITERATION
clear model
model.name = 'policy iteration';
model.P = P;
model.R = R;
model.d = 0.99; % discount factor
model.Ix = Ix;
options = struct('print',2); % show full summary of solution
results = mdpsolve(model,options);
policy = X(results.Ixopt,2:-1:1)

% examine long-term frequency of states under the optimal policy
Sfreq = longrunP(results.pstar)

% -------------  passive adaptive management ----------------------------------------------------------------

% alternative model transitions under no harvest
A1a = [0.4 0.0 0.0
           0.6 0.5 0.4
           0.0 0.5 0.6];

P1 = [A1 A2];
P2 = [A1a A2];
bigP = cell(1,2);
bigP{1} = P1; bigP{2} = P2;

bigR = cell(1,2);
bigR{1} = R; bigR{2 }= R;

model.name = 'passive adaptive management';
[v,a,B] = amdppassive(1,bigP,bigR,model,options)
[X(a(:,1))  X(a(:,2))]

i = 10;
[v,a,B] = amdppassive(i,bigP,bigR,model,options);
for j=1:(i+1)
    B(j,1)
    X(a(:,j),[2 1])
end

% -------------  active adaptive management ----------------------------------------------------------------
[b,Pb,Rb,Svalsb,Xb,Ixb]=amdp(10,bigP,bigR,S,X,Ix);
clear modelA
modelA.name='Active Adaptive Management';
modelA.R=Rb;
modelA.P=Pb;
modelA.discount = 0.99;
modelA.Ix=Ixb;
% call basic solver
results=mdpsolve(modelA,options);
Xopt=Xb(results.Ixopt,[3 2 1])

