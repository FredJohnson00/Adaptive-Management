
% Ptarmigan harvest management in East Iceland
% for background, see Johnson, F. A., and Ó. K. Nielsen. 2024. Regional demography of Icelandic
% rock ptarmigan and its implications for harvest management. Ecological Solutions and Evidence 5:e12390.

clear all;

% Discretization of state and control variables
Nmin = 0; Nmax = 35; Ninc = 1; % population size in spring (thousands)
dmin = 0; dmax = 45; dinc = 5; % hunting season length (days)

% Shocks
%... s (season natural survival - adults) beta method-of-moments for mean of mean s posteriors
[e1,w1] = qnwbeta(5, 21.86, 3.97); % quadrature nodes and weights from a beta distribution
Es = e1'*w1; % expectation

%... w (winter naturual survival - young) beta method-of-moments for mean of mean w posteriors
[e2,w2] = qnwbeta(5, 4.33, 6.47); % quadrature nodes and weights from a beta distribution
Ew = e2'*w2; % expectation

% relation between (log) fall age ratio and spring pop
e34 = [1.7019 -0.5372]; w34 = 1; % deterministic coefficients (i.e., weight = 1)
Ebeta0 = e34(1); Ebeta1 = e34(2); % expectations

%... random-year effect
[e5,w5] = qnwnorm(5, 0, 0.1137^2); % quadrature nodes and weights from a normal distribution
Esigmar = w5'*e5; % expectation

e = rectgrid(e1,e2,e34,e5); % grid of all possible shocks
grid = rectgrid(w1,w2,w34,w5); % grid of weights
w = (prod(grid,2));  % vector of weights for each combination of shocks
sum(w); % ensure weights sum to 1

% Define state and action vectors
N = (Nmin:Ninc:Nmax)'; % spring pop
D = (dmin:dinc:dmax)'; % season days
X = rectgrid(N,D);  % state-action pairs

% define alternative harvest-rate models
a1 = 1; b1 = 58.24; % null model parameters for harvest-rate function, hfcn
a2 = 0.27; b2 = 1.43; % alternative model parameters for hfcn

% state-action functions
rfcn = @(N,beta0,beta1,sigmar) exp(beta0 + beta1.*N./10 + sigmar); % fall age ratio function
hfcn = @(D,a,b) a.*D./(b+D); % harvest rate function

Ptran1 = @(N,D,s,w,beta0,beta1,sigmar) max(0, N.*s.*(1-hfcn(D,a1,b1)) .* (rfcn(N,beta0,beta1,sigmar).*w + s)); % spring pop null model
Ptran2 = @(N,D,s,w,beta0,beta1,sigmar) max(0, N.*s.*(1-hfcn(D,a2,b2)) .* (rfcn(N,beta0,beta1,sigmar).*w + s)); % spring pop alternative model
Ptran = {Ptran1,Ptran2};

Hfcn1 = @(N,D,s,beta0,beta1,sigmar) N.*s.*(1 + rfcn(N,beta0,beta1,sigmar)).*hfcn(D,a1,b1);  % harvest null model
Hfcn2 = @(N,D,s,beta0,beta1,sigmar) N.*s.*(1 + rfcn(N,beta0,beta1,sigmar)).*hfcn(D,a2,b2);  % harvest alternative model

% compute transition matrices
% extrapolation beyond state bounds can produce negative transition probabilities
% type "help g2P" to see options for handling this (note: in most cases neg. transition probs. will not be a problem)
coptions = struct('cleanup',0); % take no action for probs outside the interval [0,1]
Stran=cell(1,2); % set up transitions matrices for 2 models
for i=1:2 % transitions for each of 2 models
    G = @(X,e) Ptran{i}(X(:,1),X(:,2),e(:,1),e(:,2),e(:,3),e(:,4),e(:,5)) ;
    Stran{i} = g2P(G,N,X,e,w,coptions);
end
bigP = Stran;

% compute rewards
 penalty = @(N,goal) min(1,N./goal); % harvest penalty for a population < goal
 goal = [0.001 12.7 14.6 16.0 17.2]; % Population goal: none, min, 25%, 50%, mean of population estimates
 goalndx = 2; % select minimum population as goal
 
Rfcn1 = log(max(0.001, Hfcn1(X(:,1),X(:,2),Es,Ebeta0,Ebeta1,Esigmar))) ; 
R1 = Rfcn1 .* penalty(X(:,1),goal(goalndx)); % returns null model
Rfcn2 = log(max(0.001, Hfcn2(X(:,1),X(:,2),Es,Ebeta0,Ebeta1,Esigmar))) ; 
R2 = Rfcn2 .* penalty(X(:,1),goal(goalndx)); % returns alternative model
bigR =cell(1,2); % set up return vectors for 2 models
bigR{1}=R1; bigR{2}=R2;

% compute policy
I = getI(X,1); % get columns of X with the state variable population size
modeloptions=struct('d',0.99,'Ix',I); % discount factor = 0.99
options=struct('print',2); % full print output
intervals = 5; % intervals of information state probabilities
[v,a,B]=amdppassive(intervals,bigP,bigR,modeloptions,options);
for j=1:(intervals+1)
    B(j,1)
    X(a(:,j),:)
end

