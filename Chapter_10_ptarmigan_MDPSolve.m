
% Ptarmigan harvest management in East Iceland
clear all;

% Discretization of state and control variables
Nmin = 0; Nmax = 35; Ninc = 1; % population size in spring (thousands)
dmin = 0; dmax = 45; dinc = 5; % hunting season length (days)

% Shocks
%... s (season natural survival - adults) beta method-of-moments for mean of mean s posteriors
[e1,w1] = qnwbeta(5, 21.86, 3.97); % total variation
Es = e1'*w1; % expectation

%... w (winter naturual survival - young) beta method-of-moments for mean of mean w posteriors
[e2,w2] = qnwbeta(5, 4.33, 6.47); % total variation
Ew = e2'*w2;

% relation between (log) fall age ratio and spring pop
e34 = [1.7019 -0.5372]; w34 = 1; % deterministic coefficients
Ebeta0 = e34(1); Ebeta1 = e34(2);
%... random-year effect
[e5,w5] = qnwnorm(5, 0, 0.1137^2);
Esigmar = w5'*e5;

e = rectgrid(e1,e2,e34,e5);
grid = rectgrid(w1,w2,w34,w5);
w = (prod(grid,2));  
sum(w); 

% Define state and action vectors
N = (Nmin:Ninc:Nmax)'; % spring pop
D = (dmin:dinc:dmax)'; % season days
X = rectgrid(N,D);  % state-action pairs

% define alternative models
a1 = 1; b1 = 58.24; % null model parameters for hfcn
a2 = 0.27; b2 = 1.43; % alternative model parameters for hfcn
weights = [0.27 0.73]; % current probabilities of the two models

% state-action functions
rfcn = @(N,beta0,beta1,sigmar) exp(beta0 + beta1.*N./10 + sigmar); % fall age ratio
hfcn = @(D,a,b) a.*D./(b+D); % harvest rate function

Ptran1 = @(N,D,s,w,beta0,beta1,sigmar) max(0, N.*s.*(1-hfcn(D,a1,b1)) .* (rfcn(N,beta0,beta1,sigmar).*w + s)); % spring pop
Ptran2 = @(N,D,s,w,beta0,beta1,sigmar) max(0, N.*s.*(1-hfcn(D,a2,b2)) .* (rfcn(N,beta0,beta1,sigmar).*w + s)); % spring pop
Ptran = {Ptran1,Ptran2};

Hfcn1 = @(N,D,s,beta0,beta1,sigmar) N.*s.*(1 + rfcn(N,beta0,beta1,sigmar)).*hfcn(D,a1,b1);  % harvest
Hfcn2 = @(N,D,s,beta0,beta1,sigmar) N.*s.*(1 + rfcn(N,beta0,beta1,sigmar)).*hfcn(D,a2,b2);  % harvest

% compute transition matrices
coptions = struct('cleanup',0); 
Stran=cell(1,2);
for i=1:2
    G = @(X,e) Ptran{i}(X(:,1),X(:,2),e(:,1),e(:,2),e(:,3),e(:,4),e(:,5)) ;
    Stran{i} = g2P(G,N,X,e,w,coptions);
end
bigP = Stran;

% compute rewards
 penalty = @(N,goal) min(1,N./goal);
 goal = [0.001 12.7 14.6 16.0 17.2]; % none, min, 25%, 50%, mean
 goalndx = 2; 
 
Rfcn1 = log(max(0.001, Hfcn1(X(:,1),X(:,2),Es,Ebeta0,Ebeta1,Esigmar))) ; 
R1 = Rfcn1 .* penalty(X(:,1),goal(goalndx));
Rfcn2 = log(max(0.001, Hfcn2(X(:,1),X(:,2),Es,Ebeta0,Ebeta1,Esigmar))) ; 
R2 = Rfcn2 .* penalty(X(:,1),goal(goalndx));
bigR =cell(1,2);
bigR{1}=R1; bigR{2}=R2;

% compute policy
I = getI(X,1);
modeloptions=struct('d',0.99,'Ix',I);
options=struct('print',2);
intervals = 5;
[v,a,B]=amdppassive(intervals,bigP,bigR,modeloptions,options);
for j=1:(intervals+1)
    B(j,1)
    X(a(:,j),:)
end

