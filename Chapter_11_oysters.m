% OysterAM.m
clear all;
close all;
font="Arial";

% SELECT TYPE OF RUN ------------------------------------------------------
% density dependence
d = [1 2]; % form of DD [1 2] gives nice contrast
% habitat mgmt (CHOOSE ONLY ONE)
behavior = 0;
enhance = 0;
restore = 0; 
% -------------------------------------------------------------------------

% define state and action vectors
Nmin = 0.0; Nmax = 1.01; Ninc = 0.01;
N = (Nmin:Ninc:Nmax)'; % oyster abundance 
Hmin = 0.0; Hmax = 1.01; Hinc = 0.01; % habitat abundance
H = (Hmin:Hinc:Hmax)';
f = (0 : 5 : 250)'; % dynamic # of fishers, or...

S = rectgrid(N,H);
X = rectgrid(S,f);

% model parameters
%...oysters
r = 0.4; % oyster growth rate
%... consumption 
g = 0.001; % max season catch/fisher
h = 0.1; % oyster abundance at g/2 (half-saturation)
s = 0.05; % oyster abundance when fishing abandoned 
%... habitat 
T = 1; % habitat carrying capacity
q_all = [0.3 0.51]; % lo-hi relief production rate
delta_all = [0.16 0.05]; % lo-hi relief loss rate - natural
rho_all = [0.01 0.09]; % scalar to convert fishing rate to habitat loss rate

% habitat scenario
if behavior
    q = q_all(1);    
    delta = delta_all(1);
    rho = rho_all(1);
else
    q = q_all(1);
    delta = delta_all(1);
    rho = rho_all(2);
end
if enhance
    q = q_all(1);    
    delta = delta_all(2);
    rho = rho_all(2);
end
if restore
    q = q_all(2);    
    delta = delta_all(2);
    rho = rho_all(2);
end
    
% define shocks (deterministic)
x=0; y=1;
e1 = x; w1 = y; % oyster shock
e2 = x; w2 = y; % shell shock
e3 = x; w3 = y; % fisher shock
e = rectgrid(e1,e2,e3); w = y;

% define shocks (stochastic)
% cv = 0.05;  
% n = 5; % number of shocks 
% [x,y] = qnwnorm(n,0-cv^2/2,cv^2);
% e1 = x; w1 = y; % oyster shock
% e2 = x; w2 = y; % shell shock
% e3 = x; w3 = y; % catch shock
% e = rectgrid(e1,e2,e3); w = prod(rectgrid(w1,w2,w3),2);
 
% state transitions
Ntran1 = @(N,H,r,d,f,g,s,h,eN,eC)   max(0,     N.* exp(r.*(1-(N./H).^d(1))+eN) - min(N,max(0,exp(eC).*f.*g.*(N-s)./(N-s+h))) );
Ntran2 = @(N,H,r,d,f,g,s,h,eN,eC)   max(0,     N.* exp(r.*(1-(N./H).^d(2))+eN) - min(N,max(0,exp(eC).*f.*g.*(N-s)./(N-s+h))) );
Ntran = {Ntran1,Ntran2};

Htran1 = @(N,H,q,d,delta,rho,f,g,s,h,T,eH,eC) max(0, H + N.*(exp(q.*(1-(H./T).^d(1))+eH)-1) -...
                                                     H.*(exp(delta+rho.*log(min(N,max(0,exp(eC).*f.*g.*(N-s)./(N-s+h)))./N + 1))-1) );
Htran2 = @(N,H,q,d,delta,rho,f,g,s,h,T,eH,eC) max(0, H + N.*(exp(q.*(1-(H./T).^d(2))+eH)-1) -...
                                                     H.*(exp(delta+rho.*log(min(N,max(0,exp(eC).*f.*g.*(N-s)./(N-s+h)))./N + 1))-1) );
Htran = {Htran1,Htran2};

options = struct('cleanup',0); 
Stran=cell(1,2);
for i=1:2
    G = @(X,e) [Ntran{i}(X(:,1),X(:,2),r,d,X(:,3),g,s,h,e(:,1),e(:,3))...
                Htran{i}(X(:,1),X(:,2),q,d,delta,rho,X(:,3),g,s,h,T,e(:,2),e(:,3))];
  Stran{i} = g2P(G,{N,H},X,e,w,options);
end

% Objective function
Rfcn = @(N,f,g,s,h) min(N, max(0,f.*g.*(N-s)./(N-s+h)));
R = Rfcn(X(:,1),X(:,3),g,s,h); 
R = log(R+0.000001); 
U = cell(1,2);
U{1} = R; U{2} = R;

% Computation paramaters
bigP = Stran;
bigR = U;
Ix = getI(X,1:2);

% INFINITE TIME HORIZON, DELTA = 1; RELATIVE (AVERAGE) VALUE
% model=struct('discount',1.0,'Ix',I);
% options=struct('algorithm','f','print',2,'tol',5e-6,'relval',1);

% INFINITE TIME HORIZON, DISCOUNTED; 
%   model=struct('discount',0.9999,'Ix',Ix);
%   options=struct('tol',0.05,'print',2);
 
% INFINITE TIME HORIZON VANISHING DISCOUNT
 model=struct('d',1.0,'Ix',Ix);
 options=struct('algorithm','i','vanish',0.99999,'print',2);

bintervals = 10; % number of belief intervals


%%
% Passive AM
[v,a,B]=amdppassive(bintervals,bigP,bigR,model,options);
%... assemble policy for ease of interpretation
b = B(:,1);
policy=[];
for i=1:size(b,1)
    A=X(a(:,i),3);
    policy=[policy; A];
end
prob = repmat(b',size(S,1),1); prob = prob(:);
policy=[repmat(S,size(b,1),1) prob policy];

dlmwrite('PD.opt',policy);

%%
% Active AM
[b,Pb,Rb,Svalsb,Xb,Ixb]=amdp(bintervals,bigP,bigR,S,X,Ix);
modelA.name='Active Adaptive Management Model';
modelA.R=Rb;
modelA.P=Pb;
modelA.discount=1;
modelA.Ix=Ixb;
%... call basic solver
results=mdpsolve(modelA,options);
v=results.v; Xopt=Xb(results.Ixopt,:);
optval = [Xopt v];
dlmwrite('AD.opt',Xopt);

