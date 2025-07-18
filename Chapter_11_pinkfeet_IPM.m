
clear all;
close all;
format long;

popwt = 1.00; % weight on population objective; weight on harvest obj is 1-popwt
goal = 60;

% Discretization of state and control variables
Nmin = 0; Nmax = 100; Ninc = 1;  % population size in May
Dmin = 0; Dmax = 30;  Dinc = 2;   % thaw days  
Hmin = 0; Hmax = 50;  Hinc = 1;   % harvest quotas (optimal harvests implicitly include crippling loss)

% Shocks
%... thaw days 2000 ->
e1 = (Dmin:Dinc:Dmax)';   
w1 = [0.0308 0.0758 0.1035 0.1174 0.1208 0.1161 0.1055 0.0914 0.0753 0.0589 0.0432 0.0293 0.0179 0.0094 0.0038 0.0009]';
%sum(w1)

%... theta (natural survival)
[e2,w2] = qnwbeta(5, 51.272375,  6.865491);
Etheta = e2'*w2;

% relation between preseason age ratio and thaw days
[e34,w34] = qnwnorm([5 5],[-1.8721632  0.4990896],[0.005059301 -0.003919600; -0.003919600  0.003566221]);
Ebeta0 = w34'*e34(:,1); Ebeta1 = w34'*e34(:,2);

% diff vulnerability
[e5,w5] = qnwgamma(5, 1.69619095, 0.02041295, 1);
Edv = e5'*w5;

e = rectgrid(e1,e2,e34,e5);
grid = rectgrid(w1,w2,w34,w5);
w = (prod(grid,2));  
sum(w); 

% Define state and action vectors
N = (Nmin:Ninc:Nmax)'; 
D = (Dmin:Dinc:Dmax)';
S = rectgrid(N,D);
H = (Hmin:Hinc:Hmax)';
X = rectgrid(S,H);  % State,Action combination

% state transition functions
rfcn = @(D,beta0,beta1) (1./(1+exp(-(beta0+beta1.*D./10)))) ./ (1 - (1./(1+exp(-(beta0+beta1.*D./10)))));
PFtran = @(N,D,H,theta,beta0,beta1,dv) max(0,( ...
  N.*theta.*( ...
  ( 1-H./(N.*theta.^(4/12).*(1+dv.*rfcn(D,beta0,beta1)))) + ...
  rfcn(D,beta0,beta1).*(1-dv.* ...
      H./(N.*theta.^(4/12).*(1+dv.*rfcn(D,beta0,beta1)))) ...
   ) ...
));

Dtran = @(N,D,H) ones(size(N,1),1); % for every state/action return "1"; next year's days are given by 1 x e1, with probabilities w1

% compute transition matrices
pf = @(X,e) [max(0, PFtran(X(:,1),X(:,2),X(:,3),e(:,2),e(:,3),e(:,4),e(:,5))) Dtran(X(:,1),X(:,2),X(:,3)).*e(:,1)];
cleanup = 0; % 0 = none; 1 = neg values set to zero, renormalized; 2 = snap to grid
Stran = g2P(pf,{N,D},X,e,w,cleanup);

% objective function (constrain harvest <= N)
 ndx = (X(:,3)>X(:,1) & X(:,3)>0);
 Hutil = X(:,3)./max(H);
 Hutil(ndx) = -inf;
 
  EN = PFtran(X(:,1),X(:,2),X(:,3),Etheta,Ebeta0,Ebeta1,Edv); % expectation of next pop
  poputil = 1./(1+exp(-(10-(abs(EN-goal))))); % utility of next pop
%  poputil = 1./(1+exp(-(10-(abs(X(:,1)-goal))))); % utility of current pop (makes little difference)
  R = popwt*poputil + (1-popwt)*Hutil; 

% compute policy
I = getI(X,1:2);

% INFINITE TIME HORIZON, DELTA = 1; RELATIVE (AVERAGE) VALUE
% modelfull=struct('discount',1.0,'Ix',I);
% options=struct('algorithm','f','print',2,'tol',5e-6,'relval',1);

% INFINITE TIME HORIZON, DISCOUNTED; MODIFIED POLICY ITERATION
 modelfull=struct('d',0.99999,'Ix',I);
 options=struct('print',2);
 
% INFINITE TIME HORIZON VANISHING DISCOUNT
%   modelfull=struct('d',1.0,'Ix',I);
%   options=struct('algorithm','i','vanish',0.99999,'print',2);

modelfull.R = R;
modelfull.P = Stran;
results = mdpsolve(modelfull,options);

pstar = modelfull.P(:,results.Ixopt);
Sfreq=longrunP(pstar);
neg = Sfreq<0; i = find(neg); size(i);
Sfreq(i)=0; Sfreq = Sfreq/sum(Sfreq);
M=marginals(Sfreq,S);
Epop = M{1}'*N
SDpop = sqrt(var(N,M{1}))
Eh = Sfreq'*X(results.Ixopt,3)
SDh = sqrt(var(X(results.Ixopt,3),Sfreq))
Eval = Sfreq'*R(results.Ixopt)
SDval = sqrt(var(R(results.Ixopt),Sfreq))

policy=[X(results.Ixopt,:) results.v];

font="Arial";
x = N;
y = D;
z = reshape(policy(:,3),size(y,1),size(x,1));
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
v=0:2.5:50;
contourf(y,x,z',v)
xlabel('D')
ylabel('N')
colormap(jet(21))
colorbar

% interpolate a harvest quota not located on a grid point
n95 = [51.1 74.3];
n80 = [55.4 70.1];

x = [62.8 9]; % state to be interpolated
g50 = [62 63];
s={g50;[8 10]}; % grid points; this must be a cell array
B = rectbas(x,s,1);
ndx = match([g50(1) 8; g50(1) 10; g50(2) 8; g50(2) 10],policy(:,1:2))
H = [policy(ndx(1),3) policy(ndx(2),3) policy(ndx(3),3) policy(ndx(4),3)]
Hopt = H * B  


gl95 = [51 52];
gu95 = [74 75];
gl80 = [55 56];
gu80 = [70 71];

s={gu95;[8 10]}; % grid points; this must be a cell array
y = [n95(2) 9];
B = rectbas(y,s,1);

lo = gu95(1);
up = gu95(2);
ndx = match([lo 8; lo 10; up 8; up 10],policy(:,1:2))
H = [policy(ndx(1),3) policy(ndx(2),3) policy(ndx(3),3) policy(ndx(4),3)]
Hopt = H * B  
  
