% PFG20172June2017.m
% Pink-footed goose AHM
% with clean code & annotation

close all
clear variables
clc
disp('Pink-footed goose AHM 1-yr cycle')

% discretization of state and control variables
Ymin = 0; Ymax = 16;  Yinc = 2;   % young
Amin = 0; Amax = 80;  Ainc = 4;   % adults
% Ymin = 1; Ymax = 17;  Yinc = 2;   % young
% Amin = 2; Amax = 80;  Ainc = 3;   % adults
Dmin = 0; Dmax = 28;  Dinc = 4;   % temperature days  ! MAKE SURE THIS CORRESPONDS WITH e3 !
Hmin = 0; Hmax = 50;  Hinc = 2;   % harvest quotas
H = (Hmin:Hinc:Hmax)';                         
Y = (Ymin:Yinc:Ymax)';     
A = (Amin:Ainc:Amax)';
D = (Dmin:Dinc:Dmax)';

% set up transition matrix
S = rectgrid(Y,A,D);
X = rectgrid(S,H);
ns = size(S,1);     
na = length(H);

% define constants that appear in functions
% dv = 2.0;  % differential vulnerability
% z = 1.1;   % assumed value for (1-h)/(1-dh)
% zd = z*dv = 2.2; 

% define random shocks and probability weights 
w = ones(5,1)/5;
e1 = [0.4718; 0.7839; 1.0000; 1.2161; 1.5282]; w1 = w; % young transition shocks (based on MSE of comparing pred & obs)
e2 = [0.8538; 0.9402; 1.0000; 1.0598; 1.1462]; w2 = w; % adult transition shocks (based on MSE of comparing pred & obs)
e3 = [0; 4; 8; 12; 16; 20; 24; 28];                    % temperature-days shocks (beta-binomial distribution)
w3 = [0.0892; 0.3563; 0.3112; 0.1663; 0.0607; 0.0144; 0.0018; 0.0001];
% e3 = [0:3:30]'; % in case you want days in increments of 3
% w3 = [0.066573 0.241292 0.269270 0.206290 0.124173 0.060490 0.023499 0.006901 0.001367 0.000141 0.000004]';

% alternative reproductive functions
r1 = @(Y,A,D,H) (1./(1+exp(-(-1.6874 + 0.0482.*D - 0.0142.*A))))./(1-(1./(1+exp(-(-1.6874 + 0.0482.*D - 0.0142.*A)))));
r2 = @(Y,A,D,H) (1./(1+exp(-(-1.9893 + 0.0268.*D))))            ./(1-(1./(1+exp(-(-1.9893 + 0.0268.*D)))));
% simply use mean r3 = 0.147/(1-0.147) = 0.172 in transition functions

% alternative survival functions
% simply use mean s1 = 0.951 in transition functions
s2 = @(Y,A,D,H) 1./(1+exp(-(2.7382 + 0.0488.*D)));
s3 = @(Y,A,D,H) 1./(1+exp(-(4.2934 + 0.0531.*D - 0.0437.*(Y+A))));

% model-specific harvest functions for each population model(partition total harvest into that of young and adults)
% ensure that harvest of each age class does not exceed the pre-harvest abundance
% young
hy0 = @(Y,A,D,H) min(H.*2.2 .*r1(Y,A,D,H)./(1+2.2 .*r1(Y,A,D,H)), (Y+A).*0.951      .*1.1 .*r1(Y,A,D,H));
hy1 = @(Y,A,D,H) min(H.*2.2 .*r1(Y,A,D,H)./(1+2.2 .*r1(Y,A,D,H)), (Y+A).*s2(Y,A,D,H).*1.1 .*r1(Y,A,D,H));
hy2 = @(Y,A,D,H) min(H.*2.2 .*r1(Y,A,D,H)./(1+2.2 .*r1(Y,A,D,H)), (Y+A).*s3(Y,A,D,H).*1.1 .*r1(Y,A,D,H));
hy3 = @(Y,A,D,H) min(H.*2.2 .*r2(Y,A,D,H)./(1+2.2 .*r2(Y,A,D,H)), (Y+A).*0.951      .*1.1 .*r2(Y,A,D,H));
hy4 = @(Y,A,D,H) min(H.*2.2 .*r2(Y,A,D,H)./(1+2.2 .*r2(Y,A,D,H)), (Y+A).*s2(Y,A,D,H).*1.1 .*r2(Y,A,D,H));
hy5 = @(Y,A,D,H) min(H.*2.2 .*r2(Y,A,D,H)./(1+2.2 .*r2(Y,A,D,H)), (Y+A).*s3(Y,A,D,H).*1.1 .*r2(Y,A,D,H));
hy6 = @(Y,A,D,H) min(H.*2.2 .*0.172      ./(1+2.2 .*0.172),       (Y+A).*0.951      .*1.1 .*0.172);
hy7 = @(Y,A,D,H) min(H.*2.2 .*0.172      ./(1+2.2 .*0.172),       (Y+A).*s2(Y,A,D,H).*1.1 .*0.172);
hy8 = @(Y,A,D,H) min(H.*2.2 .*0.172      ./(1+2.2 .*0.172),       (Y+A).*s3(Y,A,D,H).*1.1 .*0.172);
% adults
ha0	= @(Y,A,D,H) min(H-hy0(Y,A,D,H), (Y+A).*0.951);
ha1	= @(Y,A,D,H) min(H-hy1(Y,A,D,H), (Y+A).*s2(Y,A,D,H));
ha2	= @(Y,A,D,H) min(H-hy2(Y,A,D,H), (Y+A).*s3(Y,A,D,H));
ha3	= @(Y,A,D,H) min(H-hy3(Y,A,D,H), (Y+A).*0.951);
ha4	= @(Y,A,D,H) min(H-hy4(Y,A,D,H), (Y+A).*s2(Y,A,D,H));
ha5	= @(Y,A,D,H) min(H-hy5(Y,A,D,H), (Y+A).*s3(Y,A,D,H));
ha6	= @(Y,A,D,H) min(H-hy6(Y,A,D,H), (Y+A).*0.951);
ha7	= @(Y,A,D,H) min(H-hy7(Y,A,D,H), (Y+A).*s2(Y,A,D,H));
ha8	= @(Y,A,D,H) min(H-hy8(Y,A,D,H), (Y+A).*s3(Y,A,D,H));

% model-specific state transitions
%... young
y0 = @(Y,A,D,H) (Y+A).*r1(Y,A,D,H).*0.951      .*1.1-hy0(Y,A,D,H);
y1 = @(Y,A,D,H) (Y+A).*r1(Y,A,D,H).*s2(Y,A,D,H).*1.1-hy1(Y,A,D,H);
y2 = @(Y,A,D,H) (Y+A).*r1(Y,A,D,H).*s3(Y,A,D,H).*1.1-hy2(Y,A,D,H);
y3 = @(Y,A,D,H) (Y+A).*r2(Y,A,D,H).*0.951      .*1.1-hy3(Y,A,D,H);
y4 = @(Y,A,D,H) (Y+A).*r2(Y,A,D,H).*s2(Y,A,D,H).*1.1-hy4(Y,A,D,H);
y5 = @(Y,A,D,H) (Y+A).*r2(Y,A,D,H).*s3(Y,A,D,H).*1.1-hy5(Y,A,D,H);
y6 = @(Y,A,D,H) (Y+A).*0.172      .*0.951      .*1.1-hy6(Y,A,D,H);
y7 = @(Y,A,D,H) (Y+A).*0.172      .*s2(Y,A,D,H).*1.1-hy7(Y,A,D,H);
y8 = @(Y,A,D,H) (Y+A).*0.172      .*s3(Y,A,D,H).*1.1-hy8(Y,A,D,H);
Ytran = {y0,y1,y2,y3,y4,y5,y6,y7,y8};
%... adults
a0 = @(Y,A,D,H) (Y+A).*0.951      -ha0(Y,A,D,H);
a1 = @(Y,A,D,H) (Y+A).*s2(Y,A,D,H)-ha1(Y,A,D,H);
a2 = @(Y,A,D,H) (Y+A).*s3(Y,A,D,H)-ha2(Y,A,D,H);
a3 = @(Y,A,D,H) (Y+A).*0.951      -ha3(Y,A,D,H);
a4 = @(Y,A,D,H) (Y+A).*s2(Y,A,D,H)-ha4(Y,A,D,H);
a5 = @(Y,A,D,H) (Y+A).*s3(Y,A,D,H)-ha5(Y,A,D,H);
a6 = @(Y,A,D,H) (Y+A).*0.951      -ha6(Y,A,D,H);
a7 = @(Y,A,D,H) (Y+A).*s2(Y,A,D,H)-ha7(Y,A,D,H);
a8 = @(Y,A,D,H) (Y+A).*s3(Y,A,D,H)-ha8(Y,A,D,H);
Atran = {a0,a1,a2,a3,a4,a5,a6,a7,a8};
%... TempDays (single model)
Dtran = @(Y,A,D,H) 1; % for every state/action simply return "1"
                      % next year's days are given by 1 x e3,
                      % with probabilities w3
                      
% COMPUTE REWARDS
%... deterministic next total population size
EN = {a0(X(:,1),X(:,2),X(:,3),X(:,4))+y0(X(:,1),X(:,2),X(:,3),X(:,4))...
      a1(X(:,1),X(:,2),X(:,3),X(:,4))+y1(X(:,1),X(:,2),X(:,3),X(:,4))...
      a2(X(:,1),X(:,2),X(:,3),X(:,4))+y2(X(:,1),X(:,2),X(:,3),X(:,4))...
      a3(X(:,1),X(:,2),X(:,3),X(:,4))+y3(X(:,1),X(:,2),X(:,3),X(:,4))...
      a4(X(:,1),X(:,2),X(:,3),X(:,4))+y4(X(:,1),X(:,2),X(:,3),X(:,4))...
      a5(X(:,1),X(:,2),X(:,3),X(:,4))+y5(X(:,1),X(:,2),X(:,3),X(:,4))...
      a6(X(:,1),X(:,2),X(:,3),X(:,4))+y6(X(:,1),X(:,2),X(:,3),X(:,4))...
      a7(X(:,1),X(:,2),X(:,3),X(:,4))+y7(X(:,1),X(:,2),X(:,3),X(:,4))...
      a8(X(:,1),X(:,2),X(:,3),X(:,4))+y8(X(:,1),X(:,2),X(:,3),X(:,4))};

%... calculate poplation utility
alpha = 10;
utility = @(pop) 1./(1+exp(-(alpha-(abs(pop-60)))));  % pop utility
U = cell(1,9);
%... total utility is harvest x population utility
% U{1} = utility(EN{1}) .*(hy0(X(:,1),X(:,2),X(:,3),X(:,4))+ha0(X(:,1),X(:,2),X(:,3),X(:,4)));
% U{2} = utility(EN{2}) .*(hy1(X(:,1),X(:,2),X(:,3),X(:,4))+ha1(X(:,1),X(:,2),X(:,3),X(:,4)));
% U{3} = utility(EN{3}) .*(hy2(X(:,1),X(:,2),X(:,3),X(:,4))+ha2(X(:,1),X(:,2),X(:,3),X(:,4)));
% U{4} = utility(EN{4}) .*(hy3(X(:,1),X(:,2),X(:,3),X(:,4))+ha3(X(:,1),X(:,2),X(:,3),X(:,4)));
% U{5} = utility(EN{5}) .*(hy4(X(:,1),X(:,2),X(:,3),X(:,4))+ha4(X(:,1),X(:,2),X(:,3),X(:,4)));
% U{6} = utility(EN{6}) .*(hy5(X(:,1),X(:,2),X(:,3),X(:,4))+ha5(X(:,1),X(:,2),X(:,3),X(:,4)));
% U{7} = utility(EN{7}) .*(hy6(X(:,1),X(:,2),X(:,3),X(:,4))+ha6(X(:,1),X(:,2),X(:,3),X(:,4)));
% U{8} = utility(EN{8}) .*(hy7(X(:,1),X(:,2),X(:,3),X(:,4))+ha7(X(:,1),X(:,2),X(:,3),X(:,4)));
% U{9} = utility(EN{9}) .*(hy8(X(:,1),X(:,2),X(:,3),X(:,4))+ha8(X(:,1),X(:,2),X(:,3),X(:,4)));

%... total utility is wp(pop) + (1-wp)(harvest)
wp = 1.0;
U{1} = wp*utility(EN{1}) + (1-wp)*(hy0(X(:,1),X(:,2),X(:,3),X(:,4))+ha0(X(:,1),X(:,2),X(:,3),X(:,4)));
U{2} = wp*utility(EN{2}) + (1-wp)*(hy1(X(:,1),X(:,2),X(:,3),X(:,4))+ha1(X(:,1),X(:,2),X(:,3),X(:,4)));
U{3} = wp*utility(EN{3}) + (1-wp)*(hy2(X(:,1),X(:,2),X(:,3),X(:,4))+ha2(X(:,1),X(:,2),X(:,3),X(:,4)));
U{4} = wp*utility(EN{4}) + (1-wp)*(hy3(X(:,1),X(:,2),X(:,3),X(:,4))+ha3(X(:,1),X(:,2),X(:,3),X(:,4)));
U{5} = wp*utility(EN{5}) + (1-wp)*(hy4(X(:,1),X(:,2),X(:,3),X(:,4))+ha4(X(:,1),X(:,2),X(:,3),X(:,4)));
U{6} = wp*utility(EN{6}) + (1-wp)*(hy5(X(:,1),X(:,2),X(:,3),X(:,4))+ha5(X(:,1),X(:,2),X(:,3),X(:,4)));
U{7} = wp*utility(EN{7}) + (1-wp)*(hy6(X(:,1),X(:,2),X(:,3),X(:,4))+ha6(X(:,1),X(:,2),X(:,3),X(:,4)));
U{8} = wp*utility(EN{8}) + (1-wp)*(hy7(X(:,1),X(:,2),X(:,3),X(:,4))+ha7(X(:,1),X(:,2),X(:,3),X(:,4)));
U{9} = wp*utility(EN{9}) + (1-wp)*(hy8(X(:,1),X(:,2),X(:,3),X(:,4))+ha8(X(:,1),X(:,2),X(:,3),X(:,4)));

% COMPUTE MODEL-SPECIFIC TRANSITIONS
cleanup = 0; % no adjustments to P 
e = rectgrid(e1,e2,e3); w = prod(rectgrid(w1,w2,w3),2);
Stran=cell(1,9);
for i=1:9   % 9 alternative models
  g = @(X,e) [Ytran{i}(X(:,1),X(:,2),X(:,3),X(:,4)).*e(:,1)  Atran{i}(X(:,1),X(:,2),X(:,3),X(:,4)).*e(:,2)  Dtran(X(:,1),X(:,2),X(:,3),X(:,4)).*e(:,3)];
  Stran{i} = g2P(g,{Y,A,D},X,e,w,cleanup);
end;

bigR = U; % rewards (all models)
bigP = Stran; % state transitions (all models)
I = getI(X,1:3);
results=cell(1,9);
lrp=cell(1,9);

% -----------------------------------------------------------
% COMPUTE MODEL-SPECIFIC POLICIES
% -----------------------------------------------------------
%... different ways of handling an infinite time horizon
% INFINITE TIME HORIZON, DELTA = 1; RELATIVE (AVERAGE) VALUE
% model=struct('discount',1.0,'Ix',I);
% options=struct('algorithm','f','print',2,'tol',1e-6,'relval',1);
% 
% INFINITE TIME HORIZON, DISCOUNTED; MODIFIED POLICY ITERATION
%   model=struct('discount',0.90,'Ix',I);
%   options=struct('print',2);
 
% INFINITE TIME HORIZON, function iteration, stopping rule
%  model=struct('discount',1,'Ix',I);
%  options=struct('algorithm','f','print',2,'nochangelim',50);

% INFINITE TIME HORIZON VANISHING DISCOUNT (PREFERRED METHOD)
 model=struct('d',1.0,'Ix',I);
 options=struct('algorithm','i','vanish',0.99999,'print',2);

% % solve for model-specific policies and display
% %... note convergence issue with exponential model M06: S(.)R(.)
 lroptions=struct('fast',1);
 for i=1:9
   model.R=bigR{i};
   model.P=bigP{i};
%    results{i}=mdpsolve(model,options);
%    pstar=results{i}.pstar;
%    pstar(pstar<0)=0;
%    pstar=mxv(pstar,1./sum(pstar,1));
%    lrp{i}=longrunP(pstar,lroptions);
 end
% 
% RS=rectgrid([1;2;3],[1;2;3]);
% figure(1); 
% CM=linspace(0.95,0,64)'*[1 1 1];
% set(gcf,'units','normalized','position',[0.2 0.1 0.65 0.75])
% for i=1:9
%     subplot(3,3,i)
%     ind=X(results{i}.Ixopt,3)==8;
%     policy = X(results{i}.Ixopt,:);
%     patchplot(policy(ind,2),policy(ind,1),policy(ind,4),[0 Hmax]);
%     title(['Model ' num2str(i-1) ' (R' num2str(RS(i,1)) 'S' num2str(RS(i,2)) ')']);
%     xlabel('Adults'); 
%     ylabel('Young'); 
% end
% h=colorbar;
% pos=get(h,'position');
% pos(1)=1-pos(3)*3;
% pos(2)=0.25; pos(4)=0.5;
% set(h,'position',pos)
% set(gcf,'Name','Optimal Actions when TempDays=8')
% 

% -----------------------------------------------------------
% COMPUTE PASSIVELY ADAPTIVE POLICY for specified model weights
% -----------------------------------------------------------
% final weights 2016
% modelwts = [0.02001 0.02272 0.00000 0.24811 0.20749 0.00000 0.21812 0.28355 0.00000]; 
% 2 June 2017 model weights
modelwts = [0.01178 0.02390 0.00000 0.30114 0.40264 0.00000 0.07444 0.18610 0.00000]; 
[v,a,B]=amdppassive2(modelwts,bigP,bigR,model,options);
policy = X(a(:,1),:); % construct policy for model state

return;

% SPECIFY POLICY OUTPUT FILE
%dlmwrite('policy 2017 H50 oldobj vanish 3June2017.dat',policy)

% SIMULATE ANY MODEL-AVERAGED, OPTIMAL POLICY WITH EACH MODEL
% i.e., What if we use a model-averaged policy and yet a single model is most appropriate?
% We could then average simulations across models based on current model
%  weights to get results for the model-averaged policy under an average
%  model 
%... inputs to mdpsim are:
%                      (1) the state index to initialize the simulation
%                          (which rows of the state matrix)
%                      (2) the ns x nx transition model P to use for simulation
%                      (3) the columns of P to use for the policy (a is
%                          defined as a particular model-weighted policy from
%                          amdpassive)
%                      (4) the length of the simulation in years

N0index =  find(ismember(S,[16.0 72.0 4.0],'rows')); % define a starting state
sims = 100;

[states,actions] = mdpsim(N0index,bigP{1},a(:,1),sims);
simresult = X(actions,:);
M0Pop = mean(simresult(:,1)+simresult(:,2));
M0Dif = mean(abs(simresult(:,1)+simresult(:,2)-60));
M0PopSD = std(simresult(:,1)+simresult(:,2)); 
M0Quota = mean(simresult(:,4));
M0QuotaSD = std(simresult(:,4));
 
[states,actions] = mdpsim(N0index,bigP{2},a(:,1),sims);
simresult = X(actions,:);
M1Pop = mean(simresult(:,1)+simresult(:,2));
M1PopSD = std(simresult(:,1)+simresult(:,2)); 
M1Dif = mean(abs(simresult(:,1)+simresult(:,2)-60));
M1Quota = mean(simresult(:,4));
M1QuotaSD = std(simresult(:,4));
 
[states,actions] = mdpsim(N0index,bigP{3},a(:,1),sims);
simresult = X(actions,:);
M2Pop = mean(simresult(:,1)+simresult(:,2));
M2PopSD = std(simresult(:,1)+simresult(:,2)); 
M2Dif = mean(abs(simresult(:,1)+simresult(:,2)-60));
M2Quota = mean(simresult(:,4));
M2QuotaSD = std(simresult(:,4));
  
[states,actions] = mdpsim(N0index,bigP{4},a(:,1),sims);
simresult = X(actions,:);
M3Pop = mean(simresult(:,1)+simresult(:,2));
M3PopSD = std(simresult(:,1)+simresult(:,2)); 
M3Dif = mean(abs(simresult(:,1)+simresult(:,2)-60));
M3Quota = mean(simresult(:,4));
M3QuotaSD = std(simresult(:,4));
  
[states,actions] = mdpsim(N0index,bigP{5},a(:,1),sims);
simresult = X(actions,:);
M4Pop = mean(simresult(:,1)+simresult(:,2));
M4PopSD = std(simresult(:,1)+simresult(:,2)); 
M4Dif = mean(abs(simresult(:,1)+simresult(:,2)-60));
M4Quota = mean(simresult(:,4));
M4QuotaSD = std(simresult(:,4));
  
[states,actions] = mdpsim(N0index,bigP{6},a(:,1),sims);
simresult = X(actions,:);
M5Pop = mean(simresult(:,1)+simresult(:,2));
M5PopSD = std(simresult(:,1)+simresult(:,2)); 
M5Dif = mean(abs(simresult(:,1)+simresult(:,2)-60));
M5Quota = mean(simresult(:,4));
M5QuotaSD = std(simresult(:,4));
 
[states,actions] = mdpsim(N0index,bigP{7},a(:,1),sims);
simresult = X(actions,:);
M6Pop = mean(simresult(:,1)+simresult(:,2));
M6PopSD = std(simresult(:,1)+simresult(:,2)); 
M6Dif = mean(abs(simresult(:,1)+simresult(:,2)-60));
M6Quota = mean(simresult(:,4));
M6QuotaSD = std(simresult(:,4));

[states,actions] = mdpsim(N0index,bigP{8},a(:,1),sims);
simresult = X(actions,:);
M7Pop = mean(simresult(:,1)+simresult(:,2));
M7PopSD = std(simresult(:,1)+simresult(:,2)); 
M7Dif = mean(abs(simresult(:,1)+simresult(:,2)-60));
M7Quota = mean(simresult(:,4));
M7QuotaSD = std(simresult(:,4));

[states,actions] = mdpsim(N0index,bigP{9},a(:,1),sims);
simresult = X(actions,:);
M8Pop = mean(simresult(:,1)+simresult(:,2));
M8PopSD = std(simresult(:,1)+simresult(:,2)); 
M8Dif = mean(abs(simresult(:,1)+simresult(:,2)-60));
M8Quota = mean(simresult(:,4));
M8QuotaSD = std(simresult(:,4));

%Overall, weighted means 
PopMean =   weighted_mean('arithmetic',[M0Pop M1Pop M2Pop M3Pop M4Pop M5Pop M6Pop M7Pop M8Pop],modelwts,'all')
PopSD =     weighted_mean('arithmetic',[M0PopSD M1PopSD M2PopSD M3PopSD M4PopSD M5PopSD M6PopSD M7PopSD M8PopSD],modelwts,'all')
DifMean =   weighted_mean('arithmetic',[M0Dif M1Dif M2Dif M3Dif M4Dif M5Dif M6Dif M7Dif M8Dif],modelwts,'all')
QuotaMean = weighted_mean('arithmetic',[M0Quota M1Quota M2Quota M3Quota M4Quota M5Quota M6Quota M7Quota M8Quota],modelwts,'all')
QuotaSD =   weighted_mean('arithmetic',[M0QuotaSD M1QuotaSD M2QuotaSD M3QuotaSD M4QuotaSD M5QuotaSD M6QuotaSD M7QuotaSD M8QuotaSD],modelwts,'all')
 
 
 
 
 
 
 
 
 
