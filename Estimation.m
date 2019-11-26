%% ########## Models estimation ##########
clear
clc
load('Data_pret30.mat')
load('Data_pret70.mat')
%% Feedback HMM (wit AFO effect on transitions)
Data = Data_pret30;
est = [];

for subj = 1:21
    
D = Data(subj);

% VM concentration
afo_vec = D.AFO;
afo_vec = afo_vec(:);
afo_neg = afo_vec(afo_vec < 0)*10;
afo_pos = afo_vec(afo_vec > 0)*10;

nchains = 10;
opts = optimoptions(@fmincon,'Algorithm','sqp'); % alternative: 'interior-point'
problem = createOptimProblem('fmincon','objective',...
    @(par)(truVM_lik(par,afo_pos,0,pi)),...
    'x0',1,'lb',0.01,'ub',200,'options',opts);
ms = MultiStart('UseParallel',true,'Display','iter');
[pars,fval,eflag,output,manymins] = run(ms,problem,nchains);

var(1) = pars;

nchains = 10;
opts = optimoptions(@fmincon,'Algorithm','sqp'); % alternative: 'interior-point'
problem = createOptimProblem('fmincon','objective',...
    @(par)(truVM_lik(par,afo_neg,-pi,0)),...
    'x0',1,'lb',0.01,'ub',200,'options',opts);
ms = MultiStart('UseParallel',true,'Display','iter');
[pars,fval,eflag,output,manymins] = run(ms,problem,nchains);

var(2) = pars;

fixpar = est30_newpar(subj,:);

nchains = 20;
low_b = [0.001 0.001];
up_b =  [0.999 0.999];
starting_points = unifrnd(low_b,up_b,1,length(low_b));
opts = optimoptions(@fmincon,'Algorithm','sqp'); % alternative: 'interior-point'
problem = createOptimProblem('fmincon','objective',...
    @(par)(tanhmodel(par,var,fixpar,D)),...
    'x0',starting_points,'lb',low_b,'ub',up_b,'options',opts);
ms = MultiStart('UseParallel',true,'Display','iter');
[pars,fval,eflag,output,manymins] = run(ms,problem,nchains);

est(subj,:) = pars;

end
%% RW model (without AFO effect on transitions)
Data = Data_pret70;
est = [];

for subj = 1:21
    
D = Data(subj);

% VM concentration

afo_vec = D.AFO;
afo_vec = afo_vec(:);
afo_neg = afo_vec(afo_vec < 0)*10;
afo_pos = afo_vec(afo_vec > 0)*10;

nchains = 10;
opts = optimoptions(@fmincon,'Algorithm','sqp'); % alternative: 'interior-point'
problem = createOptimProblem('fmincon','objective',...
    @(par)(truVM_lik(par,afo_pos,0,pi)),...
    'x0',1,'lb',0.01,'ub',200,'options',opts);
ms = MultiStart('UseParallel',true,'Display','iter');
[pars,fval,eflag,output,manymins] = run(ms,problem,nchains);

var(1) = pars;

nchains = 10;
opts = optimoptions(@fmincon,'Algorithm','sqp'); % alternative: 'interior-point'
problem = createOptimProblem('fmincon','objective',...
    @(par)(truVM_lik(par,afo_neg,-pi,0)),...
    'x0',1,'lb',0.01,'ub',200,'options',opts);
ms = MultiStart('UseParallel',true,'Display','iter');
[pars,fval,eflag,output,manymins] = run(ms,problem,nchains);

var(2) = pars;

% non-homogeneous HMM model

nchains = 10;
low_b = [0.001 0.001 0.001];
up_b =  [0.999 0.999 0.999];
starting_points = unifrnd(low_b,up_b,1,length(low_b));
opts = optimoptions(@fmincon,'Algorithm','sqp'); % alternative: 'interior-point'
problem = createOptimProblem('fmincon','objective',...
    @(par)(RWmodel(par,var,D)),...
    'x0',starting_points,'lb',low_b,'ub',up_b,'options',opts);
ms = MultiStart('UseParallel',true,'Display','iter');
[pars,fval,eflag,output,manymins] = run(ms,problem,nchains);

est(subj,:) = [var(1) var(2) pars];

end
