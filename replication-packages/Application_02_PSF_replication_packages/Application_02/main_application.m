clear all
close all
addpath('Functions')
addpath('Scripts')     
addpath('Results')
%% Load data

% cost frontier: sheet1 has w3 as numeraire
Data = readtable("Data_Application.xlsx", 'Sheet','Sheet1', ReadVariableNames=true);
time = reshape(Data.t,80,466);

% summary statistics
avex = [mean((Data.cost)), mean((Data.x1)), mean((Data.x2)), mean((Data.y1)), mean((Data.y2)), mean((Data.y3))];
stdx = [sqrt(var((Data.cost))), sqrt(var((Data.x1))), sqrt(var((Data.x2))), sqrt(var((Data.y1))), sqrt(var((Data.y2))), sqrt(var((Data.y3)))];
minx = [min((Data.cost)), min((Data.x1)), min((Data.x2)), min((Data.y1)), min((Data.y2)), min((Data.y3))];
maxx = [max((Data.cost)), max((Data.x1)), max((Data.x2)), max((Data.y1)), max((Data.y2)), max((Data.y3))];

logavex = [mean(log(Data.cost)), mean(log(Data.x1)), mean(log(Data.x2)), mean(log(Data.y1)), mean(log(Data.y2)), mean(log(Data.y3))]';
logstdx = [sqrt(var(log(Data.cost))), sqrt(var(log(Data.x1))), sqrt(var(log(Data.x2))), sqrt(var(log(Data.y1))), sqrt(var(log(Data.y2))), sqrt(var(log(Data.y3)))]';
logminx = [min(log(Data.cost)), min(log(Data.x1)), min(log(Data.x2)), min(log(Data.y1)), min(log(Data.y2)), min(log(Data.y3))]';
logmaxx = [max(log(Data.cost)), max(log(Data.x1)), max(log(Data.x2)), max(log(Data.y1)), max(log(Data.y2)), max(log(Data.y3))]';


fprintf('\n summary statistics \n');
disp([logavex, logstdx, logminx, logmaxx])

%% Set structure Params
Params.T = size(time,1);       % # time
Params.N = size(time,2);      % # firms  
Params.m = floor(Params.T^0.2); % degree of approximation of parameters in sieves
Params.p = 5;                  % # regressors:  5     
Params.K_group = 4;

Params.lambda_s3 = (Params.N*Params.T)^(1/2)*log(Params.N*Params.T)/2;  % for step 3

% inital points of sigu, c, tau
Params.c1    = 1;   
Params.sigu1 = 3;   
Params.c2    = 1;
Params.sigu2 = 3;
Params.tau   = 0.5;

% initial points
%pi0 = ones(17,1);       % 17 = (mLbar-1 for intercept, xBmLbar)
Params.omegabar0 = [Params.c1, sqrt(Params.sigu1)];
Params.omegabar0_prime = [Params.c1, sqrt(Params.sigu1), Params.c2, sqrt(Params.sigu2), Params.tau];

% bounds
Params.LB = [-100, 0];                      % lower bound for parameter search region for omegabar
Params.UB = [100, 100];                     % upper bound for parameter search region for omegabar
Params.LB_prime = [-100, 0, -100, 0, 0];    % lower bound for parameter search region for omegabar_prime
Params.UB_prime = [100, 200, 100, 200, 1];  % upper bound for parameter search region for omegabar_prime

% set size of step 3 m, or m lower bar
% if = NaN, it will vary by group,
% if = 2, it will set step 3 to =2 for every group
Params.mLbar_S3 = NaN;

% Set parameters used 
T = Params.T;
N = Params.N;
p = Params.p;
m = Params.m;
K_group = Params.K_group;
lambda_s3 = Params.lambda_s3;
mLbar_S3 = Params.mLbar_S3;

%% Set data used
% Set variables
y = log(reshape(Data.cost,T,N))';
x1 = log(reshape(Data.x1,T,N))';
x2 = log(reshape(Data.x2,T,N))';
x3 = log(reshape(Data.y1,T,N))';
x4 = log(reshape(Data.y2,T,N))';
x5 = log(reshape(Data.y3,T,N))';

% Demean and standardize data
x1_cl = (x1 - mean(x1,2)) ./ mean(sqrt(var(x1,0,2)));
x2_cl = (x2 - mean(x2,2)) ./ mean(sqrt(var(x2,0,2)));
x3_cl = (x3 - mean(x3,2)) ./ mean(sqrt(var(x3,0,2)));
x4_cl = (x4 - mean(x4,2)) ./ mean(sqrt(var(x4,0,2)));
x5_cl = (x5 - mean(x5,2)) ./ mean(sqrt(var(x5,0,2)));

%% Algorithm -- run the estimation (see Functions/Estimation.m)
% All Step 1-5' estimation lives in Estimation.m. It returns every variable
% it creates in the struct Est; unpack Est into this workspace so the scripts
% below find the variables they expect.

Est = Estimation(Params, y, x1, x2, x3, x4, x5, x1_cl, x2_cl, x3_cl, x4_cl, x5_cl);

est_fields = fieldnames(Est);
for kk = 1:numel(est_fields)
    eval([est_fields{kk} ' = Est.' est_fields{kk} ';']);
end
clear Est est_fields kk

%% run scripts in this order
run getFrontiers
run Asymptotics
run Frontiers_SE
run Figures_Application
run sigv_omegabar_estimates

