clear all
close all
clc

%% example 1. time-varying fc (using with uniform null model)

% load data
load ../data/fmri_ts_data_singlesubject.mat
[numframes,n] = size(ts);

% window size for sliding window analysis
windowsize = 46;

% total number of non-overlapping windows
numwindows = numframes/windowsize;

% structural resolution parameter
gamma = 0.1;

% interlayer coupling parameter
omega = 0.5; 

b = zeros(n*numwindows);
for iwindow = 1:numwindows
    idx = (1:windowsize) + (iwindow - 1)*windowsize;
    rho = corr(ts(idx,:));
    
    jdx = (1:n) + (iwindow - 1)*n;
    b(jdx,jdx) = rho - gamma;
    b(jdx,jdx) = b(jdx,jdx).*~eye(n);
end
b = b + omega*spdiags(ones(n*numwindows,2),[-n,n],n*numwindows,n*numwindows);

%% example 2. multisubject categorical (using fc with uniform null model)

% load data
load ../data/sc_fc_data_multisubject.mat

% structural resolution parameter
gamma = 0.1;

% interlayer coupling parameter
omega = 0.1;

[n,~,nsub] = size(fc);
b = zeros(n*nsub);
for isub = 1:nsub
    idx = (1:n) + (isub - 1)*n;
    b(idx,idx) = fc(:,:,isub) - gamma;
    b(idx,idx) = b(idx,idx).*~eye(n);
end
all2all = n*[(-nsub+1):-1,1:(nsub-1)];
b = b + omega*spdiags(ones(n*nsub,2*nsub-2),all2all,n*nsub,n*nsub);

%% example 3. multisubject/multimodal coupling

% load data
load ../data/sc_fc_data_multisubject.mat

% structural resolution parameter
gamma_fc = 0.3;
gamma_sc = 0;

% interlayer coupling parameter
omega = 1;

% structure-function coupling
eta = 0.05;

[n,~,nsub] = size(fc);
bfc = zeros(n*nsub);
bsc = zeros(n*nsub);
for isub = 1:nsub
    idx = (1:n) + (isub - 1)*n;
    bfc(idx,idx) = fc(:,:,isub) - gamma_fc;
    bfc(idx,idx) = bfc(idx,idx).*~eye(n);
    
    bsc(idx,idx) = corr(sc(:,:,isub)) - gamma_sc;
    bsc(idx,idx) = bsc(idx,idx).*~eye(n);
end
all2all = n*[(-nsub+1):-1,1:(nsub-1)];
bsc = bsc + omega*spdiags(ones(n*nsub,2*nsub-2),all2all,n*nsub,n*nsub);
bfc = bfc + omega*spdiags(ones(n*nsub,2*nsub-2),all2all,n*nsub,n*nsub);
b = [bsc, eta*eye(nsub*n); eta*eye(nsub*n), bfc];

