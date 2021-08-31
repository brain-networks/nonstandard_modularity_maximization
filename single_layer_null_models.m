clearvars
close all
clc

addpath('./fcn/')

% visualize the modularity matrix?
viz = 1 ;

%% load data

load ./data/sc_fc_data_singlesubject.mat
d = squareform(pdist(coor));
n = length(sc);

%% example 1. newman girvan

% degree-/strength-preserving null model but allows for self-loops

deg = sum(sc,2);
twom = sum(deg);
p = (deg*deg')/twom;
b = sc - p;

if viz
    imagesc(b)
    title('degree-/strength-preserving null modularity matrix')
    waitforbuttonpress
    close
end

% reference:
%   
%   Newman, M. E., & Girvan, M. (2004). Finding and evaluating community 
%   structure in networks. Physical review E, 69(2), 026113.

%% example 2. newman girvan (no loops)

% degree-/strength-preserving null model but excludes self-loops. (NOTE:
% this function estimates the expected weight by sampling individual null
% connectivity matrices. it is time consuming.

numrand = 1000;
numiter = 32;
deg = sum(sc,2);
temp = 25;
dec = 0.995;
energyfcn = 'maxpercentchange';
tolerance = 0.05;
maxiter = 1e5;

wr = zeros(n,n,numrand);
for i = 1:numrand
    ar = randmio_und(sc,numiter);
    wr(:,:,i) = fcn_match_strength(ar,nonzeros(triu(sc,1)),deg,temp,dec,energyfcn,tolerance,maxiter);
    subplot(1,2,1);
    imagesc(nanmean(wr(:,:,1:i),3));
    subplot(1,2,2);
    plot(deg,sum(nanmean(wr(:,:,1:i),3),2),'k.');
    drawnow;
end
p = nanmean(wr,3);
b = sc - p;

if viz
    imagesc(b)
    title('newman girvan (no loops) null modularity matrix')
    waitforbuttonpress
    close
end

%% example 3. uniform

% all connections (except for self-loops) have an expected weight equal to
% that of ``gamma''

gamma = 0.01;
p = ~eye(n)*gamma;
b = sc - p;

if viz
    imagesc(b)
    title('uniform null modularity matrix')
    waitforbuttonpress
    close
end

% reference(s):
%   
%   Bazzi, M., Porter, M. A., Williams, S., McDonald, M., Fenn, D. J., & 
%   Howison, S. D. (2016). Community detection in temporal multilayer 
%   networks, with an application to correlation networks. 
%   Multiscale Modeling & Simulation, 14(1), 1-41.
%
%   Traag, V. A., Van Dooren, P., & Nesterov, Y. (2011). Narrow scope for 
%   resolution-limit-free community detection. Physical Review E, 84(1), 
%   016114.

%% example 4. geometric (exponential and probabilistic)

% binary connections are drawn stochastically from a decaying exponential
% distribution. weights are then added inversely proportional to distance.

eta = 0.01;
numiter = 100;

p = zeros(n,n,numiter);
for iter = 1:numiter
    edges = find(triu(ones(n),1) > 0);
    prob = exp(-eta*d(edges));
    m = nnz(sc)/2;
    b = zeros(n);
    for i = 1:m
        c = [0; cumsum(prob)];
        r = sum(c(end)*rand > c);
        b(edges(r)) = 1;
        prob(r) = 0;
    end
    w = sort(sc(triu(sc > 0)),'ascend');
    idx = find(b);
    di = d(idx);
    [~,jdx] = sort(di,'descend');
    q = zeros(n);
    q(idx(jdx)) = w;
    q = q + q';
    p(:,:,iter) = q;
end
p = nanmean(p,3);
b = sc - p;

if viz
    imagesc(b)
    title('geometric null modularity matrix')
    waitforbuttonpress
    close
end

% reference(s):
%   
%   Betzel, R. F., Avena-Koenigsberger, A., Goñi, J., He, Y., De Reus, M. 
%   A., Griffa, A., ... & Sporns, O. (2016). Generative models of the human
%   connectome. Neuroimage, 124, 1054-1064.

%% example 5. minimum wiring

% binary connections are placed in the minimum wiring cost configuration. 
% weights are then added inversely proportional to distance.

m = nnz(sc)/2;
edges = find(triu(ones(n),1) > 0);
di = d(edges);
[~,ddx] = sort(di,'ascend');
p = zeros(n);
w = sort(sc(triu(sc > 0)),'descend');
p(edges(ddx(1:m))) = w;
p = p + p';
b = sc - p;

if viz
    imagesc(b)
    title('minimum wiring null modularity matrix')
    waitforbuttonpress
    close
end

% reference(s):
%   
%   Samu, D., Seth, A. K., & Nowotny, T. (2014). Influence of wiring cost 
%   on the large-scale architecture of human cortical connectivity. PLoS 
%   Comput Biol, 10(4), e1003557.

%% example 6. degree + space

% a degree-preserving null model in which edge swaps are restricted to
% different classes based on their distances.

numiter = 100;
nbins = 31;
nswap = 1e4;
q = zeros(n,n,numiter);
for iter = 1:numiter
    [~,ar] = fcn_match_length_degree_distribution(sc,d,nbins,nswap,true);
    q(:,:,iter) = ar + ar';
end
p = nanmean(q,3);
b = sc - p;

if viz
    imagesc(b)
    title('degree + space null modularity matrix')
    waitforbuttonpress
    close
end

% reference(s):
%   
%   Betzel, R. F., & Bassett, D. S. (2018). Specificity and robustness of 
%   long-distance connections in weighted, interareal connectomes. 
%   Proceedings of the National Academy of Sciences, 115(21), E4880-E4889.

%% example 7. binary null model

% preserves the exact binary structure of the observed network but assigns
% each edge the mean weight across all edges in the observed network.

avg = mean(nonzeros(sc));
p = (sc > 0)*avg;
b = sc - p;

if viz
    imagesc(b)
    title('binary null modularity matrix')
    waitforbuttonpress
    close
end

% reference(s):
%   
%   Bassett, D. S., Owens, E. T., Porter, M. A., Manning, M. L., & Daniels, 
%   K. E. (2015). Extraction of force-chain network architecture in 
%   granular materials using community detection. Soft Matter, 11(14), 
%   2731-2744.

%% example 8. signed and weighted (equal weighting of pos/neg contributions)

% handles positive and negative weights separately using newman-girvan null
% model. combines the positive and negative modularity matrices, weighting
% each equally.

fc = fc.*~eye(n);
fcpos = fc.*(fc > 0);
fcneg = -fc.*(fc < 0);
kpos = sum(fcpos,2);
twompos = sum(kpos);
kneg = sum(fcneg,2);
twomneg = sum(kneg);
ppos = (kpos*kpos')/twompos;
pneg = (kneg*kneg')/twomneg;
% wpos = twompos/(twompos + twomneg);
% wneg = twomneg/(twompos + twomneg);
bpos = (fcpos - ppos);
bneg = (fcneg - pneg);

b = bpos/(twompos + twomneg) - bneg/(twompos + twomneg);

if viz
    imagesc(b)
    title('signed and weighted null modularity matrix')
    waitforbuttonpress
    close
end

% reference(s):
% 
%   Gómez, S., Jensen, P., & Arenas, A. (2009). Analysis of community 
%   structure in networks of correlated data. Physical Review E, 80(1), 
%   016114.

%% example 9. signed and weighted (disproportionate weighting of pos/neg contributions)

% handles positive and negative weights separately using newman-girvan null
% model. combines the positive and negative modularity matrices, but
% weights the contribution from positive weights more heavily than that of
% negative weights.

fc = fc.*~eye(n);
fcpos = fc.*(fc > 0);
fcneg = -fc.*(fc < 0);
kpos = sum(fcpos,2);
twompos = sum(kpos);
kneg = sum(fcneg,2);
twomneg = sum(kneg);
ppos = (kpos*kpos')/twompos;
pneg = (kneg*kneg')/twomneg;
wpos = twompos/(twompos + twomneg);
wneg = twomneg/(twompos + twomneg);
bpos = (fcpos - ppos);
bneg = (fcneg - pneg);

b = bpos/(twompos) - bneg/(twompos + twomneg);

if viz
    imagesc(b)
    title('signed and weighted (unequal wei) null modularity matrix')
    waitforbuttonpress
    close
end

% reference(s):
% 
%   Rubinov, M., & Sporns, O. (2011). Weight-conserving characterization of
%   complex functional brain networks. Neuroimage, 56(4), 2068-2079.