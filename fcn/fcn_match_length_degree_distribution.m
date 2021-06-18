function [A,W] = fcn_match_length_degree_distribution(A,Dist,nbins,nswap,verbose)
% [A,W] = fcn_match_length_degree_distribution(A,Dist,nbins,nswap) takes a
% weighted, symmetric connectivity matrix, A, and Euclidean/fiber length
% matrix, Dist, and generates a randomized network with:
% 1) exactly the same degree sequence
% 2) approximately the same edge length distribution
% 3) exactly the same edge weight distribution
% 4) approximately the same weight-length relationship
%
% Inputs:
%   A,      weighted/symmetric connectivity matrix
%   Dist,   symmetric distance matrix
%   nbins,  number of distance bins (edge length matrix is performed by
%           swapping connections in the same bin)
%   nswap,  total number of edge swaps to perform
%   
% Outputs:
%   A,      binary rewired matrix
%   W,      weighted rewired matrix
%
%
%   Richard Betzel, University of Pennsylvania 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 4
    verbose = false;
end
mask = A ~= 0;                          % nonzero elements
mask = triu(mask,1);                    % upper triangle
Weights = A(mask);                      % values of edge weights
Distances = Dist(mask);                 % values of edge lengths
[~,Jdx] = sort(Distances,'ascend');     % sort lengths smallest to largest

N = length(A);                          % number of nodes
bins = linspace(...                     % length/distance bins
    min(nonzeros(Dist)),...
    max(nonzeros(Dist)),...
    nbins + 1);
bins(end) = bins(end) + 1;
B = zeros(N,N,nbins);                   % 3D stack of distance bins
for i = 1:nbins
    B(:,:,i) = (Dist >= bins(i) & Dist < bins(i + 1))*i;
end
Bsum = sum(B,3);                        % matrix of distance bins

[u,v,w] = find(triu((A ~= 0).*Bsum,1)); % indices and bins of edges
m = length(u);                          % number of edges
iswap = 0;                              % swap counter
while iswap < nswap
    
    e1 = randi(m);                      % choose edge at random
    a = u(e1);
    b = v(e1);
    w1 = w(e1);
    
    indkeep = ...                       % make a list of possible swaps
        u ~= a & u ~= b & v ~= a & v ~= b;
    c = u(indkeep);                     % stubs c
    d = v(indkeep);                     % stubs d
    
    w2poss = w(indkeep);                % the bins of possible swaps
    
    x = (a - 1)*N + c;                  % edge indices
    wx = Bsum(x);                       % bins
    y = (b - 1)*N + d;                  % other set of edge indices
    wy = Bsum(y);                       % bins
    
    idx1 = (w1 == wx) & (w2poss == wy); % get good list indices
    idx2 = (w1 == wy) & (w2poss == wx); % other set of good indices
    idx12 = idx1 | idx2;                % full set
    
    c = c(idx12);                       % update subs
    d = d(idx12);
    
    idxnew = (b - 1)*N + d;             % update edge indices
    jdxnew = (a - 1)*N + c;
    
    Ai = A(idxnew);                     % get true/false on existing edges
    Aj = A(jdxnew);
    
    ind = find(Ai == 0 & Aj == 0);      % find missing edges
    if ~isempty(ind)
        
        r = ind(randi(length(ind)));    % choose random swap
        
        c = c(r);
        d = d(r);
        
        A(a,b) = 0;
        A(b,a) = 0;
        
        A(c,d) = 0;
        A(d,c) = 0;
        
        A(a,c) = 1;
        A(c,a) = 1;
        
        A(b,d) = 1;
        A(d,b) = 1;
        
        e2 = find(indkeep);
        e2 = e2(idx12);
        e2 = e2(r);
        
        u(e1) = min(a,c);
        v(e1) = max(a,c);
        
        u(e2) = min(b,d);
        v(e2) = max(b,d);
        
        w(e1) = Bsum(a,c);
        w(e2) = Bsum(b,d);
        
        iswap = iswap + 1;
        if mod(iswap,1000) == 0 && verbose
            disp(iswap)
        end
        
    end
    
end

idx = find(triu(A,1));              % get edge indices
d = Dist(idx);                      % get distances
[~,jdx] = sort(d,'ascend');         % sort distances smallest to largest
W = zeros(N);                       % output matrix
W(idx(jdx)) = Weights(Jdx);         % add weights