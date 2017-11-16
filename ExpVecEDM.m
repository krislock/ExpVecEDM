function [X, Xrefined, xinds, info] = ExpVecEDM(Dpartial, A, r, opts)
% ExpVecEDM     Solve the Euclidean Distance Matrix completion problem.
%
% X = ExpVecEDM(Dpartial,A,r) solves the Euclidean Distance Matrix
% completion problem with partial distance (squared) matrix Dpartial, and 
% anchor positions A, in dimension r.
%
% Set A=[] to solve the problem without anchors.
% 
% NOTE:  When A~=[], the distances between the anchors must be in the
% bottom right corner of the matrix Dpartial, and A must have r columns.
%
% X = ExpVecEDM(Dpartial,A,r,OPTS) replaces the default parameters by those
% in the structure OPTS:
%
%    OPTS.verbose       = {[0],1} where 0 is no output and 1 prints output 
%                         of the progress of the algorithm.
%
%    OPTS.weighted_ls   = {[0],1} where 1 uses the inverse squared 
%                         distances as weights when solving the least 
%                         squares problem.
%
%    OPTS.expvecweights = {[0],1} where 1 uses special weights when forming
%                         the exposing vector Y as a linear combination of
%                         the exposing vectors for each face
%
%    OPTS.refine        = {0,[1]} where 1 refines the solution X using
%                         gradient descent on an unconstrained nonlinear
%                         optimization problem to obtain Xrefined.
%
%    OPTS.MaxCliqueSize = {integer > 0, [50]} specifies the maximum size of
%                         the cliques used in the algorithm.
%
% [X,Xrefined] = SNLSDPclique(...) returns the refined positions Xrefined.
%
% [X,Xrefined,xinds] = SNLSDPclique(...) returns the indices of the points
% whose positions are returned in X. NOTE: nodes having degree < r+1 are
% removed before problem is solved since their location cannot be
% determined; a warning message is printed when this happens.
%
% [X,Xrefined,xinds,INFO] = SNLSDPclique(...) returns a structure INFO with
%
%    INFO.graph   = information about the graph specified by Dpartial.
%
%    INFO.time    = runnig times for different subroutines.
%
%    INFO.eig     = eigenvalues of the exposing vector Y = UZU' and Z.
%
%    INFO.resnorm = norm(residual)/sqrt(numEdges)
%
% =====================
% Citation information:
% =====================
% Dmitriy Drusvyatskiy, Nathan Krislock, Yuen-Lam Voronin, and Henry
% Wolkowicz. Noisy Euclidean distance realization: Robust facial reduction
% and the pareto frontier. SIAM Journal on Optimization, 27(4):2301-2331,
% 2017.

% ExpVecEDM, version 0.1
% Copyright (C) 2017 by Dmitriy Drusvyatskiy, Nathan Krislock, Yuen-Lam
% Voronin, and Henry Wolkowicz.
% Last Modified 2017 Nov 14

fmt = '%16s: %6.2f s\n';
divider = '==========================\n';

% Start timer
t0 = tic();

if nargin == 3
    opts.verbose       = 0;
    opts.weighted_ls   = 0;
    opts.expvecweights = 0;
    opts.refine        = 1;
    opts.MaxCliqueSize = 50;
end

% Compute graph stats and remove nodes that have degree <= r
v = full(sum(Dpartial>0));

info.graph.mindeg = min(v);
info.graph.maxdeg = max(v);
info.graph.avgdeg = sum(v)/2/length(Dpartial);

if opts.verbose
    disp(info.graph);
end

inds = find(v > r);
D = Dpartial(inds,inds);

m = size(A, 1);
n = size(D, 1);

ainds = length(Dpartial)-m+1:length(Dpartial);
xinds = setdiff(inds, ainds);

if info.graph.mindeg <= r
    warning('Minimum degree of graph <= r');
    fprintf('Removed %d nodes that have degree <= r\n', ...
        length(Dpartial) - n);
end

% Find cliques
tt = tic();
cliques = GrowCliques(D, ainds, opts);
info.time.GrowCliques = toc(tt);
if opts.verbose
    fprintf(divider);
    fprintf(fmt, 'Grow cliques', info.time.GrowCliques);
end

% Remove cliques having length <= r+1
cliques = cliques(cellfun(@length, cliques) > r+1);

% Compute exposing vector Y
tt = tic();
Y = ExposingVector(D, cliques, r, opts);
info.time.expvec = toc(tt);
if opts.verbose
    fprintf(fmt, 'Exposing vector', info.time.expvec);
end

% Compute face
tt = tic();
[U, info.eig.Y] = face(Y, r);
info.time.face = toc(tt);
if opts.verbose
    fprintf(fmt, 'Compute face', info.time.face);
end

% Solve least-squares problem
tt = tic();
[Z, info.eig.Z] = FacialLeastSquares(D, U, opts);
info.time.least_squares = toc(tt);
if opts.verbose
    fprintf(fmt, 'Least-squares', info.time.least_squares);
end

% Compute sensor locations
P = U*sqrtm(Z);

% Procrustes
if isempty(A)
    X = P;
else
    X = procrustes_anchor(A, P);
end

info.time.pre_refinement = toc(t0);
if opts.verbose
    fprintf(divider);
    fprintf(fmt, 'Pre-refinement', info.time.pre_refinement);
end

% Refine positions
if opts.refine
    tt = tic;
    Xrefined = SNLSDP_0.refinepositions(X', A', sqrt(D(1:n-m,:)))';
    info.time.refine = toc(tt);
    if opts.verbose
        fprintf(fmt, 'Refine positions', info.time.refine);
    end
else
    Xrefined = [];
end

info.time.total = toc(t0);

if opts.verbose
    fprintf(divider);
    fprintf(fmt, 'Total', info.time.total);
    fprintf(divider);
    fprintf('\n');
end

% Compute the norm of the residual
info.resnorm.X = residualnorm([X; A], D);
if opts.refine
    info.resnorm.Xrefined = residualnorm([Xrefined; A], D);
end

end