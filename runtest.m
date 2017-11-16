% Specify parameters and options

% n is the total number of sensors and anchors
% m is the number of anchors
params.n    = 1000;
params.m    = 10;
params.r    = 2;
params.R    = 0.25;
params.nf   = 0.10;
params.seed = 2017;  % use 'shuffle' or [] for no rng seeding

opts.verbose       = 1;
opts.weighted_ls   = 0;
opts.expvecweights = 0;
opts.refine        = 1;
opts.MaxCliqueSize = 50;

opts.plotsoln      = 1;

if opts.verbose
    disp('-------------------------------------------------');
    disp(params);
end

% Generate a random EDM with embdim r
tt = tic;
[Xorig, A, Dpartial] = genrandprob(params);
if opts.verbose
    fprintf('%16s: %6.2f s\n\n', 'Generate problem', toc(tt));
end

% Solve EDM problem using Exposing Vector Algorithm
[X, Xrefined, xinds, info] = ExpVecEDM(Dpartial, A, params.r, opts);

% Compute RMSD
if isempty(A)
    [info.rmsd.X, X] = procrustes_noncenter(Xorig(xinds,:), X);
    info.rmsd.X = info.rmsd.X/params.R*100;
    if opts.refine
        [info.rmsd.Xrefined, Xrefined] = ...
            procrustes_noncenter(Xorig(xinds,:), Xrefined);
        info.rmsd.Xrefined = info.rmsd.Xrefined/params.R*100;
    end
else
    info.rmsd.X = rmsd(Xorig(xinds,:), X);
    info.rmsd.X = info.rmsd.X/params.R*100;
    if opts.refine
        info.rmsd.Xrefined = rmsd(Xorig(xinds,:), Xrefined);
        info.rmsd.Xrefined = info.rmsd.Xrefined/params.R*100;
    end
end

% Output info
if opts.verbose
    disp('info.eig:');     disp(info.eig);
    disp('info.rmsd:');    disp(info.rmsd)
    disp('info.resnorm:'); disp(info.resnorm)
end

% Plot solution
if opts.plotsoln
    plotsoln(X, Xrefined, Xorig(xinds,:), A);
end
