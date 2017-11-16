function [Xorig,A,Dpartial] = genrandprob(params)
% genrandprob       Generate a random instance of the Sensor Network
%                   Localization problem in dimension r.
%
% [Xorig,A,Dpartial,stats] = genrandprob(params) generates the original
% sensor positions Xorig, the anchor positions A, and the partial distance
% (squared) matrix Dpartial where entry (i,j) is ||pi-pj||^2 if the
% distance between sensor i and senor j, ||pi-pj||, is less than the radio
% range R; otherwise Dpartial(i,j) = 0; Note that the matrix [Xorig; A] has
% n rows and r columns.
%
% The sensors and anchors are uniformly distributed in the [0,1]^r region.
%
% Optional parameters:
%
%        nf :   the noise factor applied to the exact distances to obtain
%               a noisy partial distance matrix.  The noise is not applied
%               to the distances between the anchors
%
%      seed :   an integer between 0 and 2^32 which sets the default
%               random number stream


% Extract parameters
n    = params.n;
m    = params.m;
r    = params.r;
R    = params.R;
nf   = params.nf;
seed = params.seed;

% Seed the rng
if ~isempty(seed)
    rng(seed);
end

% generate random points in dim r
Porig = rand(n,r) - 0.5;
Xorig = Porig(1:n-m,:);
A = Porig(n-m+1:end,:);

Dpartial = formDpartial(Xorig, A, R, nf);

end
