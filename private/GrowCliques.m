function cliques = GrowCliques(Dpartial, ainds, opts)
% GrowCliques Find a clique around each node.
%
%    cliques = GrowCliques(Dpartial, opts) finds a clique around each node
%    of the graph specified by the partial distance matrix Dpartial.
%    Returns the cell array cliques that contains the index sets of the
%    cliques.

% ExpVecEDM, version 0.1
% Copyright (C) 2016 by Dmitriy Drusvyatskiy, Nathan Krislock,
% Yuen Lam Voronin, and Henry Wolkowicz.

n = length(Dpartial);
MaxCliqueSize = opts.MaxCliqueSize;

assert(MaxCliqueSize > 0, 'opts.MaxCliqueSize must be a positive integer');

II = zeros(n*MaxCliqueSize + length(ainds), 1);
JJ = zeros(n*MaxCliqueSize + length(ainds), 1);
k = 0;
for jc = 1:n % For each clique jc
    
    % Add node jc to clique jc
    k = k + 1;  II(k) = jc;  JJ(k) = jc;
    NodesToAdd = MaxCliqueSize - 1;
    
    q = Dpartial(:, jc) > 0;
    while NodesToAdd && nnz(q)
        ic = find(q);
        %idx = ceil(length(ic)/2);
        ic = ic(end);
        % Add node ic to clique jc
        k = k + 1;  II(k) = ic;  JJ(k) = jc;
        NodesToAdd = NodesToAdd - 1;
        q = q & Dpartial(:, ic(end));
    end
    
end

% Add the anchor clique
for ii = 1:length(ainds)
    ic = ainds(ii);
    k = k + 1;  II(k) = ic;  JJ(k) = n + 1;
end

% Form the clique matrix Dcq
Dcq = sparse(II(1:k), JJ(1:k), 1);

% Remove duplicate cliques
done = 0;
while ~done
    % Form the clique intersection matrix M
    M = Dcq'*Dcq;
    % Determine the cliques jj that are a subset of another clique. The
    % entries on diag(M) correspond to the size of the cliques and the
    % entries on the off-diagonal correspond to the size of the
    % intersection of pairs of cliques. If the size of a clique equals the
    % size of the intersection with another clique, then it is a subset of
    % another clique.
    [~, jj] = find(triu(bsxfun(@eq, diag(M), M), 1));
    if isempty(jj)
        done = 1;
    else
        % Delete the duplicate cliques
        notjj = true(size(Dcq, 2), 1);
        notjj(jj) = 0;
        Dcq = Dcq(:, notjj);
    end
end

% Check that no clique is a subset of another clique
if nnz(bsxfun(@eq, diag(M), M) - speye(size(M)))
    error('There is a clique that is a subset of another clique\n')
end

% Form cliques cell array
numCliques = size(Dcq, 2);
cliques = cell(numCliques, 1);
for jj = 1:numCliques
    cliques{jj} = find(Dcq(:,jj));
end

end
