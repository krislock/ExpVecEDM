function Y = ExposingVector(Dpartial, cliques, r, opts)
% Y = expvec(Dpartial, cliques, r, opts)

n = length(Dpartial);
numCliques = length(cliques);

% Compute exposing vectors
csizes = cellfun(@length, cliques);
K = sum(csizes.*(csizes + 1)/2);
II = zeros(K, 1);
JJ = zeros(K, 1);
vv = zeros(K, 1);
k = 0;  % Index for II, JJ, vv vectors

if opts.expvecweights
    expvs = cell(numCliques, 1);
    infeasdist = zeros(numCliques, 1);
    for kk = 1:numCliques
        clq = cliques{kk};
        Dclique = full(Dpartial(clq, clq));
        [expvs{kk}, infeasdist(kk)] = findexpv(Dclique, r);
    end
    
    % Compute weights
    expwts = 1 - infeasdist/sum(infeasdist);
    
    % Form final exposing vector Y
    [~, sortinds] = sort(expwts, 'ascend');
    for kk = sortinds'
        clq = cliques{kk};
        U = expvs{kk};
        wt = expwts(kk);
        B = wt*(U*U');
        for jj = 1:length(clq)
            for ii = 1:jj
                k = k + 1;
                II(k) = clq(ii);
                JJ(k) = clq(jj);
                vv(k) = B(ii,jj);
            end
        end
    end
    
else
    % Form final exposing vector Y
    for kk = 1:numCliques
        clq = cliques{kk};
        Dclique = full(Dpartial(clq, clq));
        U = findexpv(Dclique, r);
        B = U*U';
        for jj = 1:length(clq)
            for ii = 1:jj
                k = k + 1;
                II(k) = clq(ii);
                JJ(k) = clq(jj);
                vv(k) = B(ii,jj);
            end
        end
    end
    
end

Y = sparse(II, JJ, vv, n, n);
Y = Y + triu(Y, 1)';

end
