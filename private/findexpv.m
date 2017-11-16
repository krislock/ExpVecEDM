function [U, infeasdist] = findexpv(Dclique, r)
% [U, infeasdist] = findexpv(Dpartial, r)

B = Kdag(Dclique);
k = length(Dclique);

% we deflate the ones eigenvector
lamr = norm(B, 'fro');  % this is a bound on the spectral radius
%[Ub,db] = eig(B - lamr*ones(k)/k, 'vector');
[Ub,db] = eig(B - lamr*ones(k)/k);  db = diag(db);

U = Ub(:,2:k-r);
infeasdist = (norm(db(2:k-r))^2 + ...
    norm(min(0,db(k-r+1:k)))^2)/(0.5*k*(k-1));

end