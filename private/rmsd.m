function res = rmsd(Xorig, X)

res = norm(X - Xorig, 'fro')/sqrt(size(Xorig,1));

end