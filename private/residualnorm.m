function res = residualnorm(P, Dpartial)

[II, JJ, dd] = find(triu(Dpartial));

numEdges = length(dd);

residual = zeros(size(dd));
for k = 1:numEdges
    residual(k) = norm(P(II(k),:) - P(JJ(k),:))^2 - dd(k);
end
res = norm(residual)/sqrt(numEdges);

end
