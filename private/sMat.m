function M = sMat(v)

p = length(v);
n = (sqrt(1 + 8*p) - 1)/2;

k = 1;
M = zeros(n);
for jj=1:n
    for ii=1:jj-1
        M(ii,jj) = v(k)/sqrt(2);
        M(jj,ii) = M(ii,jj);
        k = k + 1;
    end
    M(jj,jj) = v(k);
    k = k + 1;
end

end