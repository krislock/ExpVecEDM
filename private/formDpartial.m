function Dpartial = formDpartial(Xorig, A, R, nf)

Porig = [Xorig; A];

n = size(Porig, 1);
m = size(A, 1);

% This memory saving technique for computing the partial distance matrix
% is due to Franz Rendl of the University of Klagenfurt.
II = cell(n,1);
JJ = cell(n,1);
DD = cell(n,1);
for ip = 1:n-m
    I = false(n,1);
    Pdiff = Porig(ip+1:end,:) - repmat(Porig(ip,:),n-ip,1);
    I(ip+1:end) = all(abs(Pdiff) < R,2);
    Pdiff = Porig(I,:) - repmat(Porig(ip,:),nnz(I),1);
    Pdist = sum(Pdiff.^2,2);
    small_dist = Pdist < R^2;
    Pdist = Pdist(small_dist);
    I(I) = small_dist;
    II{ip} = ip*ones(sum(I),1);
    JJ{ip} = find(I);
    DD{ip} = Pdist;
end

Ainds = n-m+1:n;
AJJ = repmat(Ainds,m,1);
AII = AJJ';
ADD = triu(K(A*A'));

II = [cell2mat(II); AII(:)];
JJ = [cell2mat(JJ); AJJ(:)];

DDnoisy = cell2mat(DD);

noise = nf*randn(size(DDnoisy));
while any(abs(noise) > 1)
    inds = abs(noise) > 1;
    noise(inds) = nf*randn(nnz(inds), 1);
end

DDnoisy = DDnoisy.*(1 + noise).^2;
DD = [DDnoisy; ADD(:)];

% Form Dpartial
Dpartial = sparse(II, JJ, DD, n, n);
Dpartial = Dpartial + Dpartial';

end
