function B = Kdag(D)
% B = Kdag(D)
%
% Pseudoinverse of K - forms centred Gram matrix B if D is an EDM

nn = length(D);

D = D/2;
vn = sum(D,2)/nn;
evnn = sum(vn)/nn;

B = zeros(nn);

for ii=1:nn
    for jj=1:nn
        B(ii,jj) = vn(ii) + vn(jj);
    end
end

B = B - D;
B = B - evnn;

end