function D = K(B)
% D = K(B)
%
% Operator K of B  - forms EDM if B is psd

D = zeros(size(B));
for ii=1:length(B)
    for jj=1:length(B)
        D(ii,jj) = B(ii,ii) + B(jj,jj) - 2*B(ii,jj);
    end
end

end