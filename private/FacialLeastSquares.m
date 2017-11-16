function [Z, eigZ] = FacialLeastSquares(D, U, opts)

r = size(U, 2);

% Form least-squares problem coefficient matrix A
rr = r*(r+1)/2;
[II, JJ, dd] = find(triu(D));

AT = zeros(rr, length(dd));
UT = U';
for k = 1:length(dd)
    ii = II(k);
    jj = JJ(k);
    UTeij = UT(:,ii) - UT(:,jj);
    M = UTeij*UTeij';
    
    %AT(:,k) = svec(UTeij*UTeij');
    p = 1;
    for jj=1:r
        for ii=1:jj-1
            AT(p,k) = sqrt(2)*M(ii,jj);
            p = p + 1;
        end
        AT(p,k) = M(jj,jj);
        p = p + 1;
    end
    
end
AA = AT';

% Solve least-squares problem
if opts.weighted_ls
    w = dd.^(-1);
    z = lscov(AA, dd, w);
else
    z = lscov(AA, dd);
end
Z = sMat(z);

eigZ = eig(Z)';

if min(eigZ)<0
    % Generate a 'sketch' matrix
    Sk = randn(2*size(AA,2),size(AA,1));
    Skdd = Sk*dd; 
    SkAA = Sk*AA; %#ok<*NASGU>
    
    % Solve the SDLS problem
    cvx_begin sdp
        cvx_quiet(true)
        cvx_precision best
        variable Z(r,r) symmetric
        minimize norm(SkAA*svec(Z) - Skdd)
        Z >= 0 %#ok<NOPRT>
    cvx_end
    
    eigZ = eig(Z)';
end

% if opts.verbose
%     res = norm(H.*(D-localK(U*Z*U')),'fro')/norm(H.*D,'fro');
%     fprintf('rel. res. norm(H.* K(U Z U'')-D)/norm(H.*D) %g \n',res)
% end

end    % of function
