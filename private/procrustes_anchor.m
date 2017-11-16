function Xout = procrustes_anchor(A, P)

m = size(A, 1);
n = size(P, 1);

xinds = 1:n-m;
ainds = n-m+1:n;

% Translate P so that P2 is centred at zero
P2 = P(ainds,:);
mp2 = sum(P2)/m;  % centre of P2
Pc = P - repmat(mp2,n,1);  % translate P
P2c = Pc(ainds,:);

% Translate A so that A is centred at zero
ma = sum(A)/m;  % centre of A
Ac = A - repmat(ma,m,1);  % centre A at zero

% solve min ||A-P2*Q|| s.t. Q'*Q = I
% this is the orthogonal Procrustes problem
%   -- Golub, Van Loan, Algorithm 12.4.1
C = P2c'*Ac;
[Uc,~,Vc] = svd(C);
Q = Uc*Vc';

% Compute final X
Pc = Pc*Q;  % rotate P to align P2 with anchors A

Pout = Pc + repmat(ma,n,1);  % translate P back

Xout = Pout(xinds,:);

end