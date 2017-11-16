function [U, eigY] = face(Y, r)

state = rng;  % Get the rng state

n = length(Y);
sigma = (trace(Y) + 1)/n;

if n > 1000   % use eigs
    I = speye(n);
    
    % f(x) returns (Y + I + sigma*E)\x
    [Ry, ~, S] = chol(Y + I);
    alpha = sigma/(1 + sigma*n)*ones(n, 1);
    f = @(x) S*(Ry\(Ry'\(S'*x))) - alpha*sum(x);
    
    % Compute the smallest r eigenvalues and their eigenvectors
    eigsopts.issym = true;
    eigsopts.maxit = 20;
    eigsopts.disp = 0;
    
    k = r;
    flag = 1;
    while flag
        rng(0);  % Make eigs deterministic
        [Uy, dy, flag] = eigs(f, n, k, 'SM', eigsopts);  dy = diag(dy);
        if flag
            fprintf('eigs failed with k = %d\n', k);
        end
        k = min([10*k, n]);  % Try a larger value of k
        fprintf('trying eigs with larger k = %d\n', k);
    end
    
    % Compute face
    eigY = dy(end-r+1:end)' - 1;
    U = Uy(:,end-r+1:end);
    
else  % use eig
    [Uy, dy] = eig(full(Y)+sigma*ones(n)); % deflate ones eigenvector
    dy = diag(dy);
    eigY = dy(1:r)';
    U = Uy(:,1:r);
    
end

rng(state);  % Restore rng state

end