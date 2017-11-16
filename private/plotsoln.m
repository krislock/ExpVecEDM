function plotsoln(X, Xref, Xorig, A)

m = size(A, 1);
n = size(X, 1) + m;

% Plot solution
clf
if ~isempty(Xref)
    subplot(1, 2, 1);
    plotXsol(X, 'Before refinement');
else
    plotXsol(X, '');
end

if ~isempty(Xref)
    subplot(1, 2, 2);
    plotXsol(Xref, 'After refinement');
end

%%%%%%%% plotXsol subroutine %%%%%%%%
    function plotXsol(Xsol, titletext)
        plot(Xorig(:,1), Xorig(:,2), 'ko')
        hold on
        for ii=1:n-m
            plot([Xsol(ii,1) Xorig(ii,1)], ...
                [Xsol(ii,2) Xorig(ii,2)], 'k-');
        end
        plot(A(:,1), A(:,2), 'rs')
        plot(Xsol(:,1), Xsol(:,2), 'k.')
        plot(A(:,1), A(:,2), 'r.')
        axis('equal')
        axis([-0.6, 0.6, -0.6, 0.6])
        title(titletext)
        hold off
    end

end

