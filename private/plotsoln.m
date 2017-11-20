function plotsoln(X, Xref, Xorig, A)

m = size(A, 1);
n = size(X, 1) + m;
r = size(X, 2);

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
        if r==2
            plot(Xorig(:,1), Xorig(:,2), 'ko')
        elseif r==3
            plot3(Xorig(:,1), Xorig(:,2), Xorig(:,3), 'ko')
        end 
        hold on
        for ii=1:n-m
            if r==2
                plot([Xsol(ii,1) Xorig(ii,1)], ...
                    [Xsol(ii,2) Xorig(ii,2)], 'k-');
            elseif r==3
                plot3([Xsol(ii,1) Xorig(ii,1)], ...
                    [Xsol(ii,2) Xorig(ii,2)], ...
                    [Xsol(ii,3) Xorig(ii,3)], 'k-');
            end
        end
        if ~isempty(A)
            if r==2
                plot(A(:,1), A(:,2), 'rs')
                plot(A(:,1), A(:,2), 'r.')
            elseif r==3
                plot3(A(:,1), A(:,2), A(:,3), 'rs')
                plot3(A(:,1), A(:,2), A(:,3), 'r.')
            end
        end
        if r==2
            plot(Xsol(:,1), Xsol(:,2), 'k.')
        elseif r==3
            plot3(Xsol(:,1), Xsol(:,2), Xsol(:,3), 'k.')
        end
        title(titletext)
        hold off
    end

end
