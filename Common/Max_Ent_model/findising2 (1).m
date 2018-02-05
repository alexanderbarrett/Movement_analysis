function [hvec,jmat] = findising(sexpt, cexpt, jtol, plotting, h0, j0)
%
% Use Newton-Raphson method to solve for parameters of the Ising
% model with the specified first and second order statistics.
%


if nargin<4
    plotting = 1;
end

nICs = length(sexpt);

% Initialize
if nargin<5
    % hvec = randn(nICs,1);
    hvec = -1+0.2*randn(nICs,1);
    jmat = 0.4*rand(nICs);
    jmat = jmat - diag(diag(jmat));
    jmat = (jmat+jmat')/2;
else
    hvec = h0;
    jmat = j0;
end

% Numerical search
alph = 0;
jchange = jtol+1;
eta = .95;
tic
while (jchange(end)>jtol)
    [sising,cising] = calcmoms(hvec,jmat);
    
    [f, jf] = isingobjective(hvec,jmat, sexpt, cexpt)
    
    xold = [hvec, jmat(c2use(:)==1)];
    
    xnew = xold - inv(jf)*
    
    %     [sising,cising] = calcmoms_mc(hvec,jmat);
    hvec = (1+alph)*hvec  - eta*(sising' - sexpt');
    jmat = (1+alph)*jmat  - eta*(cising - cexpt);
    %     size([((sising./sexpt)' -1);(cising(:)./cexpt(:)-1)])
    %     size([(sising./sexpt)'-1); (cising(:)./cexpt(:)-1) ])
    dcsvec = abs([cising(:)./cexpt(:)-1; (sising./sexpt)'-1]);
    jchange = [jchange, max(dcsvec(isfinite(dcsvec)))];
    if plotting && (mod(length(jchange),100)==2)
        %         eta = min(1,10000/length(jchange));
        figure(3)
        subplot(1,2,1)
        plot(cexpt, abs(cising./cexpt-1), 'r.')
        hold on
        plot(abs(sexpt), abs(sising./sexpt-1), 'b.')
        hold off
        set(gca,'YScale','log')
        set(gca,'FontSize',22,'FontWeight','bold','FontName','Helvetica','FontAngle','normal','LineWidth',2,'TickLength',[0.01,0.02])
        xlabel('Absolute value of correlation')
        ylabel('Fit distance')
        axis tight
        
        subplot(1,2,2)
        plot(jchange, '-')
        set(gca,'YScale','log','YLim',[jtol,1])
        set(gca,'FontSize',22,'FontWeight','bold','FontName','Helvetica','FontAngle','normal','LineWidth',2,'TickLength',[0.01,0.02])
        xlabel('Monte carlo step #')
        ylabel('Change in value of estimate')
        fprintf('Round %3.0f, change %3.6f. ', length(jchange), jchange(end))
        toc
        drawnow
    end
end
jchange = [jchange, abs(max([cising(:)./cexpt(:)-1; (sising./sexpt)'-1]))];

if plotting
    figure(3)
    subplot(1,2,1)
    plot(cexpt, abs(cising./cexpt-1), 'r.')
    hold on
    plot(abs(sexpt), abs(sising./sexpt-1), 'b.')
    hold off
    set(gca,'YScale','log')
    set(gca,'FontSize',22,'FontWeight','bold','FontName','Helvetica','FontAngle','normal','LineWidth',2,'TickLength',[0.01,0.02])
    xlabel('Absolute value of correlation')
    ylabel('Fit distance')
    axis tight

    subplot(1,2,2)
    plot(jchange, '-')
    set(gca,'YScale','log','YLim',[jtol,1])
    set(gca,'FontSize',22,'FontWeight','bold','FontName','Helvetica','FontAngle','normal','LineWidth',2,'TickLength',[0.01,0.02])
    xlabel('Monte carlo step #')
    ylabel('Change in value of estimate')
    drawnow
end
fprintf('Ising fit converged in %3.0f steps.\n', length(jchange))
