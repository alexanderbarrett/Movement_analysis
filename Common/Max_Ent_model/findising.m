function [hvec,jmat] = findising(sexpt, cexpt, jtol, plotting)
%
% Use (approximate) gradient descent to solve for parameters of the Ising
% model with the specified first and second order statistics.
%

if nargin<4
    plotting = 1;
end

nICs = length(sexpt);

% hvec = randn(nICs,1);
hvec = -1+0.2*randn(nICs,1);
jmat = 0.4*rand(nICs);
jmat = jmat - diag(diag(jmat));
jmat = (jmat+jmat')/2;

% Numerical search
alph = 0;
jchange = jtol+1;
eta = .85;
tic
while (jchange(end)>jtol)
    [sising,cising] = calcmoms(hvec,jmat);
    %     [sising,cising] = calcmoms_mc(hvec,jmat);
    hvec = (1+alph)*hvec  - eta*(sising' - sexpt');
    jmat = (1+alph)*jmat  - eta*(cising - cexpt);
    %     size([((sising./sexpt)' -1);(cising(:)./cexpt(:)-1)])
    %     size([(sising./sexpt)'-1); (cising(:)./cexpt(:)-1) ])
    dcsvec = abs([cising(:)./cexpt(:)-1; (sising./sexpt)'-1]);
    jchange = [jchange, max(dcsvec(isfinite(dcsvec)))];
    if plotting && (mod(length(jchange),10000)==2)
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
        ylabel('Max. difference b/w expt and model')
        fprintf('Finished %4.0f MC steps, current criterion = %4.4g; ', length(jchange), jchange(end))
        toc
        tic
        drawnow
    end
%             fprintf('Round %3.0f, change %3.6f.\n', length(jchange), jchange(end))
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
