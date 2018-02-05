function [divhist, i2hist, inhist] = schneidman2c(spmat, dt, nperms, nICsuse, jtol, plotting)
%
% Plot ratio of I_2 to I_N for many independently chosen subnetworks. Also
% plot the histogram of the Jensen-Shannon divergences between the
% distribution of activity in the subnetwork and the corresponding model,
% P2 or P1.
%

nICs = size(spmat,2);

if nargin<6
    plotting = 0;
end
if nargin<4
    nICsuse = 5; % Number of neurons to include in each group
end
if nargin<5
    jtol = 1e-6;
end

if plotting
    handle = findobj('Tag','Schneidman2bc');
    if isempty(handle)
        figure('Tag','Schneidman2bc','Name','Fig. 2b, 2C of Schneidman et al.',...
            'NumberTitle','on');
    else
        figure(handle)
        clf
    end
end

divhist=[]; i2hist=[]; inhist=[]; i1hist=[];
icperms = nchoosek(nICs,10);
plotting = 1;
icord_sigpairs = zeros(nperms,nICsuse);
hvec = []; jmat = [];
for jperms = 1:nperms
%     for jperms = 101:1000
    icord_rand = randperm(nICs);
    icord_sigpairs(jperms,:) = icord_rand(1:nICsuse);
    [divhist(jperms,:), i2hist(jperms), inhist(jperms), i1hist(jperms), hvec, jmat] = isingjs(spmat(:,icord_sigpairs(jperms,:)), dt, jtol, 0, hvec, jmat);
%     if (plotting)&(jperms>2)&(mod(jperms,20)==1)
    if (plotting)&(jperms>2)
        figure(handle)
        divhistx = logspace(log10(min(divhist(:))), log10(max(divhist(:))), max(20, jperms/2));
        subplot(1,2,1)
        y = histc(divhist,divhistx);
        bar(log10(divhistx), y, 3)
%         hold on
%         bar(log10(divhistx), y(:,2), .9, 'b')
%         hold off
        set(gca,'FontSize',22,'FontWeight','bold','FontName','Helvetica','FontAngle','normal','LineWidth',2,'TickLength',[0.01,0.02])
        legend('D_{JS}(P_2,P)','D_{JS}(P_1,P)')
        set(gca,'XTick', [-5:0])
        xt = get(gca,'XTick');
        set(gca,'XTickLabel',num2str(10.^(xt')))
        set(gca,'XScale','log')
        xlabel('D_{JS} (bits)')
        ylabel('No. of subnetworks')

        subplot(1,2,2)
        plot(inhist/dt, i2hist./inhist, 'bo', 'MarkerFaceColor', 'b')
%         hold on
%         plot(inhist/dt, i1hist./inhist, 'ko', 'MarkerFaceColor', 'k')
%         hold off
        set(gca,'FontSize',22,'FontWeight','bold','FontName','Helvetica','FontAngle','normal','LineWidth',2,'TickLength',[0.01,0.02])
        xlabel('Full multi-information I_N (bits/s)')
        ylabel('I_{(2)}/I_N')
        set(gcf,'Color','w','PaperPositionMode','auto')
    end
end
