function [constitcells, totent, entthresh] = entredcorr_clustplot(spmat, entthresh, idx, maxrounds)
% Schnitzer & Meister entropy-reduction algorithm
%
% constitcells - matrix containing indices of cell1, cell2, \DeltaH
%

% if nargin<2
%     plotting = 0;
% end

plotting = 1; % Current version requires plotting

nsamp = size(spmat,1);
nICs = size(spmat,2);
fprintf('Analyzing %0.0f spike trains, %0.0f time bins...\n', nICs, nsamp)

nICsorig = nICs;
constitcells(:,1) = [1:nICs];
constitcells(:,2) = 0;

[spt,spc] = find(spmat);
[nspikes, spcu] = hist(spc, [1:size(spmat,2)]);
sprate = nspikes/nsamp;  % Probability of a spike per bin for each cell

spmat = logical( spmat );
spmatred = zeros(size(spmat,1), nICsorig);

% Determine entropy threshold, taking into account the number of cells
if (nargin<2)|(isempty(entthresh))
    fprintf('Finding entropy threshold using shuffled data...\n')
    % Shuffled data:
    shuffmin = 100; % Minimum number of time-bins for shuffle;
    nshuff=[500, min(1000,nsamp)]; % [Min,Max] number of shuffles
    entred = zeros(nICs*(nICs-1)/2, nshuff(2));

%     threshtol = min(5e-3, nICs*(nICs-1)/2) % Fraction of shuffled cell pairs above tolerance
    threshtol = 0;
    
    convergeFlag=0;
    jshuff=0;
    [jmat,kmat] = meshgrid([1:nICs]);
    spratej = sprate(jmat);
    spratek = sprate(kmat);
    entspratej = ent(spratej);
    entspratek = ent(spratek);
    jgreaterk = jmat>kmat;
    tshuff = [];
    for j=1:nICsorig
        tshuff = [tshuff,randperm(nsamp)]; % Choosing shuffles this way ensures no two shuffles will be identical
    end
    tic
    while ((~convergeFlag)||(jshuff<nshuff(1)))&&(jshuff<nshuff(2))
        jshuff=jshuff+1;
        %     tstart = round(rand(nICsorig,1)*(nsamp - 2*shuffmin))+1;
        for j=1:nICsorig
            spmatred(:,j) = circshift( spmat(:,j), tshuff(j+nICsorig*(jshuff-1)));
        end

        % Compute change in entropy matrix
        sprate2 = (spmatred' * spmatred)/nsamp;
        entredmat = entspratej + entspratek - ...
            (ent(spratej - sprate2) + ent(spratek - sprate2) + ent(sprate2));
        entred(:,jshuff) = entredmat(jgreaterk);

        % Plot
        if plotting&&(mod(jshuff, 200)==0)
            [entthresh, convergeFlag] = plothist(entred(:,1:jshuff), threshtol);
            %         fprintf('Finished %4.0f shuffles; tolerance is %3.5g; ', jshuff, centh(1))
            fprintf('Finished %4.0f shuffles;  ', jshuff)
            toc
        end
    end
    fprintf('Finished %4.0f shuffles;  ', jshuff)
    toc
    [enth,entx]=hist(entred(:),100);
    dx = entx(2)-entx(1);
    if plotting
        entthresh = plothist(entred(:,1:jshuff), threshtol);
        fprintf('DeltaH thresh = %3.6f.\n', entthresh)
    end
end


if entthresh<0
    fprintf('Error DeltaH threshold is negative; setting to zero.\n')
    entthresh = 0;
end

% Original data:
spmatred = double(spmat);

sprate = sum(spmatred,1)/nsamp;
l=0;
nICs = size(spmat,2);
nICsorig = nICs;
maxentred(1) = 1;
totent = sum(ent(sprate))*ones(nICs,1);
maxentred = zeros(nICs,1);
maxentred(end) = 1;

% Compute initial entropy matrix
sprate2 = (spmatred' * spmatred) / nsamp;
[spratej,spratek] = meshgrid(sprate);
entred = ent(spratej) + ent(spratek) - ...
    (ent(spratej - sprate2) + ent(spratek - sprate2) + ent(sprate2));
entred = triu(entred, 1) - tril(ones(size(entred)),0);
totent(nICs) = sum(ent(sprate));

while (maxentred(end)>entthresh)&(l<maxrounds)
    l=l+1;    
    %     for j=1:nICs
    %         totent(nICs) = totent(nICs) + ent(sprate(j));
    %         for k=j+1:nICs
    %             sprate2(j,k) = sum(spmatred(:,j).*spmatred(:,k))/nsamp;
    %             entred(j,k) = (ent(sprate(j))+ent(sprate(k))) - ...
    %                 (ent(sprate(j) - sprate2(j,k))+ent(sprate(k) - sprate2(j,k))+ent(sprate2(j,k)));
    %         end
    %     end

    maxentred(nICs+1) = max(entred(:));
    if maxentred(nICs+1)>entthresh
        nICs = nICs+1;
        [x1,x2] = find(entred==maxentred(nICs), 1);
        constitcells(nICs,1) = x1;
        constitcells(nICs,2) = x2;
        constitcells(nICs,3) = maxentred(nICs);
%         fprintf('Combined cells %3.0f, %3.0f; DeltaH = %3.8f.\n', constitcells(nICs,1), constitcells(nICs,2), maxentred(nICs))
        spmatred(:,nICs) = spmatred(:,constitcells(nICs,1)).*spmatred(:,constitcells(nICs,2));
        spmatred(:,constitcells(nICs,1)) = spmatred(:,constitcells(nICs,1))-spmatred(:,nICs);
        spmatred(:,constitcells(nICs,2)) = spmatred(:,constitcells(nICs,2))-spmatred(:,nICs);
        sprate = sum(spmatred,1)/nsamp;
        totent(nICs) = sum(ent(sprate));
        
        % Update entropy reduction matrix
        sprate2 = (spmatred' * spmatred) / nsamp;
        entred = padarray(entred, [1,1], -1, 'post');
        entred(nICs,:) = ent(sprate(nICs)) + ent(sprate) - ...
            (ent(sprate(nICs) - sprate2(nICs,:))+ent(sprate - sprate2(nICs,:))+ent(sprate2(nICs,:)));
        entred(x1,:) = ent(sprate(x1)) + ent(sprate) - ...
            (ent(sprate(x1) - sprate2(x1,:)) + ent(sprate - sprate2(x1,:)) + ent(sprate2(x1,:)));
        entred(x2,:) = ent(sprate(x2)) + ent(sprate) - ...
            (ent(sprate(x2) - sprate2(x2,:)) + ent(sprate - sprate2(x2,:)) + ent(sprate2(x2,:)));
        entred(:,nICs) = entred(nICs,:);
        entred(:,x1) = entred(x1,:);
        entred(:,x2) = entred(x2,:);
        entred = triu(entred, 1) - tril(ones(size(entred)),0);
        
        if (plotting)&(mod(l,100)==0)
%             [spt, spc] = find(spmatred);
%             subplot(2,1,2)
%             scatter(spt(spc<=nICsorig), spc(spc<=nICsorig), 'b.')
%             hold on
%             scatter(spt(spc>nICsorig), spc(spc>nICsorig), 'ro')
%             hold off
%             set(gca,'FontSize',22,'FontWeight','bold','FontName','Helvetica','FontAngle','i','LineWidth',2,'TickLength',[0.01,0.02])
%             xlabel('Time bin')
%             ylabel('Cell group #')
%             box on
%             axis tight
%             drawnow
            
entplot_zones(constitcells, totent, entthresh, idx);

        end
    end
end

% [spt, spc] = find(spmatred);
% subplot(2,1,2)
% scatter(spt(spc<=nICsorig), spc(spc<=nICsorig), 'b.')
% hold on
% scatter(spt(spc>nICsorig), spc(spc>nICsorig), 'ro')
% hold off
% set(gca,'FontSize',22,'FontWeight','bold','FontName','Helvetica','FontAngle','i','LineWidth',2,'TickLength',[0.01,0.02])
% xlabel('Time bin')
% ylabel('Cell group #')
% box on
% axis tight
% ylim([0,max(spc(:))+1])
% drawnow

function [entthresh, convergeFlag] = plothist(entred, threshtol)
% Plot histogram of the density of DeltaH, return the entropy threshold

nICs = size(entred,1);
jshuff = size(entred,2);
entred = entred(:);
entx = linspace(min(entred(:))*1.1, max(entred(:))+0.1*abs(max(entred(:))), 50);
enth = hist(reshape(entred, 1, []), entx);
dx = entx(2)-entx(1);
centx = sort(entred, 'ascend');
centh = [nICs*jshuff:-1:1]/jshuff;
% centh = cumsum( flipud(hist(entred, centx)) ) / jshuff;
% centx = fliplr(centx);
% centx = [max(entred(:)), centx];
% centh = [0, centh];
[centh, i] = unique(centh);
centx = centx(i);

subplot(2,2,1)
plot(entx, enth/(dx*length(entred(:))),'b','LineWidth',2)
% set(gca,'YScale','log')
set(gca,'FontSize',22,'FontWeight','bold','FontName','Helvetica','FontAngle','i','LineWidth',2,'TickLength',[0.01,0.02])
xlabel({'\DeltaH_{shuffle} (bits)'})
ylabel('Density (1/bits)')
axis tight

subplot(2,2,2)
plot(centx, centh, '-','LineWidth',2)
hold on
plot(entx, ones(size(entx)) * threshtol, 'r--','LineWidth',2)
if centh(1)<=threshtol
    entthresh = interp1(centh, centx, threshtol, 'linear');
%     fprintf('Using entropy threshold DeltaH = %3.3g.\n', entthresh)
    convergeFlag = 1;
else
    entthresh = max(entred(:));
%     fprintf('Shuffle statistics failed; using maximum entropy threshold DeltaH = %3.3g.\n', entthresh)
    convergeFlag = 0;
end
plot([1,1]*entthresh, [threshtol/10, abs(centh(end))*1.1], 'r--','LineWidth',2)
hold off
set(gca,'FontSize',22,'FontWeight','bold','FontName','Helvetica','FontAngle','i','LineWidth',2,'TickLength',[0.01,0.02])
set(gca,'YTick',[1e-4,1e-2,1e0,1e2])
ylabel({'Avg. # of cell','pairs > \DeltaH'})
xlabel({'\DeltaH_{shuffle} (bits)'})
set(gca,'YScale','log')
axis tight
set(gcf,'Color','w','PaperPositionMode','auto')
drawnow