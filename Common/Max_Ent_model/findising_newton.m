function [hvec, jmat] = findising_newton(spmat, errtol, plotting, h0, j0)
%
% Calculate moments of Ising model with local fields given by hvec,
% interactions given by jmat.
%
% Note that sm, ssm are defined in the interval [-1,1]
%

if (nargin<2)|isempty(errtol)
    errtol = 1e-5;
end
if nargin<3
    plotting = 1;
end

if plotting
    handle = findobj('Tag','findising');
    if isempty(handle)
        handle = figure('Tag','findising','Name','Progress of Ising model finding',...
            'NumberTitle','on');
    else
        figure(handle)
        clf
    end
end

nICs = size(spmat,2);
nsamp = size(spmat,1);
spmat = 2*spmat - 1;
cexpt = spmat'*spmat/size(spmat,1);
sexpt = mean(spmat,1);

nICs = length(sexpt);
pattvec = [0:2^nICs-1];
spatt = (dec2bin(pattvec, nICs)=='1')*2 - 1;
c2use = triu(ones(nICs),1);
cexpt = cexpt(c2use(:)==1);

% tic
spatt2full = zeros(length(pattvec), nICs, nICs);
spatt2 = zeros(length(pattvec), nchoosek(nICs, 2));
spatt3 = zeros(length(pattvec), nchoosek(nICs, 2), nICs);
spatt4 = zeros(length(pattvec), nchoosek(nICs, 2), nchoosek(nICs, 2));
for j=1:length(pattvec)
    spatt2curr = spatt(j,:)' * spatt(j,:);
    spatt2full(j,:,:) = spatt2curr;
    spatt2(j,:) = spatt2curr(c2use(:)==1);
    spatt3(j,:,:) = spatt2(j,:)' * spatt(j,:);
    spatt4(j,:,:) = spatt2(j,:)' * spatt2(j,:);
end
% toc

% Initialize
if (nargin<4)|isempty(h0)
    % hvec = randn(nICs,1);
    hvec = -1+0.2*randn(nICs,1);
    jmat = 0.4*rand(nICs);
    jmat = jmat - diag(diag(jmat));
    jmat = (jmat+jmat')/2;
    jmat = jmat(c2use==1);
else
    hvec = h0;
    jmat = j0;
end

tic
k=0;
maxdparams = 1; maxfvec = Inf;
while (maxfvec > errtol)
    k=k+1;

    % Compute correlations up to 4th order
    x = exp( spatt*hvec + spatt2*jmat )'; % Relative probability of each pattern
    z = sum(x); % Partition function

    sm = (x*spatt) / z;
    c2m = (x*spatt2) / z;
    c2mfull = squeeze(sum(repmat(x',[1,nICs,nICs]).*spatt2full,1) / z);
    c3m = squeeze(sum(repmat(x',[1,nchoosek(nICs, 2),nICs]).*spatt3,1) / z);
    c4m = squeeze(sum(repmat(x',[1,nchoosek(nICs, 2),nchoosek(nICs,2)]).*spatt4,1) / z);

    % Objective function:
    fvec = [(sm - sexpt),(c2m - cexpt')];
    f = sum(fvec.^2);

    % Jacobian:
    %     jf = diag(fvec) * [ c2mfull - sm'*sm, c3m' - sm'*c2m; c3m - c2m'*sm, c4m - c2m'*c2m ];
    jf = [ c2mfull - sm'*sm, c3m' - sm'*c2m; c3m - c2m'*sm, c4m - c2m'*c2m ];

    params = [ hvec; jmat ];

    dparams = jf \ fvec';
    maxfvec = max(fvec(:));
    %     params = params - dparams/(max([3,.2*sqrt(k)]) * max([1;abs(dparams(:))]));
    params = params - dparams/(max([1;abs(dparams(:))]));
    hvec = params(1:nICs);
    jmat = params(nICs+1:end);

    if plotting
        figure(handle)
        mincorr = min(abs([cexpt(:); sexpt(:)]));
        plot( abs(sexpt), abs(sm), 'ro', abs(cexpt), abs(c2m), 'bo', 'LineWidth', 2)
        hold on
        plot( [mincorr,1], [mincorr,1], 'k--', 'LineWidth',2)
        hold off
        set(gca,'FontSize',22,'FontWeight','bold','FontName','Helvetica','FontAngle','normal','LineWidth',2,'TickLength',[0.01,0.02])
        xlabel('Observed correlation')
        ylabel('Model correlation')
    end
    drawnow
%     fprintf('Step %3.0f, error %3.0g; ', k, maxfvec)
%     toc
end
fprintf('Model finding converged in %3.0f steps; ', k)
toc