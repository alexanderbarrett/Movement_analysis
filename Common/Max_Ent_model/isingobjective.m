function [f, jf] = isingobjective(hvec, jmat, sexpt, cexpt)
%
% Calculate moments of Ising model with local fields given by hvec,
% interactions given by jmat.
%
% Note that sm, ssm are defined in the interval [-1,1]
%

% tic
jmat = jmat(:);
nICs = length(hvec);
npatts = 2^nICs;

pattvec = [0:2^nICs-1];
spatt = (dec2bin(pattvec, nICs)=='1')*2 - 1;

c2use = triu(ones(nICs),1);
cexpt = cexpt(c2use(:)==1);

% [i3,j3,k3] = ndgrid([1:nICs]);
% c3use = (k3>j3).*(j3>i3);
% [i4,j4,k4,l4] = ndgrid([1:nICs]);
% c4use = (l4>k4).*(k4>j4).*(j4>i4);
% clear i3 j3 k3 i4 j4 k4 c4use c3use

tic
spatt2full = zeros(length(pattvec), nICs, nICs);
spatt2 = zeros(length(pattvec), nchoosek(nICs, 2));
spatt3 = zeros(length(pattvec), nchoosek(nICs, 2), nICs);
spatt4 = zeros(length(pattvec), nchoosek(nICs, 2), nchoosek(nICs, 2));
% spatt3 = zeros(length(pattvec), nchoosek(nICs, 3));
% spatt4 = zeros(length(pattvec), nchoosek(nICs, 4));
for j=1:length(pattvec)
    spatt2curr = spatt(j,:)' * spatt(j,:);
    spatt2full(j,:,:) = spatt2curr;
    spatt2(j,:) = spatt2curr(c2use(:)==1);
    spatt3(j,:,:) = spatt2(j,:)' * spatt(j,:);
    %     spatt3(j,:) = spatt3curr(c3use==1);
    spatt4(j,:,:) = spatt2(j,:)' * spatt2(j,:);
    %     spatt4(j,:) = spatt4curr(c4use==1);
end
toc

hvec = -1+0.2*randn(nICs,1);
jmat = 0.4*rand(nICs);
jmat = jmat - diag(diag(jmat));
jmat = (jmat+jmat')/2;
jmat = jmat(c2use==1);

tic
k=0;
maxdparams = 1;
while (maxdparams>1e-4)
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
    maxdparams = mean(abs(dparams(:)));
    maxfvec = max(fvec(:));
    params = params - dparams/(max([3,.2*sqrt(k)]) * max(1,maxdparams));
    hvec = params(1:nICs);
    jmat = params(nICs+1:end);
    fprintf('Step %3.0f, error %3.0g; ', k, maxfvec)
    
    plot( abs(sexpt), abs(sm), 'r.', abs(cexpt), abs(c2m), 'b.', [.85,1],[.85,1], 'k--')
    drawnow
    toc
end
