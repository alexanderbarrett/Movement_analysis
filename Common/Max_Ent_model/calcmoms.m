function [sm,ssm,z] = calcmoms(hvec, jmat)
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

spatt2 = zeros(nICs, nICs, length(pattvec));
tic
for j=1:length(pattvec)
    spatt2(:,:,j) = spatt(j,:)' * spatt(j,:);
end
spatt2 = reshape(spatt2, nICs^2, length(pattvec))';

x = exp( spatt*hvec + spatt2*jmat/2 )'; % Relative probability of each pattern
z = sum(x); % Partition function

sm = (x*spatt) / z;
ssm = (x*spatt2) / z;
ssm = reshape(ssm, nICs, nICs);
z = z/npatts;
% fprintf('Finished calculating Ising moments; ')
% toc
