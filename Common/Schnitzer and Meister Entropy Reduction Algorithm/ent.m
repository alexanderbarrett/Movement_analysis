function entout = ent(p)
% Compute entropy of a binary variable with rate p
% entout - entropy (bits)

pzero = logical((p<=0) + (p>=1));
p(pzero) = 0.5;
entout = -(p.*log2(p) + (1-p).*log2(1-p));
entout(pzero) = 0;
