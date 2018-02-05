function [xc]=contributions(A,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the true contributions
%
% Inputs
% A : mixing matrix whose A(i,j) element is the scale coefficient 
% associated to the j^th source and to the i^th observation
% s : source vector
%
% Output
% xc : 3D matrix whose xc(:,i,j) element is the true contribution of the
% j^th source in the i^th observation

xc=zeros(length(s(1,:)),size(A,1),size(A,2));
for i=1:size(A,1)
    for j=1:size(A,2)
        xc(:,i,j)=filter(A(i,j),1,s(j,:)');
    end
end